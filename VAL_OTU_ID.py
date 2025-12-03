#!/usr/bin/env python3
"""
VAL_OTU_ID.py - Unified eDNA Taxonomy Validation Pipeline

Validates species identifications from eDNA/metabarcoding by:
1. Checking geographic plausibility against GBIF/OBIS occurrence databases
2. Validating database coverage against MIDORI2 reference sequences
3. Combining results into final taxonomy decisions

Author: Glenn Dunshea / Claude
Version: 1.0
"""

import argparse
import asyncio
import aiohttp
import pandas as pd
import subprocess
import sys
import tempfile
import re
from pathlib import Path
from collections import defaultdict
from math import radians, sin, cos, sqrt, atan2
from tqdm import tqdm

# =============================================================================
# CONFIGURATION - Adjust these parameters as needed
# =============================================================================

CONFIG = {
    # --- Geographic Validation Parameters ---
    'buffer_km': 500,              # Buffer around study area for "in range" (km)
    'distance_threshold_km': 2000, # Max distance for PLAUSIBLE/POSSIBLE categories
    'congener_search_buffer_km': 2000,  # Search area expansion for finding congeners
    
    # --- Sequence Identity Thresholds ---
    'min_species_pct': 97.0,       # Minimum % identity for confident species-level ID
    'min_genus_pct': 90.0,         # Minimum % identity for genus-level ID
    'min_family_pct': 80.0,        # Below this, drop to order
    'cf_threshold_pct': 97.0,      # Minimum % for "sp. cf." notation
    'reassign_diff_pct': 0.5,      # Minimum % difference to reassign with multiple congeners
    
    # --- API and Performance Settings ---
    'max_concurrent_requests': 10, # Concurrent API requests
    'worms_rate_limit': 0.1,       # Seconds between WoRMS requests
    
    # --- Sequence Comparison ---
    'kmer_size': 15,               # K-mer size for sequence comparison fallback
    'vsearch_path': 'vsearch',     # Path to vsearch executable (if installed)
}

# =============================================================================
# PHYLOSEQ EXTRACTION (calls R)
# =============================================================================

def extract_from_phyloseq(rds_path, output_dir):
    """
    Extract taxonomy and sequences from a phyloseq RDS file by calling R.
    
    Returns tuple of (taxonomy_csv_path, sequences_fasta_path)
    """
    rds_path = Path(rds_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Generate output filenames based on input
    base_name = rds_path.stem
    taxonomy_csv = output_dir / f"{base_name}_taxonomy.csv"
    sequences_fasta = output_dir / f"{base_name}_sequences.fasta"
    
    # R script to extract data
    r_script = f'''
    suppressPackageStartupMessages({{
      library(phyloseq)
      library(Biostrings)
    }})
    
    ps <- readRDS("{rds_path}")
    
    cat("  Loaded phyloseq:", ntaxa(ps), "ASVs,", nsamples(ps), "samples\\n")
    
    # Extract taxonomy
    tax <- as.data.frame(tax_table(ps))
    
    # Create ASV IDs
    if (all(grepl("^[ACGT]+$", rownames(tax)))) {{
      tax$ASV_ID <- paste0("ASV", seq_len(nrow(tax)))
    }} else {{
      tax$ASV_ID <- rownames(tax)
    }}
    
    # Move ASV_ID to first column
    tax <- tax[, c("ASV_ID", setdiff(names(tax), "ASV_ID"))]
    write.csv(tax, "{taxonomy_csv}", row.names = FALSE)
    
    # Extract sequences
    seqs <- refseq(ps)
    if (is.null(seqs)) {{
      if (all(grepl("^[ACGT]+$", rownames(tax_table(ps))))) {{
        seqs <- DNAStringSet(rownames(tax_table(ps)))
        names(seqs) <- tax$ASV_ID
      }} else {{
        stop("No sequences found in phyloseq object")
      }}
    }} else {{
      if (all(grepl("^[ACGT]+$", names(seqs)))) {{
        names(seqs) <- tax$ASV_ID
      }}
    }}
    writeXStringSet(seqs, "{sequences_fasta}")
    
    cat("  Extracted taxonomy and sequences\\n")
    '''
    
    # Run R script
    try:
        result = subprocess.run(
            ['Rscript', '-e', r_script],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            print(f"  R error: {result.stderr}")
            raise RuntimeError(f"Failed to extract from phyloseq: {result.stderr}")
        
        # Print R output
        if result.stdout.strip():
            print(result.stdout.strip())
        
    except FileNotFoundError:
        raise RuntimeError("Rscript not found. Please ensure R is installed and in your PATH.")
    except subprocess.TimeoutExpired:
        raise RuntimeError("R extraction timed out. The phyloseq object may be too large.")
    
    # Verify outputs exist
    if not taxonomy_csv.exists():
        raise RuntimeError(f"Failed to create taxonomy file: {taxonomy_csv}")
    if not sequences_fasta.exists():
        raise RuntimeError(f"Failed to create sequences file: {sequences_fasta}")
    
    return taxonomy_csv, sequences_fasta


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculate distance between two points in km."""
    R = 6371
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    return R * c


def clean_species_name(name):
    """Clean species name for querying."""
    if pd.isna(name) or not name:
        return None
    name = str(name).strip()
    name = re.sub(r'\s+', ' ', name)
    name = re.sub(r'\bcf\.?\s*', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\baff\.?\s*', '', name, flags=re.IGNORECASE)
    name = re.sub(r'(?<![a-zA-Z])(sp\.|spp\.?)\s*$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\s*\(.*?\)\s*', ' ', name)
    name = re.sub(r'\s+', ' ', name).strip()
    if not name or len(name) < 3:
        return None
    return name


def get_genus(species_name):
    """Extract genus from species name."""
    if not species_name:
        return None
    parts = species_name.split()
    if parts:
        return parts[0]
    return None


def point_in_bbox(lat, lon, bbox):
    """Check if point is within bounding box."""
    min_lon, min_lat, max_lon, max_lat = bbox
    return min_lat <= lat <= max_lat and min_lon <= lon <= max_lon


def expand_bbox(bbox, km):
    """Expand bounding box by km in all directions."""
    min_lon, min_lat, max_lon, max_lat = bbox
    lat_expand = km / 111.0
    center_lat = (min_lat + max_lat) / 2
    lon_expand = km / (111.0 * cos(radians(center_lat)))
    return (
        min_lon - lon_expand,
        min_lat - lat_expand,
        max_lon + lon_expand,
        max_lat + lat_expand
    )


# =============================================================================
# GEOGRAPHIC OCCURRENCE VALIDATION
# =============================================================================

async def query_gbif_in_bbox(session, species_name, bbox, semaphore, limit=500):
    """Query GBIF for species occurrences within a bounding box."""
    async with semaphore:
        try:
            min_lon, min_lat, max_lon, max_lat = bbox
            url = "https://api.gbif.org/v1/occurrence/search"
            params = {
                'scientificName': species_name,
                'hasCoordinate': 'true',
                'decimalLatitude': f'{min_lat},{max_lat}',
                'decimalLongitude': f'{min_lon},{max_lon}',
                'limit': limit
            }
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=30)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    return data.get('count', 0), data.get('results', [])
        except Exception:
            pass
        return 0, []


async def query_obis_in_bbox(session, species_name, bbox, semaphore, limit=500):
    """Query OBIS for species occurrences within a bounding box."""
    async with semaphore:
        try:
            min_lon, min_lat, max_lon, max_lat = bbox
            url = "https://api.obis.org/v3/occurrence"
            # OBIS uses geometry parameter with WKT
            wkt = f"POLYGON(({min_lon} {min_lat},{max_lon} {min_lat},{max_lon} {max_lat},{min_lon} {max_lat},{min_lon} {min_lat}))"
            params = {
                'scientificname': species_name,
                'geometry': wkt,
                'size': limit
            }
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=30)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    results = data.get('results', [])
                    total = data.get('total', len(results))
                    return total, results
        except Exception:
            pass
        return 0, []


async def query_gbif_global(session, species_name, semaphore, limit=100):
    """Query GBIF globally for nearest occurrence (used when no local records found)."""
    async with semaphore:
        try:
            url = "https://api.gbif.org/v1/occurrence/search"
            params = {
                'scientificName': species_name,
                'hasCoordinate': 'true',
                'limit': limit
            }
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=30)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    return data.get('count', 0), data.get('results', [])
        except Exception:
            pass
        return 0, []


async def query_obis_global(session, species_name, semaphore, limit=100):
    """Query OBIS globally for nearest occurrence (used when no local records found)."""
    async with semaphore:
        try:
            url = "https://api.obis.org/v3/occurrence"
            params = {
                'scientificname': species_name,
                'size': limit
            }
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=30)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    results = data.get('results', [])
                    total = data.get('total', len(results))
                    return total, results
        except Exception:
            pass
        return 0, []


async def query_gbif_congeners(session, genus, expanded_bbox, semaphore):
    """Query GBIF for congeners in expanded region."""
    async with semaphore:
        congeners = set()
        genus_key = None
        
        try:
            # Approach 1: Get genus key with rank=GENUS
            url = "https://api.gbif.org/v1/species/match"
            params = {'name': genus, 'rank': 'GENUS'}
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    genus_key = data.get('usageKey')
            
            # Approach 2: Try without rank restriction
            if not genus_key:
                params = {'name': genus}
                async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        if data.get('rank') == 'GENUS':
                            genus_key = data.get('usageKey')
            
            # Approach 3: Search for a species in the genus and extract genusKey
            if not genus_key:
                # Try common species patterns
                for suffix in ['us', 'a', 'um', 'is', 'flesus', 'harengus', 'morhua']:
                    test_species = f'{genus} {suffix}'
                    params = {'name': test_species}
                    async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=10)) as resp:
                        if resp.status == 200:
                            data = await resp.json()
                            if data.get('genusKey') and data.get('matchType') != 'NONE':
                                genus_key = data.get('genusKey')
                                break
            
            # Approach 4: Direct species database search
            if not genus_key:
                search_url = "https://api.gbif.org/v1/species/search"
                params = {'q': genus, 'rank': 'GENUS', 'limit': 5}
                async with session.get(search_url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        for result in data.get('results', []):
                            if result.get('rank') == 'GENUS':
                                canonical = result.get('canonicalName', '')
                                if canonical.lower() == genus.lower():
                                    genus_key = result.get('key')
                                    break
            
            if not genus_key:
                return congeners
            
            # Get children (species in genus)
            url = f"https://api.gbif.org/v1/species/{genus_key}/children"
            params = {'limit': 200}
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    for child in data.get('results', []):
                        if child.get('rank') == 'SPECIES':
                            name = child.get('canonicalName') or child.get('scientificName', '').split()[0:2]
                            if isinstance(name, list):
                                name = ' '.join(name)
                            if name:
                                congeners.add(name)
            
            # Also try faceted search in region
            min_lon, min_lat, max_lon, max_lat = expanded_bbox
            url = "https://api.gbif.org/v1/occurrence/search"
            params = {
                'genusKey': genus_key,
                'hasCoordinate': 'true',
                'decimalLatitude': f'{min_lat},{max_lat}',
                'decimalLongitude': f'{min_lon},{max_lon}',
                'facet': 'speciesKey',
                'facetLimit': 500,
                'limit': 0
            }
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    for facet in data.get('facets', []):
                        if facet.get('field') == 'SPECIES_KEY':
                            for count in facet.get('counts', []):
                                species_key = count.get('name')
                                if species_key:
                                    # Get species name from key
                                    sp_url = f"https://api.gbif.org/v1/species/{species_key}"
                                    async with session.get(sp_url, timeout=aiohttp.ClientTimeout(total=10)) as sp_resp:
                                        if sp_resp.status == 200:
                                            sp_data = await sp_resp.json()
                                            name = sp_data.get('canonicalName') or sp_data.get('species')
                                            if name:
                                                congeners.add(name)
        except asyncio.TimeoutError:
            pass  # Timeout - return what we have
        except aiohttp.ClientError:
            pass  # Network error - return what we have
        except Exception as e:
            # Log unexpected errors (optional: could add logging here)
            pass
        return congeners


async def validate_single_species(session, species, bbox, buffer_bbox, threshold_bbox, congener_bbox, center_lat, center_lon, semaphore):
    """
    Validate a single species using tiered geographic queries.
    
    Strategy:
    1. Check study area first → CONFIRMED
    2. Check buffer area → PLAUSIBLE  
    3. Check threshold area → POSSIBLE
    4. Query globally for nearest record → DOUBTFUL or NO_RECORDS
    """
    cleaned = clean_species_name(species)
    if not cleaned:
        return species, {
            'query_name': None,
            'n_records': 0,
            'assessment': 'INVALID_NAME',
            'min_distance_km': None,
            'congeners_in_region': ''
        }
    
    # Tier 1: Check study area (CONFIRMED)
    gbif_count, gbif_results = await query_gbif_in_bbox(session, cleaned, bbox, semaphore)
    obis_count, obis_results = await query_obis_in_bbox(session, cleaned, bbox, semaphore)
    
    if gbif_count > 0 or obis_count > 0:
        return species, {
            'query_name': cleaned,
            'n_records': gbif_count + obis_count,
            'n_gbif': gbif_count,
            'n_obis': obis_count,
            'min_distance_km': 0,
            'in_study_area': True,
            'in_buffer': True,
            'assessment': 'CONFIRMED',
            'congeners_in_region': ''
        }
    
    # Tier 2: Check buffer area (PLAUSIBLE)
    gbif_count, gbif_results = await query_gbif_in_bbox(session, cleaned, buffer_bbox, semaphore)
    obis_count, obis_results = await query_obis_in_bbox(session, cleaned, buffer_bbox, semaphore)
    
    if gbif_count > 0 or obis_count > 0:
        # Calculate approximate distance (edge of study area to nearest record)
        coords = []
        for r in gbif_results:
            lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
            if lat and lon:
                coords.append((lat, lon))
        for r in obis_results:
            lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
            if lat and lon:
                coords.append((lat, lon))
        
        min_dist = float('inf')
        for lat, lon in coords:
            dist = haversine_distance(center_lat, center_lon, lat, lon)
            min_dist = min(min_dist, dist)
        
        # Get congeners for context
        genus = get_genus(cleaned)
        congeners_str = ''
        if genus:
            congeners = await query_gbif_congeners(session, genus, congener_bbox, semaphore)
            congeners.discard(cleaned)
            if congeners:
                congeners_str = '; '.join(sorted(congeners))
        
        return species, {
            'query_name': cleaned,
            'n_records': gbif_count + obis_count,
            'n_gbif': gbif_count,
            'n_obis': obis_count,
            'min_distance_km': min_dist if min_dist != float('inf') else None,
            'in_study_area': False,
            'in_buffer': True,
            'assessment': 'PLAUSIBLE',
            'congeners_in_region': congeners_str
        }
    
    # Tier 3: Check threshold area (POSSIBLE)
    gbif_count, gbif_results = await query_gbif_in_bbox(session, cleaned, threshold_bbox, semaphore)
    obis_count, obis_results = await query_obis_in_bbox(session, cleaned, threshold_bbox, semaphore)
    
    if gbif_count > 0 or obis_count > 0:
        coords = []
        for r in gbif_results:
            lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
            if lat and lon:
                coords.append((lat, lon))
        for r in obis_results:
            lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
            if lat and lon:
                coords.append((lat, lon))
        
        min_dist = float('inf')
        for lat, lon in coords:
            dist = haversine_distance(center_lat, center_lon, lat, lon)
            min_dist = min(min_dist, dist)
        
        genus = get_genus(cleaned)
        congeners_str = ''
        if genus:
            congeners = await query_gbif_congeners(session, genus, congener_bbox, semaphore)
            congeners.discard(cleaned)
            if congeners:
                congeners_str = '; '.join(sorted(congeners))
        
        return species, {
            'query_name': cleaned,
            'n_records': gbif_count + obis_count,
            'n_gbif': gbif_count,
            'n_obis': obis_count,
            'min_distance_km': min_dist if min_dist != float('inf') else None,
            'in_study_area': False,
            'in_buffer': False,
            'assessment': 'POSSIBLE',
            'congeners_in_region': congeners_str
        }
    
    # Tier 4: Query globally to find nearest occurrence (DOUBTFUL or NO_RECORDS)
    gbif_count, gbif_results = await query_gbif_global(session, cleaned, semaphore)
    obis_count, obis_results = await query_obis_global(session, cleaned, semaphore)
    
    total_records = gbif_count + obis_count
    
    if total_records == 0:
        return species, {
            'query_name': cleaned,
            'n_records': 0,
            'n_gbif': 0,
            'n_obis': 0,
            'assessment': 'NO_RECORDS',
            'min_distance_km': None,
            'congeners_in_region': ''
        }
    
    # Calculate distance to nearest global record
    coords = []
    for r in gbif_results:
        lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
        if lat and lon:
            coords.append((lat, lon))
    for r in obis_results:
        lat, lon = r.get('decimalLatitude'), r.get('decimalLongitude')
        if lat and lon:
            coords.append((lat, lon))
    
    min_dist = float('inf')
    nearest_coord = None
    for lat, lon in coords:
        dist = haversine_distance(center_lat, center_lon, lat, lon)
        if dist < min_dist:
            min_dist = dist
            nearest_coord = (lat, lon)
    
    # Get congeners
    genus = get_genus(cleaned)
    congeners_str = ''
    if genus:
        congeners = await query_gbif_congeners(session, genus, congener_bbox, semaphore)
        congeners.discard(cleaned)
        if congeners:
            congeners_str = '; '.join(sorted(congeners))
    
    return species, {
        'query_name': cleaned,
        'n_records': total_records,
        'n_gbif': gbif_count,
        'n_obis': obis_count,
        'min_distance_km': min_dist if min_dist != float('inf') else None,
        'nearest_lat': nearest_coord[0] if nearest_coord else None,
        'nearest_lon': nearest_coord[1] if nearest_coord else None,
        'in_study_area': False,
        'in_buffer': False,
        'assessment': 'DOUBTFUL',
        'congeners_in_region': congeners_str
    }


async def validate_species_geography(species_list, bbox):
    """Validate species against GBIF/OBIS occurrence data using tiered geographic queries."""
    results = {}
    semaphore = asyncio.Semaphore(CONFIG['max_concurrent_requests'])
    
    # Calculate bounding boxes for each tier
    min_lon, min_lat, max_lon, max_lat = bbox
    center_lat = (min_lat + max_lat) / 2
    center_lon = (min_lon + max_lon) / 2
    
    buffer_bbox = expand_bbox(bbox, CONFIG['buffer_km'])
    threshold_bbox = expand_bbox(bbox, CONFIG['distance_threshold_km'])
    congener_bbox = expand_bbox(bbox, CONFIG['congener_search_buffer_km'])
    
    total = len(species_list)
    
    print(f"  Validating {total} unique species against GBIF/OBIS...")
    
    async with aiohttp.ClientSession() as session:
        # Create progress bar
        pbar = tqdm(total=total, desc="  Geographic validation", unit="species")
        
        # Process species in batches for better progress tracking
        batch_size = 10
        for i in range(0, total, batch_size):
            batch = species_list[i:i+batch_size]
            tasks = [
                validate_single_species(
                    session, species, bbox, buffer_bbox, threshold_bbox, 
                    congener_bbox, center_lat, center_lon, semaphore
                )
                for species in batch
            ]
            
            batch_results = await asyncio.gather(*tasks)
            for species, result in batch_results:
                results[species] = result
            
            pbar.update(len(batch))
            
            # Update postfix with current stats
            confirmed = sum(1 for r in results.values() if r.get('assessment') == 'CONFIRMED')
            plausible = sum(1 for r in results.values() if r.get('assessment') == 'PLAUSIBLE')
            doubtful = sum(1 for r in results.values() if r.get('assessment') == 'DOUBTFUL')
            pbar.set_postfix({'CONF': confirmed, 'PLAUS': plausible, 'DOUBT': doubtful})
        
        pbar.close()
    
    return results


# =============================================================================
# DATABASE COVERAGE VALIDATION
# =============================================================================

def parse_midori_fasta(fasta_path, target_genera):
    """Parse MIDORI2 FASTA and index by species for target genera only."""
    db_index = defaultdict(list)
    
    # Count lines first for progress
    print(f"    Counting sequences...")
    n_seqs = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                n_seqs += 1
    print(f"    Found {n_seqs:,} sequences in database")
    
    print(f"    Indexing sequences for {len(target_genera)} target genera...")
    indexed = 0
    
    with open(fasta_path, 'r') as f:
        current_header = None
        current_seq = []
        
        pbar = tqdm(total=n_seqs, desc="    Parsing MIDORI2", unit="seq")
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header and current_seq:
                    species = extract_species_from_midori_header(current_header)
                    if species:
                        genus = species.split()[0] if species else None
                        if genus in target_genera:
                            db_index[species].append(''.join(current_seq))
                            indexed += 1
                
                pbar.update(1)
                pbar.set_postfix({'indexed': indexed})
                
                current_header = line[1:]  # Remove >
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget last sequence
        if current_header and current_seq:
            species = extract_species_from_midori_header(current_header)
            if species:
                genus = species.split()[0] if species else None
                if genus in target_genera:
                    db_index[species].append(''.join(current_seq))
        
        pbar.close()
    
    return db_index


def extract_species_from_midori_header(header):
    """Extract species name from MIDORI2 header format."""
    # Format: Kingdom;Phylum;...;Genus_ID;Species name_ID;
    parts = header.split(';')
    if len(parts) >= 2:
        # Last non-empty part should be species
        species_part = parts[-2] if parts[-1] == '' else parts[-1]
        # Remove trailing ID number
        species = re.sub(r'_\d+$', '', species_part)
        # Clean up
        species = species.replace('_', ' ').strip()
        if ' ' in species:  # Should be binomial
            return species
    return None


def calculate_sequence_identity(seq1, seq2):
    """Calculate sequence identity using k-mer similarity as fallback."""
    # Try vsearch first if available
    try:
        result = subprocess.run(
            [CONFIG['vsearch_path'], '--version'],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            return calculate_identity_vsearch(seq1, seq2)
    except FileNotFoundError:
        pass
    
    # Try biopython
    try:
        from Bio import pairwise2
        alignments = pairwise2.align.localms(seq1, seq2, 2, -1, -2, -1, one_alignment_only=True)
        if alignments:
            alignment = alignments[0]
            matches = sum(1 for a, b in zip(alignment.seqA, alignment.seqB) if a == b and a != '-')
            length = max(len(seq1), len(seq2))
            return 100 * matches / length
    except ImportError:
        pass
    
    # Fallback to k-mer similarity
    k = CONFIG['kmer_size']
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    if not kmers1 or not kmers2:
        return 0
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return 100 * intersection / union if union > 0 else 0


def calculate_identity_vsearch(seq1, seq2):
    """Calculate identity using vsearch."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as q:
        q.write(f">query\n{seq1}\n")
        query_file = q.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as t:
        t.write(f">target\n{seq2}\n")
        target_file = t.name
    
    try:
        result = subprocess.run(
            [CONFIG['vsearch_path'], '--usearch_global', query_file, '--db', target_file,
             '--id', '0.5', '--maxaccepts', '1', '--userout', '/dev/stdout',
             '--userfields', 'id', '--quiet'],
            capture_output=True, text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            return float(result.stdout.strip())
    except Exception:
        pass
    finally:
        Path(query_file).unlink(missing_ok=True)
        Path(target_file).unlink(missing_ok=True)
    
    return 0


async def validate_database_coverage(taxonomy_df, sequences, midori_path, bbox, species_column):
    """Validate taxonomy assignments against MIDORI2 database."""
    results = []
    
    # Get unique genera to index - from BOTH species assignments and genus column
    genera = set()
    for species in taxonomy_df[species_column].dropna():
        genus = get_genus(str(species))
        if genus:
            genera.add(genus)
    
    # Also add genera from genus-only entries (where species is NA but genus exists)
    if 'genus' in taxonomy_df.columns:
        for idx, row in taxonomy_df.iterrows():
            species = row.get(species_column, '')
            if pd.isna(species) or not str(species).strip():
                genus = row.get('genus', '')
                if genus and pd.notna(genus) and str(genus).strip():
                    genera.add(str(genus).strip())
    
    print(f"  Found {len(genera)} unique genera to check")
    
    # Parse MIDORI2 for relevant genera
    print(f"  Parsing MIDORI2 database (this may take a minute)...")
    db_index = parse_midori_fasta(midori_path, genera)
    print(f"  Indexed {len(db_index)} species from MIDORI2")
    
    # Get congeners for each genus from GBIF
    print(f"  Querying GBIF for regional congeners...")
    congener_bbox = expand_bbox(bbox, CONFIG['congener_search_buffer_km'])
    semaphore = asyncio.Semaphore(CONFIG['max_concurrent_requests'])
    
    genus_congeners = {}
    genera_list = list(genera)
    
    async with aiohttp.ClientSession() as session:
        # Process in batches with progress bar
        pbar = tqdm(total=len(genera_list), desc="  Congener lookup", unit="genera")
        batch_size = 10
        for i in range(0, len(genera_list), batch_size):
            batch = genera_list[i:i+batch_size]
            tasks = [query_gbif_congeners(session, genus, congener_bbox, semaphore) for genus in batch]
            batch_results = await asyncio.gather(*tasks)
            for genus, congeners in zip(batch, batch_results):
                genus_congeners[genus] = congeners
            pbar.update(len(batch))
        pbar.close()
    
    # Validate each ASV
    print(f"  Validating ASV sequences...")
    total_asvs = len(taxonomy_df)
    
    for idx, row in tqdm(taxonomy_df.iterrows(), total=total_asvs, desc="  Sequence validation", unit="ASV"):
        asv_id = row.get('ASV_ID', f'ASV{idx}')
        species = row.get(species_column, '')
        genus_from_tax = row.get('genus', '')
        asv_seq = sequences.get(asv_id, '')
        
        # Handle genus-only sequences (no species epithet)
        if pd.isna(species) or not str(species).strip():
            # Try to get best match using genus from taxonomy
            if genus_from_tax and pd.notna(genus_from_tax) and str(genus_from_tax).strip() and asv_seq:
                genus = str(genus_from_tax).strip()
                
                # Find all species in this genus from MIDORI2
                genus_species = [sp for sp in db_index.keys() if get_genus(sp) == genus]
                
                if genus_species:
                    # Find best matching species
                    best_species = None
                    best_pct = 0
                    for sp in genus_species:
                        for ref_seq in db_index[sp][:3]:
                            pct = calculate_sequence_identity(asv_seq, ref_seq)
                            if pct > best_pct:
                                best_pct = pct
                                best_species = sp
                    
                    results.append({
                        'ASV_ID': asv_id,
                        'assigned_species': None,
                        'genus': genus,
                        'assigned_in_db': False,
                        'assigned_n_refs': 0,
                        'assigned_pct_identity': 0,
                        'n_congeners_in_region': 0,
                        'n_congeners_in_db': len(genus_species),
                        'n_congeners_missing': 0,
                        'congeners_in_db': '; '.join(sorted(genus_species)[:10]),
                        'congeners_missing': '',
                        'best_congener': best_species,
                        'best_congener_pct_identity': round(best_pct, 1) if best_species else None,
                        'validation_flag': 'GENUS_ONLY'
                    })
                else:
                    results.append({
                        'ASV_ID': asv_id,
                        'assigned_species': None,
                        'genus': genus,
                        'validation_flag': 'GENUS_ONLY_NO_REF'
                    })
            else:
                results.append({
                    'ASV_ID': asv_id,
                    'assigned_species': None,
                    'validation_flag': None
                })
            continue
        
        species = str(species).strip()
        genus = get_genus(species)
        asv_seq = sequences.get(asv_id, '')
        
        if not asv_seq:
            results.append({
                'ASV_ID': asv_id,
                'assigned_species': species,
                'validation_flag': 'NO_SEQUENCE'
            })
            continue
        
        # Check if assigned species is in MIDORI2
        assigned_in_db = species in db_index
        assigned_refs = db_index.get(species, [])
        
        # Calculate identity to assigned species
        assigned_pct = 0
        if assigned_refs:
            for ref_seq in assigned_refs[:5]:  # Check up to 5 refs
                pct = calculate_sequence_identity(asv_seq, ref_seq)
                assigned_pct = max(assigned_pct, pct)
        
        # Get regional congeners
        regional_congeners = genus_congeners.get(genus, set())
        congeners_in_db = [c for c in regional_congeners if c in db_index]
        congeners_missing = [c for c in regional_congeners if c not in db_index and c != species]
        
        # Find best matching congener
        best_congener = None
        best_congener_pct = 0
        for congener in congeners_in_db:
            if congener == species:
                continue
            for ref_seq in db_index[congener][:3]:
                pct = calculate_sequence_identity(asv_seq, ref_seq)
                if pct > best_congener_pct:
                    best_congener_pct = pct
                    best_congener = congener
        
        # Determine validation flag
        if not assigned_in_db:
            flag = 'NO_REF'
        elif not regional_congeners:
            flag = 'NO_CONGENERS'
        elif best_congener_pct > assigned_pct + 0.5:
            flag = 'REASSIGN'
        elif len(congeners_missing) > 0 and assigned_pct < 100:
            flag = 'UNCERTAIN'
        elif assigned_pct >= 99:
            flag = 'CONFIDENT'
        else:
            flag = 'LIKELY'
        
        results.append({
            'ASV_ID': asv_id,
            'assigned_species': species,
            'genus': genus,
            'assigned_in_db': assigned_in_db,
            'assigned_n_refs': len(assigned_refs),
            'assigned_pct_identity': round(assigned_pct, 1),
            'n_congeners_in_region': len(regional_congeners),
            'n_congeners_in_db': len(congeners_in_db),
            'n_congeners_missing': len(congeners_missing),
            'congeners_in_db': '; '.join(sorted(congeners_in_db)[:10]),
            'congeners_missing': '; '.join(sorted(congeners_missing)[:10]),
            'best_congener': best_congener,
            'best_congener_pct_identity': round(best_congener_pct, 1) if best_congener else None,
            'validation_flag': flag
        })
    
    return pd.DataFrame(results)


# =============================================================================
# COMBINED DECISION LOGIC
# =============================================================================

def determine_final_decision(row, geo_results=None):
    """Combine geographic and database assessments into final decision."""
    occ = row.get('occ_assessment')
    db = row.get('validation_flag')
    
    # Handle None and NaN values
    pct_id = row.get('assigned_pct_identity')
    pct_id = 0 if pd.isna(pct_id) else pct_id
    
    best_congener = row.get('best_congener')
    best_congener = '' if pd.isna(best_congener) else str(best_congener)
    
    best_congener_pct = row.get('best_congener_pct_identity')
    best_congener_pct = 0 if pd.isna(best_congener_pct) else best_congener_pct
    
    min_dist = row.get('occ_min_distance_km')
    min_dist = 0 if pd.isna(min_dist) else min_dist
    
    congeners_in_region = row.get('occ_congeners_in_region')
    congeners_in_region = '' if pd.isna(congeners_in_region) else str(congeners_in_region)
    
    # No species-level ID
    if pd.isna(occ) and pd.isna(db):
        return ('ALREADY_GENUS', 'N/A', 'No species-level identification')
    if pd.isna(db):
        return ('ALREADY_GENUS', 'N/A', 'No species-level identification')
    
    # Count local congeners - use BOTH sources
    # occ_congeners_in_region comes from geographic validation (GBIF/OBIS for study area)
    # n_congeners_in_db comes from database validation (congeners in MIDORI2 that are also in regional GBIF)
    n_local_congeners = 0
    local_congener_list = []
    
    # First try geographic validation congeners
    if congeners_in_region and str(congeners_in_region).strip():
        local_congener_list = [c.strip() for c in str(congeners_in_region).split('; ') if c.strip()]
        n_local_congeners = len(local_congener_list)
    
    # If geographic didn't find congeners, use database validation count
    # (this happens when assigned species is CONFIRMED but congener matches better)
    if n_local_congeners == 0:
        n_congeners_db = row.get('n_congeners_in_db')
        if pd.notna(n_congeners_db) and n_congeners_db > 0:
            n_local_congeners = int(n_congeners_db)
            # Try to get the congener list from database validation
            congeners_in_db = row.get('congeners_in_db', '')
            if pd.notna(congeners_in_db) and str(congeners_in_db).strip():
                local_congener_list = [c.strip() for c in str(congeners_in_db).split('; ') if c.strip()]
    
    # Check if best_congener is geographically confirmed/plausible using geo_results
    def is_congener_geographically_local(congener_name):
        """Check if a congener is in the study area (CONFIRMED or PLAUSIBLE)."""
        if not congener_name or not geo_results:
            return False
        geo_info = geo_results.get(congener_name, {})
        assessment = geo_info.get('assessment', '')
        return assessment in ('CONFIRMED', 'PLAUSIBLE')
    
    MIN_SPECIES_PCT = CONFIG['min_species_pct']
    MIN_GENUS_PCT = CONFIG['min_genus_pct']
    MIN_FAMILY_PCT = CONFIG['min_family_pct']
    CF_THRESHOLD_PCT = CONFIG['cf_threshold_pct']
    REASSIGN_DIFF_PCT = CONFIG['reassign_diff_pct']
    
    # === FIRST CHECK: Sequence identity thresholds ===
    # If sequence identity is too low, we cannot make a confident species-level ID
    # regardless of geographic status
    if pct_id and pct_id < MIN_SPECIES_PCT:
        # Check if we can use cf. notation for single local congener
        # This applies when: assigned species is wrong/doubtful, but only ONE local congener exists
        if n_local_congeners == 1 and best_congener and best_congener_pct:
            single_congener = local_congener_list[0] if local_congener_list else best_congener
            
            if best_congener_pct >= CF_THRESHOLD_PCT:
                # Good enough match to the single local congener for cf. notation
                return ('DROP_TO_GENUS', 'MEDIUM',
                        f'Only 1 local congener ({single_congener}); sequence {best_congener_pct:.1f}% - cf. notation')
            elif best_congener_pct >= MIN_GENUS_PCT:
                # Moderate match - genus level, note the single congener
                return ('DROP_TO_GENUS', 'MEDIUM',
                        f'Only 1 local congener ({single_congener}) but match only {best_congener_pct:.1f}%')
            elif best_congener_pct >= MIN_FAMILY_PCT:
                return ('DROP_TO_FAMILY', 'HIGH',
                        f'Sequence only {best_congener_pct:.1f}% - below genus threshold ({MIN_GENUS_PCT}%)')
            else:
                return ('DROP_TO_ORDER', 'HIGH',
                        f'Sequence only {best_congener_pct:.1f}% - below family threshold ({MIN_FAMILY_PCT}%)')
        
        # No single local congener - apply standard thresholds
        if pct_id >= MIN_GENUS_PCT:
            # Genus-level confidence
            if occ in ('CONFIRMED', 'PLAUSIBLE'):
                return ('DROP_TO_GENUS', 'MEDIUM',
                        f'Geographic {occ} but sequence only {pct_id:.1f}% (below {MIN_SPECIES_PCT}% threshold)')
            else:
                return ('DROP_TO_GENUS', 'HIGH',
                        f'Geographic {occ} + sequence only {pct_id:.1f}%')
        elif pct_id >= MIN_FAMILY_PCT:
            return ('DROP_TO_FAMILY', 'HIGH',
                    f'Sequence only {pct_id:.1f}% - below genus threshold ({MIN_GENUS_PCT}%)')
        else:
            return ('DROP_TO_ORDER', 'HIGH',
                    f'Sequence only {pct_id:.1f}% - below family threshold ({MIN_FAMILY_PCT}%)')
    
    # === DATABASE SAYS REASSIGN ===
    if db == 'REASSIGN':
        if best_congener and best_congener_pct:
            # CRITICAL: Only reassign if the target congener is geographically present
            # Use the actual geographic validation results, not just the 2000km search zone
            congener_is_local = is_congener_geographically_local(best_congener)
            
            if best_congener_pct >= MIN_SPECIES_PCT:
                # Calculate difference between best congener and assigned species
                pct_diff = best_congener_pct - pct_id if pct_id else best_congener_pct
                
                if not congener_is_local:
                    # Best matching congener is NOT in study area - don't reassign to it
                    if occ in ('CONFIRMED', 'PLAUSIBLE') and pct_id >= MIN_SPECIES_PCT:
                        # Keep original - it's geographically OK and matches well enough
                        return ('KEEP', 'MEDIUM',
                                f'Better DB match ({best_congener} {best_congener_pct:.1f}%) not local; keeping {pct_id:.1f}% match')
                    else:
                        # Neither option is great - drop to genus
                        return ('DROP_TO_GENUS', 'MEDIUM',
                                f'Better DB match ({best_congener}) not local; assigned only {pct_id:.1f}%')
                
                # Congener IS geographically local - proceed with reassignment logic
                # Count how many congeners are actually local (not just in 2000km zone)
                n_truly_local = sum(1 for c in local_congener_list if is_congener_geographically_local(c))
                
                if n_truly_local == 1:
                    # Only one truly local congener - confident reassignment
                    if occ == 'DOUBTFUL':
                        return ('REASSIGN', 'HIGH', 
                                f'Geographic DOUBTFUL + only local congener is {best_congener} ({best_congener_pct:.1f}%)')
                    else:
                        return ('REASSIGN', 'MEDIUM',
                                f'Local congener {best_congener} ({best_congener_pct:.1f}%) matches better than assigned ({pct_id:.1f}%)')
                elif pct_diff >= REASSIGN_DIFF_PCT:
                    # Multiple local congeners BUT one is clearly better (≥threshold difference)
                    return ('REASSIGN', 'MEDIUM',
                            f'{best_congener} matches {best_congener_pct:.1f}% vs assigned {pct_id:.1f}% (+{pct_diff:.1f}%)')
                else:
                    # Multiple local congeners, too close to call
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'Multiple local congeners, matches too similar ({pct_id:.1f}% vs {best_congener_pct:.1f}%)')
            elif best_congener_pct >= CF_THRESHOLD_PCT and congener_is_local:
                # Check if this is the only local congener for cf. notation
                n_truly_local = sum(1 for c in local_congener_list if is_congener_geographically_local(c))
                if n_truly_local == 1:
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'Only 1 local congener ({best_congener}); sequence {best_congener_pct:.1f}% - cf. notation')
                else:
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'Best local match {best_congener} only {best_congener_pct:.1f}% - genus-level ID')
            elif best_congener_pct >= MIN_GENUS_PCT:
                return ('DROP_TO_GENUS', 'MEDIUM',
                        f'Best local match {best_congener} only {best_congener_pct:.1f}% - genus-level ID')
            elif best_congener_pct >= MIN_FAMILY_PCT:
                return ('DROP_TO_FAMILY', 'HIGH',
                        f'Best match only {best_congener_pct:.1f}% - too low for genus-level ID')
            else:
                return ('DROP_TO_ORDER', 'HIGH',
                        f'Best match only {best_congener_pct:.1f}% - too low for family-level ID')
        else:
            # DB says REASSIGN but we don't have congener match data
            if pct_id:
                if pct_id >= MIN_SPECIES_PCT:
                    if occ in ('CONFIRMED', 'PLAUSIBLE'):
                        return ('KEEP', 'LOW', 
                                f'DB suggested congener but assigned species matches {pct_id:.1f}% and is geographically {occ}')
                    else:
                        return ('DROP_TO_GENUS', 'MEDIUM',
                                f'Geographic {occ}; assigned species {pct_id:.1f}% but DB suggests congener')
                else:
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'Assigned species only {pct_id:.1f}%; DB suggests congener but none confirmed locally')
            else:
                if occ == 'DOUBTFUL':
                    return ('DROP_TO_GENUS', 'HIGH',
                            f'Geographic DOUBTFUL; DB suggests congener')
                else:
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'DB suggests congener but no sequence match data available')
    
    # === POSSIBLE/DOUBTFUL with local congeners matching equally ===
    if occ in ('POSSIBLE', 'DOUBTFUL') and min_dist and min_dist > 500:
        if best_congener and best_congener_pct and pct_id:
            if best_congener_pct >= pct_id - 0.5:
                if n_local_congeners == 1 and best_congener_pct >= CF_THRESHOLD_PCT:
                    # Single local congener with good match - cf. notation
                    return ('DROP_TO_GENUS', 'MEDIUM',
                            f'{min_dist:.0f}km away; only 1 local congener ({best_congener}) {best_congener_pct:.1f}% - cf. notation')
                elif n_local_congeners >= 1:
                    return ('DROP_TO_GENUS', 'HIGH',
                            f'{min_dist:.0f}km away; local congener(s) match equally - likely local species')
        if n_local_congeners > 0:
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'{min_dist:.0f}km from known range with {n_local_congeners} local congener(s)')
    
    # === DOUBTFUL ===
    if occ == 'DOUBTFUL':
        # Check 1: Single local congener with good sequence match from MIDORI2
        if n_local_congeners == 1 and best_congener and best_congener_pct and best_congener_pct >= CF_THRESHOLD_PCT:
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'Geographic DOUBTFUL; only 1 local congener ({best_congener}) {best_congener_pct:.1f}% - cf. notation')
        
        # Check 2: Single local congener exists but NOT in MIDORI2 (can't verify by sequence)
        # This is the case where GBIF shows one local congener but we can't compare sequences
        # because that congener isn't in the reference database
        elif n_local_congeners == 1 and local_congener_list and (not best_congener or not best_congener_pct):
            single_local = local_congener_list[0]
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'Geographic DOUBTFUL; only 1 local congener ({single_local}) - cf. notation (not in MIDORI2)')
        
        elif n_local_congeners > 0:
            return ('DROP_TO_GENUS', 'HIGH',
                    f'Geographic DOUBTFUL with {n_local_congeners} local congener(s) - likely misidentification')
        elif db == 'CONFIDENT' and pct_id and pct_id >= MIN_SPECIES_PCT:
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'Geographic DOUBTFUL despite {pct_id:.1f}% match - possible DB error or cryptic sp.')
        else:
            return ('DROP_TO_GENUS', 'HIGH',
                    f'Geographic DOUBTFUL + DB {db}')
    
    # === CONFIRMED + good DB match + good sequence identity ===
    # (We already filtered out low sequence identity above)
    if occ == 'CONFIRMED' and db in ('CONFIDENT', 'LIKELY'):
        if pct_id and pct_id >= 99:
            return ('KEEP', 'HIGH', f'Geographic CONFIRMED + sequence match ({pct_id:.1f}%)')
        else:
            return ('KEEP', 'MEDIUM', f'Geographic CONFIRMED + sequence {pct_id:.1f}%')
    
    # === PLAUSIBLE ===
    if occ == 'PLAUSIBLE' and db in ('CONFIDENT', 'LIKELY'):
        if pct_id and pct_id >= 99:
            return ('KEEP', 'MEDIUM', f'Geographic PLAUSIBLE + sequence match ({pct_id:.1f}%)')
        else:
            return ('KEEP', 'LOW', f'Geographic PLAUSIBLE + sequence {pct_id:.1f}%')
    
    # === UNCERTAIN DB ===
    if db == 'UNCERTAIN':
        if occ in ('CONFIRMED', 'PLAUSIBLE'):
            return ('DROP_TO_GENUS', 'LOW',
                    f'Geographic {occ} but missing congeners in DB - cannot validate')
        else:
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'Geographic {occ} + DB UNCERTAIN')
    
    # === NO_CONGENERS ===
    if db == 'NO_CONGENERS':
        if pct_id and pct_id >= MIN_SPECIES_PCT:
            if occ == 'CONFIRMED':
                return ('KEEP', 'MEDIUM', f'Geographic CONFIRMED, {pct_id:.1f}% match (monotypic/rare genus)')
            elif occ == 'PLAUSIBLE':
                return ('KEEP', 'LOW', f'Geographic PLAUSIBLE, {pct_id:.1f}% match (no congeners to compare)')
            else:
                return ('DROP_TO_GENUS', 'LOW', 
                        f'Geographic {occ}, {pct_id:.1f}% match but no congeners to compare')
        else:
            return ('DROP_TO_GENUS', 'MEDIUM',
                    f'Only {pct_id:.1f}% match with no congeners to validate')
    
    # === POSSIBLE ===
    if occ == 'POSSIBLE':
        if db == 'CONFIDENT' and pct_id and pct_id >= 99:
            if min_dist and min_dist > 500:
                return ('DROP_TO_GENUS', 'LOW',
                        f'{min_dist:.0f}km from range - sequence match may reflect DB errors')
            return ('KEEP', 'LOW', f'Geographic edge case but perfect sequence match ({pct_id:.1f}%)')
        else:
            return ('DROP_TO_GENUS', 'LOW',
                    f'Geographic POSSIBLE + DB {db} ({pct_id:.1f}%)')
    
    return ('DROP_TO_GENUS', 'LOW', f'Unhandled case: occ={occ}, db={db}')


# =============================================================================
# HABITAT LOOKUP (WoRMS)
# =============================================================================

async def query_worms_habitats(taxa_list, max_concurrent=10):
    """Query WoRMS API for habitat information."""
    import urllib.parse
    
    habitat_cache = {}
    semaphore = asyncio.Semaphore(max_concurrent)
    total = len(taxa_list)
    pbar = tqdm(total=total, desc="  WoRMS habitat lookup", unit="taxa")
    
    async def query_single(session, name):
        async with semaphore:
            try:
                encoded_name = urllib.parse.quote(name)
                url = f"https://www.marinespecies.org/rest/AphiaRecordsByName/{encoded_name}"
                params = {'like': 'false', 'marine_only': 'false'}
                
                async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=15)) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        if data and len(data) > 0:
                            record = data[0]
                            habitats = []
                            if record.get('isMarine', 0) == 1:
                                habitats.append('marine')
                            if record.get('isBrackish', 0) == 1:
                                habitats.append('brackish')
                            if record.get('isFreshwater', 0) == 1:
                                habitats.append('freshwater')
                            if record.get('isTerrestrial', 0) == 1:
                                habitats.append('terrestrial')
                            if habitats:
                                habitat_cache[name] = '/'.join(habitats)
                
                await asyncio.sleep(CONFIG['worms_rate_limit'])
            except Exception:
                pass
            finally:
                pbar.update(1)
    
    try:
        async with aiohttp.ClientSession() as session:
            tasks = [query_single(session, name) for name in taxa_list]
            await asyncio.gather(*tasks, return_exceptions=True)
    except Exception as e:
        print(f"\n  WoRMS connection error: {e}")
    finally:
        pbar.close()
    
    return habitat_cache


def lookup_habitats(df):
    """Look up habitat for each taxon via WoRMS."""
    taxa_to_query = set()
    for idx, row in df.iterrows():
        species = row.get('species_raw', '')
        genus = row.get('genus', '')
        family = row.get('family', '')
        order = row.get('order', '')
        
        if species and pd.notna(species):
            parts = str(species).strip().split()[0:2]
            if len(parts) >= 2:
                taxa_to_query.add(' '.join(parts))
        if genus and pd.notna(genus):
            taxa_to_query.add(str(genus).strip())
        if family and pd.notna(family):
            taxa_to_query.add(str(family).strip())
        if order and pd.notna(order):
            taxa_to_query.add(str(order).strip())
    
    print(f"  Querying WoRMS for {len(taxa_to_query)} unique taxa...")
    
    try:
        habitat_cache = asyncio.run(query_worms_habitats(list(taxa_to_query)))
        n_found = sum(1 for v in habitat_cache.values() if v is not None)
        print(f"  Found habitat info for {n_found}/{len(taxa_to_query)} taxa")
    except Exception as e:
        print(f"  WoRMS query failed: {e}")
        habitat_cache = {}
    
    # Map to rows
    habitats = []
    for idx, row in df.iterrows():
        species = row.get('species_raw', '')
        genus = row.get('genus', '')
        family = row.get('family', '')
        order = row.get('order', '')
        decision = row.get('final_decision', '')
        cleaned = row.get('cleaned_species', '')
        
        habitat = None
        
        if decision == 'KEEP' and species and pd.notna(species):
            parts = str(species).strip().split()[0:2]
            if len(parts) >= 2:
                name = ' '.join(parts)
                habitat = habitat_cache.get(name)
        
        elif decision == 'REASSIGN' and cleaned and pd.notna(cleaned):
            habitat = habitat_cache.get(str(cleaned).strip())
        
        elif decision == 'DROP_TO_ORDER':
            if order and pd.notna(order):
                habitat = habitat_cache.get(str(order).strip())
        
        elif decision == 'DROP_TO_FAMILY':
            if family and pd.notna(family):
                habitat = habitat_cache.get(str(family).strip())
        
        # Fallback chain
        if not habitat and genus and pd.notna(genus):
            habitat = habitat_cache.get(str(genus).strip())
        if not habitat and family and pd.notna(family):
            habitat = habitat_cache.get(str(family).strip())
        if not habitat and order and pd.notna(order):
            habitat = habitat_cache.get(str(order).strip())
        
        habitats.append(habitat or 'unknown')
    
    n_unknown = habitats.count('unknown')
    if n_unknown > 0:
        print(f"  Warning: {n_unknown} taxa have 'unknown' habitat")
    
    return habitats


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_pipeline(taxonomy_csv, sequences_fasta, midori_path, bbox, species_column, output_prefix):
    """Run the complete validation pipeline."""
    
    print(f"Input taxonomy: {taxonomy_csv}")
    print(f"Input sequences: {sequences_fasta}")
    print(f"MIDORI2 database: {midori_path}")
    print(f"Output prefix: {output_prefix}")
    print()
    
    # Load taxonomy
    print("Loading taxonomy...")
    tax_df = pd.read_csv(taxonomy_csv)
    if 'Unnamed: 0' in tax_df.columns:
        tax_df = tax_df.rename(columns={'Unnamed: 0': 'ASV_ID'})
    elif 'ASV_ID' not in tax_df.columns:
        tax_df['ASV_ID'] = [f'ASV{i+1}' for i in range(len(tax_df))]
    
    print(f"  Loaded {len(tax_df)} ASVs")
    
    # Load sequences
    print("Loading sequences...")
    sequences = {}
    current_id = None
    current_seq = []
    with open(sequences_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq)
    print(f"  Loaded {len(sequences)} sequences")
    
    # Get species list
    species_list = tax_df[species_column].dropna().unique().tolist()
    print(f"  Found {len(species_list)} unique species to validate")
    
    # Step 1: Geographic validation
    print()
    print("=" * 70)
    print("STEP 1: Geographic Occurrence Validation")
    print("=" * 70)
    
    geo_results = asyncio.run(validate_species_geography(species_list, bbox))
    
    # Map results to taxonomy
    geo_cols = ['occ_query_name', 'occ_n_records', 'occ_n_gbif', 'occ_n_obis',
                'occ_min_distance_km', 'occ_nearest_lat', 'occ_nearest_lon',
                'occ_in_study_area', 'occ_in_buffer', 'occ_assessment', 
                'occ_congeners_in_region']
    
    for col in geo_cols:
        tax_df[col] = None
    
    for idx, row in tax_df.iterrows():
        species = row.get(species_column)
        if species and species in geo_results:
            result = geo_results[species]
            tax_df.at[idx, 'occ_query_name'] = result.get('query_name')
            tax_df.at[idx, 'occ_n_records'] = result.get('n_records')
            tax_df.at[idx, 'occ_n_gbif'] = result.get('n_gbif')
            tax_df.at[idx, 'occ_n_obis'] = result.get('n_obis')
            tax_df.at[idx, 'occ_min_distance_km'] = result.get('min_distance_km')
            tax_df.at[idx, 'occ_nearest_lat'] = result.get('nearest_lat')
            tax_df.at[idx, 'occ_nearest_lon'] = result.get('nearest_lon')
            tax_df.at[idx, 'occ_in_study_area'] = result.get('in_study_area')
            tax_df.at[idx, 'occ_in_buffer'] = result.get('in_buffer')
            tax_df.at[idx, 'occ_assessment'] = result.get('assessment')
            tax_df.at[idx, 'occ_congeners_in_region'] = result.get('congeners_in_region')
    
    # Print summary
    print("\nGeographic Assessment Summary:")
    for cat in ['CONFIRMED', 'PLAUSIBLE', 'POSSIBLE', 'DOUBTFUL', 'NO_RECORDS', 'INVALID_NAME']:
        count = (tax_df['occ_assessment'] == cat).sum()
        if count > 0:
            print(f"  {cat}: {count}")
    
    # Save geographic results
    geo_output = f"{output_prefix}_geographic_validated.csv"
    tax_df.to_csv(geo_output, index=False)
    print(f"\nSaved geographic validation: {geo_output}")
    
    # Step 2: Database coverage validation
    print()
    print("=" * 70)
    print("STEP 2: Database Coverage Validation")
    print("=" * 70)
    
    db_results = asyncio.run(validate_database_coverage(tax_df, sequences, midori_path, bbox, species_column))
    
    # Save database results
    db_output = f"{output_prefix}_database_validated.csv"
    db_results.to_csv(db_output, index=False)
    print(f"\nSaved database validation: {db_output}")
    
    # Print summary
    print("\nDatabase Validation Summary:")
    for flag in ['CONFIDENT', 'LIKELY', 'UNCERTAIN', 'REASSIGN', 'NO_REF', 'NO_CONGENERS']:
        count = (db_results['validation_flag'] == flag).sum()
        if count > 0:
            print(f"  {flag}: {count}")
    
    # Step 3: Combine validations
    print()
    print("=" * 70)
    print("STEP 3: Combined Taxonomy Decisions")
    print("=" * 70)
    
    # Merge geographic and database results
    merged = tax_df.merge(
        db_results[['ASV_ID', 'assigned_pct_identity', 'validation_flag',
                   'best_congener', 'best_congener_pct_identity',
                   'n_congeners_in_db', 'congeners_in_db',
                   'n_congeners_missing', 'congeners_missing']],
        on='ASV_ID',
        how='left'
    )
    
    # Apply decision logic
    print("\nApplying decision logic...")
    decisions = merged.apply(lambda row: determine_final_decision(row, geo_results), axis=1)
    merged['final_decision'] = [d[0] for d in decisions]
    merged['decision_confidence'] = [d[1] for d in decisions]
    merged['decision_reason'] = [d[2] for d in decisions]
    
    # Create cleaned species column
    def get_cleaned_species(row):
        decision = row['final_decision']
        genus = row.get('genus', '')
        species_raw = row.get(species_column, row.get('species', ''))
        pct_id = row.get('assigned_pct_identity', 0)
        occ = row.get('occ_assessment', '')
        db = row.get('validation_flag', '')
        reason = row.get('decision_reason', '')
        best_congener = row.get('best_congener', '')
        congeners_in_region = row.get('occ_congeners_in_region', '')
        
        if decision == 'KEEP':
            return species_raw if pd.notna(species_raw) else ''
        elif decision == 'REASSIGN':
            return row.get('best_congener', '')
        elif decision == 'DROP_TO_GENUS':
            # Case 1: Single local congener cf. notation with best_congener from MIDORI2
            if 'cf. notation' in str(reason).lower():
                # Try best_congener first (from MIDORI2)
                if best_congener and pd.notna(best_congener):
                    parts = str(best_congener).split()
                    if len(parts) >= 2:
                        return f"{genus} sp. cf. {parts[1]}"
                
                # Try to extract congener from occ_congeners_in_region (from GBIF)
                if congeners_in_region and pd.notna(congeners_in_region):
                    congener_list = [c.strip() for c in str(congeners_in_region).split(';') if c.strip()]
                    if len(congener_list) == 1:
                        parts = congener_list[0].split()
                        if len(parts) >= 2:
                            return f"{genus} sp. cf. {parts[1]}"
                
                # Try to extract from decision reason (fallback)
                import re
                match = re.search(r'local congener \(([^)]+)\)', str(reason))
                if match:
                    congener_name = match.group(1)
                    parts = congener_name.split()
                    if len(parts) >= 2:
                        return f"{genus} sp. cf. {parts[1]}"
            
            # Case 2: Geographic OK but DB gaps (original species, not congener)
            geographically_ok = occ in ('CONFIRMED', 'PLAUSIBLE')
            dropped_for_db_gaps = db == 'UNCERTAIN' or 'missing congeners' in str(reason).lower()
            high_match = pct_id and pct_id >= CONFIG['cf_threshold_pct']
            if geographically_ok and dropped_for_db_gaps and high_match:
                if species_raw and pd.notna(species_raw):
                    parts = str(species_raw).split()
                    if len(parts) >= 2:
                        return f"{genus} sp. cf. {parts[1]}"
            return ''
        else:
            return ''
    
    merged['cleaned_species'] = merged.apply(get_cleaned_species, axis=1)
    merged['species_raw'] = merged[species_column]
    
    # Print decision summary
    print("\nDecision Summary:")
    for dec in ['KEEP', 'REASSIGN', 'DROP_TO_GENUS', 'DROP_TO_FAMILY', 'DROP_TO_ORDER', 'ALREADY_GENUS']:
        count = (merged['final_decision'] == dec).sum()
        if count > 0:
            pct = 100 * count / len(merged)
            symbol = {"KEEP": "✓", "REASSIGN": "→", "DROP_TO_GENUS": "↑", 
                     "DROP_TO_FAMILY": "↑↑", "DROP_TO_ORDER": "↑↑↑", "ALREADY_GENUS": "·"}.get(dec, "?")
            print(f"  {symbol} {dec}: {count} ({pct:.1f}%)")
    
    n_cf = merged['cleaned_species'].str.contains(' cf\\. ', regex=True, na=False).sum()
    if n_cf > 0:
        print(f"  ? sp. cf. (tentative): {n_cf}")
    
    # Save combined results
    combined_output = f"{output_prefix}_combined_decisions.csv"
    merged.to_csv(combined_output, index=False)
    print(f"\nSaved combined decisions: {combined_output}")
    
    # Step 4: Create phyloseq-ready output
    print()
    print("=" * 70)
    print("STEP 4: Creating Phyloseq-Ready Output")
    print("=" * 70)
    
    # Build species_final with unique naming
    taxon_counts = {}
    species_final_list = []
    
    for idx, row in merged.iterrows():
        cleaned = row.get('cleaned_species', '')
        decision = row.get('final_decision', '')
        
        if cleaned and pd.notna(cleaned) and str(cleaned).strip():
            taxon = str(cleaned).strip()
            if taxon not in taxon_counts:
                taxon_counts[taxon] = 0
            taxon_counts[taxon] += 1
            species_final_list.append((taxon, taxon_counts[taxon]))
        else:
            if decision == 'DROP_TO_ORDER':
                rank_order = ['order', 'class', 'phylum', 'kingdom']
            elif decision == 'DROP_TO_FAMILY':
                rank_order = ['family', 'order', 'class', 'phylum', 'kingdom']
            else:
                rank_order = ['genus', 'family', 'order', 'class', 'phylum', 'kingdom']
            
            found = False
            for rank in rank_order:
                val = row.get(rank, '')
                if val and pd.notna(val) and str(val).strip():
                    taxon = str(val).strip()
                    if taxon not in taxon_counts:
                        taxon_counts[taxon] = 0
                    taxon_counts[taxon] += 1
                    species_final_list.append((taxon, taxon_counts[taxon]))
                    found = True
                    break
            
            if not found:
                if 'Unknown' not in taxon_counts:
                    taxon_counts['Unknown'] = 0
                taxon_counts['Unknown'] += 1
                species_final_list.append(('Unknown', taxon_counts['Unknown']))
    
    merged['species_final'] = [f"{t} {c}" for t, c in species_final_list]
    
    # Create collapse column
    collapse_list = []
    for idx, row in merged.iterrows():
        decision = row.get('final_decision', '')
        cleaned = row.get('cleaned_species', '')
        species_final = merged.at[idx, 'species_final']
        
        if decision in ('KEEP', 'REASSIGN'):
            if cleaned and pd.notna(cleaned) and ' ' in str(cleaned):
                if 'sp. cf.' not in str(cleaned) and 'sp.' not in str(cleaned):
                    collapse_list.append(cleaned)
                    continue
        collapse_list.append(species_final)
    
    merged['species_final_collapse'] = collapse_list
    
    # Lookup habitats
    print("\nLooking up habitat information...")
    habitats = lookup_habitats(merged)
    merged['habitat'] = habitats
    
    # Create phyloseq output
    phyloseq_cols = ['ASV_ID', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                     'species_final', 'species_final_collapse', 'species_original', 'habitat',
                     'validation_decision', 'decision_reason']
    
    phyloseq_df = pd.DataFrame()
    phyloseq_df['ASV_ID'] = merged['ASV_ID']
    for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']:
        phyloseq_df[col] = merged.get(col, '')
    phyloseq_df['species_final'] = merged['species_final']
    phyloseq_df['species_final_collapse'] = merged['species_final_collapse']
    phyloseq_df['species_original'] = merged['species_raw']
    phyloseq_df['habitat'] = merged['habitat']
    phyloseq_df['validation_decision'] = merged['final_decision']
    phyloseq_df['decision_reason'] = merged['decision_reason']
    
    # Add best database match columns
    def get_best_match(row):
        """Get the best matching species in MIDORI2 (either assigned or congener)."""
        assigned = row.get('species_raw', '')
        if pd.isna(assigned):
            assigned = ''
        assigned_pct = row.get('assigned_pct_identity') 
        assigned_pct = 0 if pd.isna(assigned_pct) else assigned_pct
        
        congener = row.get('best_congener', '')
        if pd.isna(congener):
            congener = ''
        congener_pct = row.get('best_congener_pct_identity')
        congener_pct = 0 if pd.isna(congener_pct) else congener_pct
        
        if congener_pct > assigned_pct and congener:
            return congener, congener_pct
        elif assigned_pct > 0 and assigned:
            return assigned, assigned_pct
        elif congener_pct > 0 and congener:
            return congener, congener_pct
        else:
            return '', None
    
    best_matches = merged.apply(get_best_match, axis=1)
    phyloseq_df['best_db_match'] = [m[0] for m in best_matches]
    phyloseq_df['best_db_match_pct'] = [m[1] for m in best_matches]
    
    # Add species_final_redund column (strip trailing numbers)
    def strip_trailing_number(s):
        """Remove trailing integer from species name like 'Phoxinus 2' -> 'Phoxinus'"""
        if pd.isna(s) or not s:
            return ''
        s = str(s).strip()
        # Match pattern: text followed by space and one or more digits at end
        match = re.match(r'^(.+?)\s+\d+$', s)
        if match:
            return match.group(1).strip()
        return s
    
    phyloseq_df['species_final_redund'] = phyloseq_df['species_final'].apply(strip_trailing_number)
    
    phyloseq_output = f"{output_prefix}_phyloseq.csv"
    phyloseq_df.to_csv(phyloseq_output, index=False)
    print(f"Saved phyloseq taxonomy: {phyloseq_output}")
    
    # Final summary
    print()
    print("=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    n_species_before = merged['species_raw'].notna().sum()
    n_keep = (merged['final_decision'] == 'KEEP').sum()
    n_reassign = (merged['final_decision'] == 'REASSIGN').sum()
    n_dropped = n_species_before - n_keep - n_reassign
    
    print(f"Species-level IDs before:  {n_species_before}")
    print(f"  KEEP (confident):        {n_keep}")
    print(f"  REASSIGN (to congener):  {n_reassign}")
    print(f"  Dropped (genus/family):  {n_dropped}")
    print()
    print(f"Output files:")
    print(f"  {geo_output}")
    print(f"  {db_output}")
    print(f"  {combined_output}")
    print(f"  {phyloseq_output}")


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="VAL_OTU_ID - Unified eDNA Taxonomy Validation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage (simplest - direct from phyloseq RDS):
  python VAL_OTU_ID.py \\
      --input ps_chordata.rds \\
      --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \\
      --bbox 5 62 13 65

Or with base name (auto-finds _taxonomy.csv and _sequences.fasta):
  python VAL_OTU_ID.py \\
      --input ps_chordata \\
      --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \\
      --bbox 5 62 13 65

This will automatically look for:
  - starting_files/ps_chordata.rds (extracts taxonomy & sequences)
  OR
  - starting_files/ps_chordata_taxonomy.csv
  - starting_files/ps_chordata_sequences.fasta

MIDORI2 is looked up in Reference_databases/ if not found.

With custom thresholds:
  python VAL_OTU_ID.py \\
      --input ps_chordata.rds \\
      --midori MIDORI2.fasta \\
      --bbox 5 62 13 65 \\
      --min-species-pct 98 \\
      --min-genus-pct 92 \\
      --buffer-km 300
        """
    )
    
    # === Input options (simplified OR explicit) ===
    input_grp = parser.add_argument_group('Input options (use --input OR --taxonomy/--sequences)')
    input_grp.add_argument("--input", "-i",
                        help="Phyloseq RDS file (e.g., 'ps_chordata.rds') OR base name "
                             "(e.g., 'ps_chordata' looks for _taxonomy.csv and _sequences.fasta)")
    input_grp.add_argument("--taxonomy", "-t",
                        help="Full path to taxonomy CSV (alternative to --input)")
    input_grp.add_argument("--sequences", "-s",
                        help="Full path to sequences FASTA (alternative to --input)")
    input_grp.add_argument("--midori", "-m", required=True,
                        help="MIDORI2 filename (looked up in Reference_databases/) or full path")
    input_grp.add_argument("--bbox", "-b", nargs=4, type=float, required=True,
                        metavar=('minLon', 'minLat', 'maxLon', 'maxLat'),
                        help="Study area bounding box (same order as bboxfinder.com)")
    
    # === Input/output options ===
    io_opts = parser.add_argument_group('Input/output options')
    io_opts.add_argument("--species-column", "-c", default="species_raw",
                        help="Column name containing species names (default: species_raw)")
    io_opts.add_argument("--output-dir", "-o", default="results_files",
                        help="Output directory (default: results_files)")
    io_opts.add_argument("--starting-dir", default="starting_files",
                        help="Directory for input files when using --input (default: starting_files)")
    io_opts.add_argument("--reference-dir", default="Reference_databases",
                        help="Directory for MIDORI2 when filename only given (default: Reference_databases)")
    
    # === Geographic validation parameters ===
    geo_opts = parser.add_argument_group('Geographic validation parameters')
    geo_opts.add_argument("--buffer-km", type=float, default=500,
                        help="Buffer around study area for 'in range' classification in km (default: 500)")
    geo_opts.add_argument("--distance-threshold-km", type=float, default=2000,
                        help="Max distance for PLAUSIBLE/POSSIBLE categories in km (default: 2000)")
    geo_opts.add_argument("--congener-search-km", type=float, default=2000,
                        help="Search area expansion for finding congeners in km (default: 2000)")
    
    # === Sequence identity thresholds ===
    seq_opts = parser.add_argument_group('Sequence identity thresholds')
    seq_opts.add_argument("--min-species-pct", type=float, default=97.0,
                        help="Minimum %% identity for confident species-level ID (default: 97.0)")
    seq_opts.add_argument("--min-genus-pct", type=float, default=90.0,
                        help="Minimum %% identity for genus-level ID (default: 90.0)")
    seq_opts.add_argument("--min-family-pct", type=float, default=80.0,
                        help="Minimum %% identity for family-level ID; below this drops to order (default: 80.0)")
    seq_opts.add_argument("--cf-threshold-pct", type=float, default=97.0,
                        help="Minimum %% identity for 'sp. cf.' notation (default: 97.0)")
    seq_opts.add_argument("--reassign-diff-pct", type=float, default=0.5,
                        help="Minimum %% difference to reassign when multiple congeners present (default: 0.5)")
    
    # === API and performance settings ===
    api_opts = parser.add_argument_group('API and performance settings')
    api_opts.add_argument("--max-concurrent", type=int, default=10,
                        help="Maximum concurrent API requests (default: 10)")
    api_opts.add_argument("--worms-rate-limit", type=float, default=0.1,
                        help="Seconds between WoRMS requests (default: 0.1)")
    
    # === Sequence comparison settings ===
    comp_opts = parser.add_argument_group('Sequence comparison settings')
    comp_opts.add_argument("--kmer-size", type=int, default=15,
                        help="K-mer size for sequence comparison fallback (default: 15)")
    comp_opts.add_argument("--vsearch-path", default="vsearch",
                        help="Path to vsearch executable if installed (default: vsearch)")
    
    args = parser.parse_args()
    
    # Resolve input file paths
    if args.input:
        input_path = Path(args.input)
        
        # Check if input is an RDS file (direct phyloseq input)
        if args.input.endswith('.rds'):
            # Look for the RDS file
            rds_file = input_path
            if not rds_file.exists():
                rds_file = Path(args.starting_dir) / args.input
            if not rds_file.exists():
                parser.error(f"RDS file not found: {args.input}")
            
            # Extract taxonomy and sequences from phyloseq
            print("Extracting data from phyloseq RDS file...")
            taxonomy_csv, sequences_fasta = extract_from_phyloseq(rds_file, args.starting_dir)
        else:
            # Use simplified input - construct paths from base name
            taxonomy_csv = Path(args.starting_dir) / f"{args.input}_taxonomy.csv"
            sequences_fasta = Path(args.starting_dir) / f"{args.input}_sequences.fasta"
            
            # Check if files exist, or if there's an RDS file to extract from
            if not taxonomy_csv.exists() or not sequences_fasta.exists():
                # Try to find an RDS file with this name
                rds_file = Path(args.starting_dir) / f"{args.input}.rds"
                if rds_file.exists():
                    print(f"Input files not found, extracting from {rds_file}...")
                    taxonomy_csv, sequences_fasta = extract_from_phyloseq(rds_file, args.starting_dir)
    
    elif args.taxonomy and args.sequences:
        # Use explicit paths - but check starting_dir if not found
        taxonomy_csv = Path(args.taxonomy)
        if not taxonomy_csv.exists():
            taxonomy_csv = Path(args.starting_dir) / args.taxonomy
        
        sequences_fasta = Path(args.sequences)
        if not sequences_fasta.exists():
            sequences_fasta = Path(args.starting_dir) / args.sequences
    else:
        parser.error("Must specify either --input OR both --taxonomy and --sequences")
    
    # Resolve MIDORI path - check reference directory if not found
    midori_path = Path(args.midori)
    if not midori_path.exists():
        midori_path = Path(args.reference_dir) / args.midori
    
    # Validate files exist
    if not taxonomy_csv.exists():
        parser.error(f"Taxonomy file not found: {taxonomy_csv}")
    if not sequences_fasta.exists():
        parser.error(f"Sequences file not found: {sequences_fasta}")
    if not midori_path.exists():
        parser.error(f"MIDORI2 file not found: {midori_path}")
    
    # Update CONFIG with all command line args
    CONFIG['buffer_km'] = args.buffer_km
    CONFIG['distance_threshold_km'] = args.distance_threshold_km
    CONFIG['congener_search_buffer_km'] = args.congener_search_km
    CONFIG['min_species_pct'] = args.min_species_pct
    CONFIG['min_genus_pct'] = args.min_genus_pct
    CONFIG['min_family_pct'] = args.min_family_pct
    CONFIG['cf_threshold_pct'] = args.cf_threshold_pct
    CONFIG['reassign_diff_pct'] = args.reassign_diff_pct
    CONFIG['max_concurrent_requests'] = args.max_concurrent
    CONFIG['worms_rate_limit'] = args.worms_rate_limit
    CONFIG['kmer_size'] = args.kmer_size
    CONFIG['vsearch_path'] = args.vsearch_path
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Generate output prefix from input filename
    input_stem = taxonomy_csv.stem
    # Remove common suffixes
    for suffix in ['_taxonomy', '_tax', '_OTUcheck', '_check']:
        if input_stem.endswith(suffix):
            input_stem = input_stem[:-len(suffix)]
    output_prefix = output_dir / input_stem
    
    # Print configuration
    print("=" * 70)
    print("VAL_OTU_ID - eDNA Taxonomy Validation Pipeline")
    print("=" * 70)
    print("\nInput files:")
    print(f"  Taxonomy:   {taxonomy_csv}")
    print(f"  Sequences:  {sequences_fasta}")
    print(f"  MIDORI2:    {midori_path}")
    print("\nConfiguration:")
    print(f"  Species column: {args.species_column}")
    print(f"  Study area bbox: {args.bbox}")
    print(f"  Buffer (in range): {CONFIG['buffer_km']} km")
    print(f"  Distance threshold: {CONFIG['distance_threshold_km']} km")
    print(f"  Congener search: {CONFIG['congener_search_buffer_km']} km")
    print(f"  Min species %%: {CONFIG['min_species_pct']}")
    print(f"  Min genus %%: {CONFIG['min_genus_pct']}")
    print(f"  Min family %%: {CONFIG['min_family_pct']}")
    print(f"  cf. threshold %%: {CONFIG['cf_threshold_pct']}")
    print(f"  Reassign diff %%: {CONFIG['reassign_diff_pct']}")
    print()
    
    # Run pipeline
    run_pipeline(
        taxonomy_csv=str(taxonomy_csv),
        sequences_fasta=str(sequences_fasta),
        midori_path=str(midori_path),
        bbox=tuple(args.bbox),
        species_column=args.species_column,
        output_prefix=str(output_prefix)
    )


if __name__ == "__main__":
    main()
