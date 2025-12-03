#!/usr/bin/env python3
"""
Debug script to test GBIF congener query for any genus.
Run this to see what the VAL_OTU_ID.py script sees when querying GBIF.

Usage:
    python debug_gbif_congeners.py Platichthys
    python debug_gbif_congeners.py Clupea
    python debug_gbif_congeners.py Gadus
"""

import sys
import requests

def test_genus_congeners(genus_name, bbox=None):
    """Test GBIF congener query for a genus."""
    
    if bbox is None:
        # Default: Norwegian waters + 2000km buffer
        bbox = (-15, 45, 35, 80)
    
    min_lon, min_lat, max_lon, max_lat = bbox
    
    print(f"{'='*60}")
    print(f"Testing GBIF API for genus: {genus_name}")
    print(f"Bounding box: {bbox}")
    print(f"{'='*60}\n")
    
    congeners = set()
    genus_key = None
    
    # Step 1: Get genus key - try multiple approaches
    print("Step 1: Getting genus key...")
    
    # Approach 1a: With rank=GENUS
    print("  Approach 1a: With rank=GENUS...")
    try:
        url = "https://api.gbif.org/v1/species/match"
        params = {'name': genus_name, 'rank': 'GENUS'}
        resp = requests.get(url, params=params, timeout=15)
        print(f"    URL: {resp.url}")
        print(f"    Status: {resp.status_code}")
        
        if resp.status_code == 200:
            data = resp.json()
            print(f"    matchType: {data.get('matchType')}")
            print(f"    scientificName: {data.get('scientificName')}")
            print(f"    rank: {data.get('rank')}")
            genus_key = data.get('usageKey')
            print(f"    usageKey: {genus_key}")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # Approach 1b: Without rank restriction
    if not genus_key:
        print("\n  Approach 1b: Without rank restriction...")
        try:
            url = "https://api.gbif.org/v1/species/match"
            params = {'name': genus_name}
            resp = requests.get(url, params=params, timeout=15)
            print(f"    URL: {resp.url}")
            print(f"    Status: {resp.status_code}")
            
            if resp.status_code == 200:
                data = resp.json()
                print(f"    matchType: {data.get('matchType')}")
                print(f"    scientificName: {data.get('scientificName')}")
                print(f"    rank: {data.get('rank')}")
                if data.get('rank') == 'GENUS':
                    genus_key = data.get('usageKey')
                    print(f"    usageKey: {genus_key}")
                else:
                    print(f"    usageKey: {data.get('usageKey')} (but rank is {data.get('rank')}, not GENUS)")
        except Exception as e:
            print(f"    ERROR: {e}")
    
    # Approach 1c: Search for a known species and extract genus key
    if not genus_key:
        test_species = f'{genus_name} flesus'
        print(f"\n  Approach 1c: Search for '{test_species}' and get genusKey...")
        try:
            url = "https://api.gbif.org/v1/species/match"
            params = {'name': test_species}
            resp = requests.get(url, params=params, timeout=15)
            print(f"    URL: {resp.url}")
            print(f"    Status: {resp.status_code}")
            
            if resp.status_code == 200:
                data = resp.json()
                print(f"    matchType: {data.get('matchType')}")
                print(f"    scientificName: {data.get('scientificName')}")
                print(f"    rank: {data.get('rank')}")
                print(f"    genus: {data.get('genus')}")
                print(f"    genusKey: {data.get('genusKey')}")
                if data.get('genusKey'):
                    genus_key = data.get('genusKey')
        except Exception as e:
            print(f"    ERROR: {e}")
    
    # Approach 1d: Direct search in GBIF species database
    if not genus_key:
        print(f"\n  Approach 1d: Search GBIF species database for '{genus_name}'...")
        try:
            url = "https://api.gbif.org/v1/species/search"
            params = {'q': genus_name, 'rank': 'GENUS', 'limit': 5}
            resp = requests.get(url, params=params, timeout=15)
            print(f"    URL: {resp.url}")
            print(f"    Status: {resp.status_code}")
            
            if resp.status_code == 200:
                data = resp.json()
                results = data.get('results', [])
                print(f"    Found {len(results)} results")
                for r in results:
                    print(f"      - {r.get('scientificName')} (key: {r.get('key')}, rank: {r.get('rank')})")
                    if r.get('rank') == 'GENUS' and r.get('canonicalName', '').lower() == genus_name.lower():
                        genus_key = r.get('key')
                        print(f"    -> Using this key: {genus_key}")
                        break
        except Exception as e:
            print(f"    ERROR: {e}")
    
    if not genus_key:
        print(f"\n  FAILED: Could not find genus key for '{genus_name}'")
        print("  The congener query will return empty set.")
        return
    
    print(f"\n  SUCCESS: Found genus key = {genus_key}")
    
    # Step 2: Get children (species in genus) - globally
    print(f"\nStep 2: Getting children of genus {genus_key}...")
    try:
        url = f"https://api.gbif.org/v1/species/{genus_key}/children"
        params = {'limit': 200}
        resp = requests.get(url, params=params, timeout=15)
        print(f"  URL: {resp.url}")
        print(f"  Status: {resp.status_code}")
        
        if resp.status_code != 200:
            print(f"  ERROR: Bad status code")
        else:
            data = resp.json()
            results = data.get('results', [])
            print(f"  Found {len(results)} children (globally)")
            
            species_children = [r for r in results if r.get('rank') == 'SPECIES']
            print(f"  Of which {len(species_children)} are SPECIES rank:")
            
            for child in species_children:
                name = child.get('canonicalName') or child.get('scientificName', '').split()[0:2]
                if isinstance(name, list):
                    name = ' '.join(name)
                key = child.get('key')
                print(f"    - {name} (key: {key})")
                if name:
                    congeners.add(name)
                    
    except Exception as e:
        print(f"  ERROR: {type(e).__name__}: {e}")
    
    # Step 3: Occurrence facet search in region
    print(f"\nStep 3: Testing occurrence facet search in bounding box...")
    try:
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
        resp = requests.get(url, params=params, timeout=15)
        print(f"  URL: {resp.url}")
        print(f"  Status: {resp.status_code}")
        
        if resp.status_code != 200:
            print(f"  ERROR: Bad status code")
        else:
            data = resp.json()
            total = data.get('count', 0)
            print(f"  Total occurrences in region: {total}")
            
            for facet in data.get('facets', []):
                if facet.get('field') == 'SPECIES_KEY':
                    counts = facet.get('counts', [])
                    print(f"  Species facets found: {len(counts)}")
                    
                    for count in counts:
                        species_key = count.get('name')
                        n_records = count.get('count')
                        
                        # Get species name from key
                        sp_url = f"https://api.gbif.org/v1/species/{species_key}"
                        sp_resp = requests.get(sp_url, timeout=10)
                        if sp_resp.status_code == 200:
                            sp_data = sp_resp.json()
                            name = sp_data.get('canonicalName') or sp_data.get('species')
                            print(f"    - {name}: {n_records} records (key: {species_key})")
                            if name:
                                congeners.add(name)
                        else:
                            print(f"    - [key {species_key}]: {n_records} records (couldn't get name)")
                            
    except Exception as e:
        print(f"  ERROR: {type(e).__name__}: {e}")
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: Found {len(congeners)} congeners for {genus_name}")
    print(f"{'='*60}")
    for c in sorted(congeners):
        print(f"  - {c}")
    
    if len(congeners) == 0:
        print("\n  WARNING: No congeners found! This explains the NO_CONGENERS flag.")
    elif len(congeners) == 1:
        single = list(congeners)[0]
        print(f"\n  NOTE: Only 1 congener found ({single}).")
        print(f"  If this is the local species and assigned species is DOUBTFUL,")
        print(f"  cf. notation should be applied.")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python debug_gbif_congeners.py <genus_name> [min_lon min_lat max_lon max_lat]")
        print("\nExamples:")
        print("  python debug_gbif_congeners.py Platichthys")
        print("  python debug_gbif_congeners.py Platichthys 5 62 13 65")
        sys.exit(1)
    
    genus = sys.argv[1]
    
    if len(sys.argv) >= 6:
        bbox = (float(sys.argv[2]), float(sys.argv[3]), 
                float(sys.argv[4]), float(sys.argv[5]))
    else:
        # Default: expanded bbox for Norwegian waters
        bbox = (-15, 45, 35, 80)
    
    test_genus_congeners(genus, bbox)
