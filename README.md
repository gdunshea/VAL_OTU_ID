# VAL_OTU_ID

**Unified Eukaryote eDNA Taxonomy Identification Validation Pipeline**

Validates species identifications from eDNA metabarcoding by cross-referencing:
1. **Geographic plausibility** — GBIF/OBIS occurrence databases
2. **Sequence identity** — MIDORI2 reference database
3. **Regional biodiversity** — Local congener availability
4. **Conservation status** — IUCN Red List (optional)
5. **Invasive species** — GBIF/GRIIS databases (optional)

Detailed intructions can be found at: https://gdunshea.github.io/VAL_OTU_ID/

## Why Use This?

DADA2's `assignTaxonomy()` is a sequence analysis exercise that ignores species ranges and reference database completeness. It can misidentify species when:
- The true species is missing from the reference database
- A closely related species has a better-represented sequence
- Geographic context is ignored (assigns Pacific species to Atlantic samples)

VAL_OTU_ID corrects these issues by validating each ASV against occurrence records and re-evaluating sequence matches within the geographic context.

## Key Features

- **Direct phyloseq input** — Reads `.rds` files directly (no manual extraction needed)
- **Intelligent reassignment** — Corrects DADA2 misidentifications when a **geographically local** congener matches better
- **Geographic verification** — Only reassigns to species that are CONFIRMED or PLAUSIBLE in your study area
- **Robust GBIF queries** — Multiple fallback methods ensure congeners are found even for problematic genus names
- **Tiered decisions** — KEEP, REASSIGN, DROP_TO_GENUS, DROP_TO_FAMILY, DROP_TO_ORDER
- **cf. notation** — Uses `Genus sp. cf. species` when confident of genus but uncertain of species
- **Conservation status** — Flags IUCN Red List threatened species (CR, EN, VU)
- **Invasive species detection** — Identifies species listed in GRIIS/GISD databases
- **Habitat lookup** — Adds marine/freshwater/brackish/terrestrial from WoRMS
- **Phyloseq-ready output** — Direct import back into R

## Installation

### Requirements

- Python 3.8+
- R with `phyloseq` and `Biostrings` packages (for `.rds` input)
- Optional but recommended: `vsearch` for accurate sequence alignment

### Python Dependencies

```bash
pip install pandas aiohttp tqdm
```

### Install vsearch (recommended)

```bash
# macOS
brew install vsearch

# conda
conda install -c bioconda vsearch

# Ubuntu/Debian
sudo apt install vsearch
```

### Reference Database

Download reference databases used for dada2 assignTaxonomy(). Examples are MIDORI2 databases (http://www.reference-midori.info/) or the SILVA databases (https://www.arb-silva.de/)

Whichever databases are used, these need to be added to the "Reference_databases" directory in the repo.

## Quick Start

### Easiest Use (from phyloseq .rds)

1. Get taxonomic assignments of ASVs/OTUs using a validated reference database (e.g., MIDORI or SILVA) and `assignTaxonomy()` in DADA2
2. Construct a phyloseq object with minimum of `otu_table`, `tax_table`, `refseq` slots filled
3. Save as `.rds` file

```bash
python VAL_OTU_ID.py \
    --input ps_chordata.rds \
    --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \
    --bbox 6 63 12 66
```

### Expected Directory Structure

```
VAL_OTU/
├── VAL_OTU_ID.py
├── starting_files/
│   ├── ps_chordata.rds          # Your phyloseq object
│   └── iucn_api_token.txt       # Optional: for --check-iucn
├── Reference_databases/
│   └── MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta
└── results_files/               # Created automatically
    ├── ps_chordata_geographic_validated.csv
    ├── ps_chordata_database_validated.csv
    ├── ps_chordata_combined_decisions.csv
    └── ps_chordata_phyloseq.csv  # Ready for R import
```

### With Custom Parameters

```bash
python VAL_OTU_ID.py \
    --input ps_chordata.rds \
    --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \
    --bbox 6 63 12 66 \
    --min-species-pct 98 \
    --min-genus-pct 92 \
    --buffer-km 300
```

### With Conservation and Invasive Species Checking

```bash
python VAL_OTU_ID.py \
    --input ps_chordata.rds \
    --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \
    --bbox 6 63 12 66 \
    --check-iucn \
    --check-invasive
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Phyloseq `.rds` file OR base name of `_taxonomy.csv`/`_sequences.fasta` |
| `--midori` | Reference database FASTA file in dada2 format (filename or full path) |
| `--bbox` | Study area bounding box: `minLon minLat maxLon maxLat` |

### Geographic Validation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--buffer-km` | 500 | Buffer around study area for "PLAUSIBLE" classification |
| `--distance-threshold-km` | 2000 | Maximum distance for "POSSIBLE" classification |
| `--congener-search-km` | 2000 | Search area expansion for finding congeners |

### Sequence Identity Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-species-pct` | 97.0 | Minimum % identity for species-level ID |
| `--min-genus-pct` | 90.0 | Minimum % identity for genus-level ID |
| `--min-family-pct` | 80.0 | Minimum % identity for family-level ID |
| `--cf-threshold-pct` | 97.0 | Minimum % identity for `sp. cf.` notation |
| `--reassign-diff-pct` | 0.5 | Minimum % difference to reassign when multiple local congeners |

### Recommended Settings by Amplicon

| Marker | Length | `--min-species-pct` | `--reassign-diff-pct` |
|--------|--------|---------------------|----------------------|
| MiFish 12S | ~170bp | 97-98% | 0.5-0.6% |
| 16S V4 | ~300bp | 97% | 0.3-0.4% |
| COI | ~650bp | 97-98% | 0.15-0.2% |

### Conservation and Invasive Species (Optional)

| Parameter | Description |
|-----------|-------------|
| `--check-iucn` | Enable IUCN Red List status checking (requires API token) |
| `--check-invasive` | Enable invasive species checking via GBIF/GRIIS |

**IUCN API Token Setup:**
1. Register for free at https://apiv3.iucnredlist.org/
2. Create file `starting_files/iucn_api_token.txt` containing only your token

**Note:** 
- Conservation status is checked against `species_final_redund` (the validated final species name)
- For cf. notation (e.g., "Platichthys sp. cf. flesus"), the implied species ("Platichthys flesus") is queried
- The original DADA2 assignment is ALWAYS checked (stored in `_original` columns) — this catches invasive/threatened species even when they're dropped to genus level
- `species_concern_flag` is set if EITHER the final OR original species is threatened/invasive

## Output Files

### `_phyloseq.csv`

Ready for R import with columns:

| Column | Description |
|--------|-------------|
| `ASV_ID` | Unique identifier for each ASV |
| `kingdom` through `genus` | Higher taxonomy (unchanged from DADA2) |
| `species_final` | Validated species with unique number for each ASV/OTU (e.g., "Gadus morhua 3") |
| `species_final_collapse` | Species without number (for aggregating replicates) |
| `species_final_redund` | Species name with trailing number stripped (for grouping) |
| `species_original` | Original DADA2 assignment |
| `habitat` | marine/freshwater/brackish/terrestrial (from WoRMS) |
| `validation_decision` | KEEP/REASSIGN/DROP_TO_GENUS/etc. |
| `decision_reason` | Human-readable explanation |
| `best_db_match` | Best matching species in MIDORI2 |
| `best_db_match_pct` | % identity to best match |
| `iucn_category` | IUCN Red List category for final species (CR/EN/VU/NT/LC/DD/NE/NA). For cf. species, queries the implied species. |
| `iucn_population_trend` | Population trend for final species (Increasing/Stable/Decreasing/Unknown) |
| `invasive_status` | Invasive status for final species: INVASIVE/INTRODUCED/NOT_LISTED/NA |
| `iucn_category_original` | IUCN category for original DADA2 assignment (all decisions with species-level assignment) |
| `invasive_status_original` | Invasive status for original DADA2 assignment (all decisions with species-level assignment) |
| `species_concern_flag` | Summary flag considering BOTH final and original: THREATENED, INVASIVE, THREATENED+INVASIVE, or NONE |

### Import to R

```r
library(phyloseq)

# Read validated taxonomy
new_tax <- read.csv("results_files/ps_chordata_phyloseq.csv", row.names = 1)

# Update phyloseq object
tax_table(ps) <- tax_table(as.matrix(new_tax))
```

## Decision Logic Overview

The pipeline makes decisions based on **three sources of evidence**:

1. **Geographic assessment** — Is the assigned species found in the study area?
2. **Sequence identity** — How well does the ASV match the reference?
3. **Local congeners** — Are there better-matching species that ARE geographically local?

### Critical Rule: Geographic Verification for Reassignment

**REASSIGN only occurs when the target congener is geographically CONFIRMED or PLAUSIBLE in your study area.** 

This prevents incorrect reassignments like:
- ASV matches *Clupea pallasii* (Pacific herring) at 99.4%
- ASV matches *Clupea harengus* (Atlantic herring) at 98.8%  
- *C. pallasii* is NOT in Norwegian waters
- **Result**: KEEP *C. harengus* (not reassign to non-local species)

### Decision Categories

| Decision | Meaning |
|----------|---------|
| **KEEP** | Original assignment retained |
| **REASSIGN** | Changed to a better-matching LOCAL congener |
| **DROP_TO_GENUS** | Species uncertain, genus retained |
| **DROP_TO_FAMILY** | Genus uncertain, family retained |
| **DROP_TO_ORDER** | Family uncertain, order retained |

See `VAL_OTU_Instructions.html` for the complete decision flowchart.

## Troubleshooting

### GBIF not finding congeners

The pipeline uses multiple fallback methods to query GBIF for genus information. Use the debug script to test:

```bash
python debug_gbif_congeners.py Platichthys
```

### "R extraction failed"

Ensure R packages are installed:
```r
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "Biostrings"))
```

### "vsearch not found"

Install vsearch for accurate sequence identity calculations:
```bash
brew install vsearch  # macOS
conda install -c bioconda vsearch  # conda
```

## Citation

If you use VAL_OTU_ID in your research, please cite:

> Dunshea, G. (2024). VAL_OTU_ID: Unified eDNA Taxonomy Validation Pipeline. https://github.com/[your-repo]

## License

MIT License

## Contributing

Issues and pull requests welcome!
