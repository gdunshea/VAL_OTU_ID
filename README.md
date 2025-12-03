# VAL_OTU_ID

**Unified eDNA Taxonomy Validation Pipeline**

Validates species identifications from eDNA metabarcoding by cross-referencing:
1. **Geographic plausibility** — GBIF/OBIS occurrence databases
2. **Sequence identity** — MIDORI2 reference database
3. **Regional biodiversity** — Local congener availability

## Why Use This?

DADA2's `assignTaxonomy()` can misidentify species when:
- The true species is missing from the reference database
- A closely related species has a better-represented sequence
- Geographic context is ignored (assigns Pacific species to Atlantic samples)

VAL_OTU_ID corrects these issues by validating each ASV against occurrence records and re-evaluating sequence matches within the geographic context.

## Key Features

- **Direct phyloseq input** — Reads `.rds` files directly (no manual extraction needed)
- **Intelligent reassignment** — Corrects DADA2 misidentifications when a **geographically local** congener matches better
- **Geographic verification** — Only reassigns to species that are CONFIRMED or PLAUSIBLE in your study area
- **Tiered decisions** — KEEP, REASSIGN, DROP_TO_GENUS, DROP_TO_FAMILY, DROP_TO_ORDER
- **cf. notation** — Uses `Genus sp. cf. species` when confident of genus but uncertain of species
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

Download MIDORI2 from: http://www.reference-midori.info/

Recommended: `MIDORI2_UNIQ_NUC_GB###_srRNA_DADA2.fasta`

## Quick Start

### Simplest Usage (from phyloseq .rds)

```bash
python VAL_OTU_ID.py \
    --input ps_chordata.rds \
    --midori MIDORI2_UNIQ_NUC_GB268_srRNA_DADA2.fasta \
    --bbox 5 62 13 65
```

### Expected Directory Structure

```
VAL_OTU/
├── VAL_OTU_ID.py
├── starting_files/
│   └── ps_chordata.rds          # Your phyloseq object
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
    --bbox 5 62 13 65 \
    --min-species-pct 98 \
    --min-genus-pct 92 \
    --buffer-km 300 \
    --reassign-diff-pct 0.6
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Phyloseq `.rds` file OR base name of `_taxonomy.csv`/`_sequences.fasta` |
| `--midori` | MIDORI2 reference FASTA (filename or full path) |
| `--bbox` | Study area bounding box: `minLon minLat maxLon maxLat` |

### Geographic Validation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--buffer-km` | 500 | Buffer around study area for "PLAUSIBLE" classification |
| `--distance-threshold-km` | 2000 | Maximum distance for "POSSIBLE" classification |
| `--congener-search-km` | 2000 | Search area expansion for finding congeners in MIDORI2 |

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

## Output Files

### `_phyloseq.csv`

Ready for R import with columns:

| Column | Description |
|--------|-------------|
| `species_final` | Validated species name with unique numbering (e.g., "Gadus morhua 3") |
| `species_final_collapse` | Species name without numbering (for aggregation) |
| `species_original` | Original DADA2 assignment |
| `habitat` | marine/freshwater/brackish/terrestrial (from WoRMS) |
| `validation_decision` | KEEP/REASSIGN/DROP_TO_GENUS/etc. |
| `decision_reason` | Human-readable explanation |
| `best_db_match` | Best matching species in MIDORI2 |
| `best_db_match_pct` | % identity to best match |

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

## Example Scenarios

### Scenario 1: Correct Assignment
- DADA2 assigned: *Gadus morhua* (99.5% match)
- Geographic: CONFIRMED in study area
- **Result**: KEEP

### Scenario 2: Wrong Species, Local Congener Available
- DADA2 assigned: *Pollachius virens* (98.8% match)
- Local congener: *Pollachius pollachius* (100% match, CONFIRMED in study area)
- **Result**: REASSIGN to *P. pollachius*

### Scenario 3: Better DB Match But Not Local
- DADA2 assigned: *Clupea harengus* (98.8% match, CONFIRMED locally)
- Better DB match: *Clupea pallasii* (99.4% match, NOT in study area)
- **Result**: KEEP *C. harengus* (better match is not local)

### Scenario 4: Geographic Impossibility
- DADA2 assigned: *Platichthys stellatus* (Pacific species)
- Geographic: DOUBTFUL (nearest record 8000km away)
- Local congener: *Platichthys flesus* (97.5% match)
- **Result**: DROP_TO_GENUS with cf. notation → "Platichthys sp. cf. flesus"

## Citation

If you use VAL_OTU_ID in your research, please cite:

> Dunshea, G. (2024). VAL_OTU_ID: Unified eDNA Taxonomy Validation Pipeline. https://github.com/[your-repo]

## License

MIT License

## Contributing

Issues and pull requests welcome!
