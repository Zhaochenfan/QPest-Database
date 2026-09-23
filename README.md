# QPD Marker Database Scripts

This repository contains scripts for building marker-specific BLAST databases and performing spatial thinning on occurrence data.

## Scripts

### 1. build_marker_blast_db.py

A Django management command that builds marker-specific BLAST databases from DNABarcode records.

**Supported Marker Families:**

| Family | Included Markers |
|--------|------------------|
| COI | COI-5P, COI-3P, COII, COXIII, COI-PSEUDO, COI-5PNMT1, COII-COI |
| ITS | ITS, ITS1, ITS2, 5-8S |
| 18S | 18S, 18S-5P, 18S-3P |
| matK | matK, matK-like |
| rbcL | rbcL, rbcLa, rbcL-like |
| 16S | 16S |
| CYTB | CYTB |
| 28S | 28S, 28S-D2, 28S-D2-D3, 28S-D1-D2, 28S-D3 |
| ND | ND1, ND2, ND3, ND4, ND4L, ND5-0, ND6, MT-NADH5 |
| trn | trnH-psbA, trnL, trnL-F, trnD-trnY-trnE, cp-trnL-intron, atpF-atpH, atpF-intron, atpB-rbcL, ndhK-ndhC, psbK-psbI, psbE-petL, rpL32-trnL, rps16-trnK |
| OTHER | Markers not belonging to any known family |

**Features:**

- Builds separate BLAST databases for each marker family
- Generates a combined ALL database (all markers merged)
- Sequence quality filtering: only retains sequences >= 200bp with valid nucleotides
- Automatic deduplication based on process_id
- Cross-database global deduplication for the ALL database

**Usage:**

```bash
python manage.py build_marker_blast_db
```

**Configuration:**

Requires the following settings in `settings.py`:

```python
BLAST_MAKEBLASTDB_BIN_PATH = '/path/to/makeblastdb'  # Path to BLAST makeblastdb binary

# Database output paths
BLAST_DB_COI = '/path/to/db/coi'
BLAST_DB_ITS = '/path/to/db/its'
BLAST_DB_18S = '/path/to/db/18s'
BLAST_DB_matK = '/path/to/db/matk'
BLAST_DB_rbcL = '/path/to/db/rbcl'
BLAST_DB_16S = '/path/to/db/16s'
BLAST_DB_CYTB = '/path/to/db/cytb'
BLAST_DB_28S = '/path/to/db/28s'
BLAST_DB_ND = '/path/to/db/nd'
BLAST_DB_trn = '/path/to/db/trn'
BLAST_DB_OTHER = '/path/to/db/other'
BLAST_DB_ALL = '/path/to/db/all'
```

---

### 2. spatial_thinned.py

A spatial thinning script based on environmental raster resolution, designed to remove overly dense sampling points from geographic occurrence data.

**Features:**

- Scientific spatial thinning using Grid Binning method
- Configurable spatial resolution (default: 2.5 arc-minutes)
- Comprehensive data quality control:
  - Removes records with missing latitude/longitude
  - Removes zero-coordinate records (0, 0)
  - Removes coordinates outside global bounds (-90 to 90 latitude, -180 to 180 longitude)
- Detailed thinning report in output

**Usage:**

```python
from spatial_thinned import scientific_spatial_thinning

scientific_spatial_thinning(
    input_csv='path/to/input.csv',
    output_csv='path/to/output.csv',
    decimalLatitude='decimalLatitude',
    decimalLongitude='decimalLongitude',
    resolution_min=2.5  # arc-minutes, default is 2.5
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| input_csv | str | Input CSV file path |
| output_csv | str | Output CSV file path |
| decimalLatitude | str | Latitude column name |
| decimalLongitude | str | Longitude column name |
| resolution_min | float | Spatial resolution in arc-minutes (default: 2.5) |

**Resolution Reference:**

| Resolution | Degrees |
|------------|---------|
| 2.5 arc-min | 0.0416667° |
| 5 arc-min | 0.0833333° |
| 10 arc-min | 0.1666667° |

**Output Report Example:**

```
--- Thinning Report (input.csv) ---
Original points: 10000
Retained points: 3500
Removed points: 6500 (65.0%)
Resolution grid: 2.5 arc-minutes (0.041667 degrees)
Output saved to: output.csv
```

---

## Background

These scripts support a DNA Barcode database project for plant quarantine applications. The database contains multiple molecular markers (COI, ITS, matK, rbcL, etc.) to support species identification and biodiversity research.
