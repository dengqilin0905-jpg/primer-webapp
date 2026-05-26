# Primer Designer (Streamlit + SQLite)

Primer Designer is a local/internal primer design tool for cloning workflows.
It keeps per-user private data, supports multi-format sequence import, batch design, Sangon Excel export, and SnapGene-compatible GenBank map export.

## Core stack

- Python
- Streamlit
- SQLite
- Entry: `app/app.py`

## Main features

- User register/login (password hashed with PBKDF2 + salt)
- Per-user backbone library:
  - `backbone name`, `HOMO_F`, `HOMO_R`, enzyme label
  - optional vector map fields:
    - `vector_sequence`
    - `insertion_mode`
    - restriction sites and junction extra bases
    - map note
- Gene import:
  - FASTA/TXT/GenBank/SnapGene `.dna`
- Design batch modes:
  - `mapping` mode (existing fine-grained control)
  - `cartesian` mode (Gene × Backbone combinations)
- Optional frame check:
  - in-frame / stop codon / junction shift warnings
- Exports:
  - Sangon Excel (`sangon_template.xlsx`, format preserved)
  - Fragment TXT (full/simple)
  - non-CDS FASTA
  - SnapGene-compatible GenBank (`.gbk`) or ZIP bundle for multi-combo runs
- History:
  - keeps each run with auto status (`pending/success/fail`)
  - supports re-export

## Local run

### 1) Create venv

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

Linux/macOS:

```bash
source .venv/bin/activate
```

### 2) Install dependencies

```bash
pip install -r requirements.txt
```

`snapgene-reader==0.1.23` is required for `.dna` parsing.

### 3) Start Streamlit

```bash
streamlit run app/app.py --server.address 127.0.0.1 --server.port 8501
```

Open:

- <http://127.0.0.1:8501>

## `.dna` import notes

- Upload `.dna` directly in Genes page (or in Backbones vector upload).
- If `.dna` parse fails, the app shows clear per-file errors and continues processing other files.
- Binary files mislabeled as `.fa/.txt` are content-checked and routed safely.

## Design workflow

1. Add backbones (with HOMO arms, optional vector sequence and restriction settings).
2. Import genes.
3. Go to Design page:
   - choose `mapping` or `cartesian` batch mode
   - set Tm/length parameters
   - optional frame-check options
4. Generate primers.
5. Export:
   - Sangon Excel
   - GenBank map (`.gbk` / ZIP)

## Basic vs Advanced design controls

- Basic mode keeps experiment-changing decisions visible:
  - stop codon removal/retention
  - Tm min / Tm max / Target Tm
  - primer length range
  - homology arm length
- Advanced mode contains GC limits, ion/Tm profile parameters, frame check, and manual vector/restriction-site details.
- When a vector sequence is uploaded and the selected restriction sites are unique, the Backbones page can auto-generate HOMO_F/HOMO_R from the selected homology-arm length.

## Canonical local QA data

Use the `test/` folder for release smoke tests:

- Backbone: `pGreenII-62SK.dna`
- Restriction sites: BamHI / EcoRI
- Expected 21 bp homology arms:
  - `CGCTCTAGAACTAGTGGATCC`
  - `GATAAGCTTGATATCGAATTC`
- CDS inserts:
  - `ZmMYBR26.dna`
  - `ZmMYBS1.dna`
- Required CDS option: keep terminal stop codon (`remove_stop=False`).

## Sangon Excel export guarantees

- Uses `app/sangon_template.xlsx` as template baseline.
- Keeps template format, formula, data validation, and headers.
- Auto-extends rows when primer count exceeds base template range.
- Writes profile fields into template header.
- Clears the full order data area before writing primers so each export contains only the current run.

## Repository layout

```text
app/
  app.py
  db.py
  design_ui.py
  history_ui.py
  backbones_ui.py
  primer_design.py
  excel_export.py
  genbank_export.py
deploy/
scripts/
```

## Quick checks

```bash
python scripts/check_indent.py
python -m py_compile app/*.py
```

