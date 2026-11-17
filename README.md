# Primer WebApp

Installable Streamlit software for cloning-primer design, annotated GenBank export, sequencing-order workbooks and independent sequencing verification. It runs without laboratory datasets. All project data and experimental resources belong to the operator and remain outside the software repository.

## Install and run

Python 3.11 is the tested runtime. From a source checkout:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements-lock.txt .
primer-webapp --server.address 127.0.0.1 --server.port 8501
```

Open http://127.0.0.1:8501, register a local account and log in. Add your own backbone and import a gene as FASTA, GenBank or SnapGene. Design primers and export the current result; History supports re-export. No bundled experimental backbone or identity registry is needed for custom-backbone design.

The default database is under the current user's `.local/share/primer-webapp/`, outside the install directory. `DB_PATH` overrides it. Keep the listener on localhost unless you have separately configured access control and a reverse proxy.

## Private laboratory resources

Optional validated presets, project catalogs, identity/annotation tables, migration registries and original provider workbooks load from an external JSON file:

```bash
export PRIMER_PRIVATE_CONFIG=/srv/primer-private/resources.json
export DB_PATH=/srv/primer-private/app.db
export VERIFICATION_STATE_DIR=/srv/primer-private/verification
primer-webapp --server.address 127.0.0.1 --server.port 8501
```

For PowerShell use `$env:PRIMER_PRIVATE_CONFIG="C:\primer-private\resources.json"`. See [private resource configuration](docs/PRIVATE_CONFIGURATION.md). Unconfigured catalogs are empty; required missing identity/reference evidence blocks the dependent operation. No fallback scans a personal home directory.

## Features

- Basic/advanced primer controls, explicit stop policy, user backbones and optional validated shared presets.
- Shared naming and canonical Contract Project assembly with Web/Batch parity.
- Design Batch, History and Contract sequencing-order exports; Contract provides ZIP and standalone Excel.
- ZIP result import, content-hash deduplication, explicit reference assignment and separate colony/clone/event histories.
- Direction/circular-aware independent alignment, regional QC and PASS/FAIL/REVIEW/INCOMPLETE decisions.
- Feature-preserving GBK, FASTA, alignment coordinates, differences, Excel and HTML with actual aligned bases.

Bundled workbooks are empty **generic software examples**, not current company-authorized submission templates. For formal orders configure the original provider workbook; its formatting, named ranges, validations and allowed antibiotic values are preserved. Existing primers are not redesigned during export.

Native SnapGene sequencing alignment output is not implemented. An ordinary .dna map is never presented as a sequencing alignment.

## Synthetic examples and checks

Everything named `synthetic` is artificial software data and must not be used as a biological reference.

```bash
python -B scripts/run_public_tests.py
python -B scripts/synthetic_smoke.py
python -B scripts/audit_publication.py
```

The synthetic smoke test exercises custom external presets, canonical assembly, Web/Batch parity, order exports and Target PASS without laboratory resources. See [sequencing verification](docs/SEQUENCING_VERIFICATION.md) and [installation](docs/INSTALLATION.md).

## Publication and deployment

GitHub CI tests software only. Main pushes do not deploy, no SSH credentials are used by CI, and remote_deploy.sh is disabled. Laboratory deployments are separately approved, versioned installations with external persistent data. See [deployment](DEPLOYMENT.md).
