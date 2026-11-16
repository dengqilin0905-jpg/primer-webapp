# Installation

Use Python 3.11 and a fresh virtual environment. Install the checkout with `python -m pip install -r requirements-lock.txt .`, then launch `primer-webapp`. The wheel includes configuration schemas, generic workbooks and artificial examples; it contains no private dataset.

You may also run `python -m streamlit run app/app.py` from a source checkout. Account creation initializes a separate local SQLite database. Set DB_PATH before launch to select another private database. Set PRIMER_STATE_DIR to relocate the default application state root. Projects and references are optional until a feature requires them.

Run the public tests with `python -B scripts/run_public_tests.py`. This clears private-resource overrides and isolates HOME, DB_PATH and verification storage. Run `python -B scripts/synthetic_smoke.py` for canonical design/export/verification. Tests never need an external laboratory file.

Use SNAPGENE_EXE only for an optional installed, licensed Desktop map converter. GenBank export works without Desktop. A native sequencing alignment .dna still needs a separate, verified software workflow.

The lock file records a tested dependency set including a PyArrow version compatible with the NumPy/Streamlit combination. Do not independently upgrade only one numerical/UI dependency in production; repeat install and UI acceptance in a fresh environment first.
