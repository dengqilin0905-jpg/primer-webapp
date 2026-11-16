"""Compatibility entrypoint: every input is explicit; no built-in project dataset."""
from run_primer_project import main
from design_config import resolve_batch_backbone
if __name__=='__main__':raise SystemExit(main())
