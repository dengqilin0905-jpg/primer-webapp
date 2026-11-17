# External private resources

Copy config/private_resources.example.json outside the checkout, adjust its paths, and set PRIMER_PRIVATE_CONFIG before starting Python. Relative paths resolve against that JSON file. Each setting can be overridden with PRIMER_<UPPERCASE_KEY>, for example PRIMER_BACKBONE_REGISTRY. Empty values are ignored. No configuration load creates or changes a project.

| Key | Meaning |
| --- | --- |
| projects_root | Read-only formal project/reference catalog |
| backbone_registry | Private preset registry in the existing JSON-compatible YAML schema |
| backbone_root | Base directory for the selected registry's vector paths |
| identity_registry | Formal identity master TSV |
| construct_master | Supplementary historical construct evidence TSV |
| rice_annotation / curated_annotation | Operator-provided annotation TSVs |
| migration_registry | Private historical project registration JSON |
| order_template | Original primer-order workbook |
| sequencing_template | Original current sequencing-order workbook |
| legacy_sequencing_template | Explicit registered legacy workbook, when required |
| skill_root | Optional external workflow/contract rules for portfolio auditing |
| state_dir | Default writable application state root |

Preset vector paths are absolute or relative to the registry's project root (registry parent directory's parent). Set backbone_root to that same root. The complete source sequence, hash, replacement boundaries and homology arms are still verified; configuration does not weaken scientific checks. Never replace an experimental preset with the synthetic example to silence a validation error.

The software's defaults have no shared presets and no real identity/annotation records. Upload your own reference plus explicit assembly provenance, or configure an authoritative project catalog. An absent optional catalog is empty; missing mandatory evidence remains blocking. Numbered identities must be checked against the supplied exact registry.

DB_PATH and VERIFICATION_STATE_DIR are separate runtime environment settings. Verification ledgers are user scoped and separate from the primer database. Use private directory permissions, include both stores in your own backup plan, and do not put them inside releases or Git.

The bundled generic sequencing workbook demonstrates the formatter and allowed-value validation. Laboratory deployments should select their existing original workbook, preserving layout and validation rules. These examples do not change any historical template, reference, primer or project file.
