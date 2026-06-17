# Agent Workflow Rules

## Default Work Cycle (Mandatory)

For every bugfix or release-oriented update in this repository, work in this order:

1. Reproduce locally with the real app or a faithful isolated QA database.
2. Fix the smallest responsible code path.
3. Verify with the `test/` folder data and export checks.
4. Update docs and this `AGENT.md` when workflow rules change.
5. Sync local code, GitHub, release notes, and deployment/distribution artifacts.

## Required Test Folder Scenario

The `test/` folder is the canonical smoke test for cloning primer design:

- `pGreenII-62SK.dna` is the empty/vector backbone.
- Use BamHI / EcoRI restriction sites.
- Expected 21 bp homology arms:
  - `HOMO_F=CGCTCTAGAACTAGTGGATCC`
  - `HOMO_R=GATAAGCTTGATATCGAATTC`
- `ZmMYBR26.dna` and `ZmMYBS1.dna` are CDS inserts.
- The required test mode is CDS with terminal stop codon retained (`remove_stop=False`).
- After each relevant change, run import -> backbone setup -> CDS selection -> primer generation -> save history -> Sangon Excel export.

## User-Facing Error Rule

Do not expose raw exceptions alone. Primer-design failures must include Chinese, actionable parameter advice. At minimum, mention the relevant knobs when applicable:

- `Tm min`, `Tm max`, `Target Tm`
- `min primer length`, `max primer length`
- `GC min`, `GC max`
- homology arm length and direction
- whether to remove or retain stop codon

## Basic / Advanced UI Rule

Basic mode is the default. Safe defaults can be simplified, but parameters that change experimental results must remain visible:

- stop codon removal/retention
- `Tm min`, `Tm max`, `Target Tm`
- primer length range
- homology arm length

Advanced mode may fold GC limits, ion/Tm profile details, frame-check settings, manual homology-arm sequences, and restriction-site details.

If a user uploads a vector map and restriction sites are unique, Basic mode should allow selecting sites and homology-arm length, then auto-generate homology arms. If sites are not unique or cannot be found, show a clear manual-entry path in Advanced mode.

## Sangon Excel Export Rule

Every Sangon Excel export must contain only the primers from the current run/export request. It must not retain old primer names, sequences, remarks, row numbers, or stale order rows from the built-in template or a user-uploaded template.

Every Sangon Excel export must reset order header cells `B4:B14` to `0` before writing current-user Profile fields. Blank Profile fields must stay `0`; never inherit personal header data from the built-in template or from a user-uploaded template.

When touching export code, verify two consecutive exports with different primer rows; the second file must not contain the first file's content.

## Construct Map Export Rule

Map export must prioritize GenBank `.gb/.gbk` as the stable source format. Do not hard-depend on native SnapGene `.dna` generation.

If uploaded backbone features can be read from GenBank or SnapGene `.dna`, preserve them as much as possible in the exported construct map.

The inserted fragment must be annotated according to the user's selected feature type: `CDS`, `promoter`, or `misc_feature`.

Forward/reverse homology arms and restriction sites must be written as features with clear labels such as `F_homology_arm`, `R_homology_arm`, `BamHI`, and `EcoRI`.

Native SnapGene `.dna` export may be offered only as an optional conversion step through a verified local SnapGene Desktop CLI. If the CLI is absent or conversion fails, keep `.gbk` export available and show a Chinese action hint.

After every construct-map export fix, validate with the `test/` folder pGreenII-62SK + BamHI/EcoRI + ZmMYBR26/ZmMYBS1 scenario.

## Branch and Release Rule (Mandatory)

For this repository, any software update intended for release must end on `main`.

Required workflow:

1. Development can happen on feature branches.
2. Before delivery/release, changes must be merged into `main`.
3. `main` must be pushed to `origin/main`.
4. Release tags and GitHub Releases must be created from `main`.

## Additional Notes

- Do not use `master` for new updates.
- If a change is completed on a feature branch, always sync it to `main` in the same update cycle.
