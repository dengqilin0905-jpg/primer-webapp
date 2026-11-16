# Software development rules

This repository contains installable software, generic templates and artificial test fixtures only. Keep all laboratory projects, real sequencing results, experimental maps, user databases, physical colony records and private configuration outside the repository. Never commit them, even temporarily; review every new reachable Git object before publication.

Preserve canonical business assembly in project_service.py and Web/Batch normalization parity. Exporters format canonical inputs and do not redesign frozen primers. Preserve minimal backbone + formal species/identity naming, default three colonies, confirmed 18-bp arms and the 1499-bp sequencing profile. Unknown identity or assembly evidence must never be guessed from names. Legacy records remain immutable.

Reference validation includes the entire design insert, non-CDS segments, junctions, linker, tag, stop and required elements. Missing coverage cannot pass. A deletion requires mapped anchors and actual supported spanning read sequences. Reference confirmation, design disagreement and actual sample identity are separate. Native SnapGene map conversion is not an alignment file.

Use external resource_config settings for private presets, maps, project migration registries, identity/annotation tables and provider templates. Default startup and custom-backbone design must work without them. Do not turn private missing inputs into false PASS or silently use the synthetic demo as an experimental preset.

Run scripts/run_public_tests.py and scripts/synthetic_smoke.py in an isolated environment. Private integration tests run externally against the same code without changing scientific assertions. Test databases and reports stay outside Git. Preserve existing worktrees and source data. No push, merge, deployment or service restart without explicit task authorization. CI must not deploy on main pushes.
