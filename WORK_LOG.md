# Software release work log

The software release is built from the public main baseline using an allowlisted, sanitized source snapshot. Private development history is not imported. Laboratory worktrees and datasets are preserved separately.

Resources are explicit external configuration: presets/maps, project registrations, identities, annotations, original provider workbooks and optional workflow rules. Public defaults are empty registries and generic workbooks; examples and tests are artificial. A private regression suite preserves the existing scientific assertions outside Git.

Publication checks inspect both the current software tree and new reachable Git objects. CI never deploys. Installation, public synthetic tests, private scientific regression, real-case comparisons and protected-source hashes are verified before delivery. Detailed private evidence stays outside this repository.
