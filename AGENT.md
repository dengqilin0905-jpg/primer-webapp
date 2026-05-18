# Agent Workflow Rules

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

