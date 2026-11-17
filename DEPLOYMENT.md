# Software publication and private deployment are separate

This repository does not update a production server on push or merge. The GitHub workflow only runs public software tests and publication checks; it contains no deployment secrets. scripts/remote_deploy.sh intentionally exits without stopping services, copying data or replacing directories.

A private deployment should install an explicitly reviewed commit/wheel into a new `releases/<SHA>` directory and a fresh environment. Keep the database, verification ledger, formal projects, experimental maps, templates, secrets and service configuration outside release directories. Point a separately maintained launcher to the selected release only after isolated acceptance and human approval. The example service file is a placeholder and is never automatically installed.

Before switching, use an isolated database and ledger for preview, check synthetic workflows and your private integration suite, and make an SQLite-consistent backup with the SQLite backup API. Do not copy a live SQLite file blindly. Record the old release and environment. Roll back the launcher to them if necessary; do not overwrite newly created database or sequencing records as part of a code rollback.

Never run a directory-replacement deployment script against a root containing laboratory data. Never run git pull/reset in a dirty laboratory working tree. Main protection, PR review and required checks should be configured by the repository administrator before publication. If the remote repository still has an older automatic deployment workflow, disable it before merging the software-only workflow.
