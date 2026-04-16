# Setting up the GitHub repository

A one-time setup script to initialise the repo, configure Git LFS for
the two large `.RData` files, and push to GitHub.

## Prerequisites

- Git installed (`git --version`)
- Git LFS installed (`git lfs --version` — install from https://git-lfs.com if missing)
- A GitHub repository already created (empty) at your account. For
  Daniel: https://github.com/dortizbarrientos/celegans-channel-operator

## Commands

From inside the `celegans_pipeline/` directory:

```bash
# Initialise
git init
git branch -M main

# Configure Git LFS for the large RData files
# (the .gitattributes file is already in the repo)
git lfs install

# Stage and commit
git add .
git commit -m "Initial commit: reproducible pipeline for C. elegans posterior analysis"

# Add remote and push
git remote add origin https://github.com/dortizbarrientos/celegans-channel-operator.git
git push -u origin main
```

## Verifying LFS is working

After pushing, on GitHub.com open `data/raw/G_matrices_high_salt.RData`.
The page should say "Stored with Git LFS" rather than trying to render
the binary. If it says the file exceeds 100 MB, LFS did not engage —
rerun `git lfs install` and re-commit.

## If the SOP team needs to clone

```bash
git lfs install           # once per machine
git clone https://github.com/dortizbarrientos/celegans-channel-operator.git
cd celegans-channel-operator
bash run_pipeline.sh
```

The LFS-tracked files download automatically during `clone`.
