# GitHub for Academics
## A Complete Guide — From Zero to Versioned Science

*Illustrated with the `additive-channels` project (Ortiz-Barrientos & Cooper, PNAS)*

---

## Quick Reference Card

```bash
# THE FOUR COMMANDS YOU NEED EVERY DAY
git pull overleaf master          # get latest from Overleaf / collaborator
git add .                         # stage all changes
git commit -m "What you changed"  # save a snapshot
git push origin main              # send to GitHub
```

```bash
# WHEN THINGS GO WRONG
git status                        # what is going on?
git log --oneline                 # see history
git pull origin main --rebase     # sync before pushing
git diff                          # see exact line changes
```

---

## Part 1 — Why GitHub?

### The problem every researcher recognises

```
thesis_final.docx
thesis_final_v2.docx
thesis_final_ACTUAL_final.docx
thesis_final_mark_comments.docx
thesis_final_mark_comments_revised.docx
thesis_SUBMIT_THIS_ONE.docx
```

Version control replaces this. Git is the tool. GitHub is the platform.

### What Git gives you

- **A complete history** of every change ever made, with who made it and when
- **The ability to undo anything** — go back to any previous state
- **Safe collaboration** — multiple people editing without overwriting each other
- **Proof of priority** — timestamped commits establish when ideas were developed
- **A backup** that is not on your laptop or your institution's drive

### What GitHub adds

- A web interface to browse your project history
- Issue tracking for to-do lists and bug reports
- A permanent URL to share with reviewers, editors, and the public
- Integration with Zenodo for a citable DOI
- Integration with Overleaf for LaTeX manuscripts

---

## Part 2 — Core concepts

### The four objects you need to understand

| Object | What it is | Analogy |
|--------|-----------|---------|
| **Repository (repo)** | The project folder git tracks | A lab notebook |
| **Commit** | A saved snapshot of the project | A dated entry in the notebook |
| **Branch** | A parallel version for experimenting | A scratch page |
| **Remote** | A copy of the repo on another server | A photocopy at the library |

### The three states of a file

```
Working directory  →  Staging area  →  Repository
  (you edit here)     (git add .)      (git commit)
```

A file is either:
- **Untracked** — git does not know about it
- **Modified** — git knows about it and it has changed
- **Staged** — ready to be committed
- **Committed** — safely saved in history

### Remotes

You can have multiple remotes — named copies of the repo on different servers.

```
origin    →  GitHub   (your archive and backup)
overleaf  →  Overleaf (your writing environment)
```

This is the setup used in the `additive-channels` project.

---

## Part 3 — Setup pipeline (complete, from scratch)

### Step 1 — Install git

```bash
# Check if git is installed
git --version

# If not, install via Homebrew (Mac)
brew install git

# Or download from https://git-scm.com
```

### Step 2 — Configure git (once, ever)

```bash
git config --global user.name "Your Name"
git config --global user.email "your@email.com"
git config --global core.editor "nano"   # or vim, code, etc.
```

### Step 3 — Create a GitHub repository

1. Go to `https://github.com/new`
2. Fill in:
   - **Repository name** — short, lowercase, hyphens not underscores (e.g. `additive-channels`)
   - **Visibility** — Private until you post a preprint; then Public
   - **Add README** — Yes
   - **Add .gitignore** — choose TeX template for LaTeX projects
3. Click **Create repository**

### Step 4 — Clone it to your computer

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
cd YOUR_REPO
```

### Step 5 — Create folder structure

```bash
mkdir -p manuscript/figures manuscript/mathematica notes
```

### Step 6 — Copy your files in and make the first commit

```bash
cp /path/to/your/files/* manuscript/

git add .
git status                          # always check before committing
git commit -m "Initial commit: manuscript, SI, figures"
git push origin main
```

### Step 7 — Link to Overleaf (optional but recommended)

```bash
# Add Overleaf as a second remote
git remote add overleaf https://git@git.overleaf.com/YOUR_PROJECT_ID

# Pull Overleaf content in, letting Overleaf win any conflicts
git pull overleaf master --allow-unrelated-histories -X theirs

# Push the merged result to GitHub
git push origin main
```

Your Overleaf project ID is in the URL:
`https://www.overleaf.com/project/69c719779a16bf9458552d24`
                                     ↑ this part

---

## Part 4 — Daily workflow

### The three-command cycle

```bash
git pull overleaf master          # sync from Overleaf
git add .
git commit -m "Describe change"
git push origin main              # archive to GitHub
```

Run this at the end of every writing session.

### Writing good commit messages

A good commit message completes the sentence: *"If applied, this commit will..."*

```bash
# Good
git commit -m "Revise Discussion paragraph 2"
git commit -m "Fix bib key CooperMessina2023"
git commit -m "Update Fig 4 caption and SI Section S6"
git commit -m "Add breeding implications paragraph"
git commit -m "Revise abstract after reviewer comments"

# Bad
git commit -m "update"
git commit -m "fix"
git commit -m "changes"
git commit -m "asdfgh"
```

### Checking what changed

```bash
git status                # summary of what is staged and what is not
git diff                  # exact line-by-line changes (unstaged)
git diff --staged         # exact line-by-line changes (staged)
git log --oneline         # compact history
git log --oneline -10     # last 10 commits only
```

---

## Part 5 — Collaboration

### Adding a collaborator

1. GitHub → Your repo → **Settings** → **Collaborators**
2. Add by GitHub username
3. They clone the repo and use the same workflow

### The collaboration cycle

```bash
# Before you start working each day
git pull origin main

# After you finish
git add .
git commit -m "What you did"
git push origin main
```

### What to do if two people edited the same file

Git will report a **merge conflict**. Open the file and you will see:

```
<<<<<<< HEAD
Your version of the line
=======
Their version of the line
>>>>>>> origin/main
```

Delete the markers, keep the correct version, then:

```bash
git add .
git commit -m "Resolve merge conflict in manuscript"
git push origin main
```

---

## Part 6 — Branches (for major revisions)

Branches let you work on a big change without touching the main version.

```bash
# Create a new branch for major revisions
git checkout -b revision-round-1

# Work on the manuscript...
git add .
git commit -m "Address reviewer 1 comments"

# When done, merge back to main
git checkout main
git merge revision-round-1

# Push to GitHub
git push origin main
```

Useful branch names for academic work:
```
revision-round-1
revision-round-2
submitted-version
accepted-version
```

---

## Part 7 — Tagging versions

Tags mark important moments in the project history — like submission dates.

```bash
# Tag the bioRxiv submission
git tag -a v1.0 -m "bioRxiv submission — March 2026"
git push origin v1.0

# Tag after journal submission
git tag -a v1.1 -m "PNAS submission"
git push origin v1.1

# Tag the accepted version
git tag -a v2.0 -m "Accepted — PNAS"
git push origin v2.0
```

View all tags:
```bash
git tag
```

---

## Part 8 — What to track and what to ignore

### Always track
- `.tex` source files
- `.bib` bibliography files
- `.py` figure scripts
- `.R` analysis scripts
- `.md` documentation
- Figure PDFs in `manuscript/figures/`

### Never track
- Compiled `.pdf` documents (generated from source)
- LaTeX artefacts (`.aux`, `.log`, `.bbl`, `.blg`, `.synctex.gz`)
- Large data files (use a data repository instead)
- Personal paths, tokens, or credentials

### The .gitignore file

```gitignore
# LaTeX build artefacts
*.aux
*.log
*.bbl
*.blg
*.out
*.toc
*.fls
*.fdb_latexmk
*.synctex.gz

# Compiled PDFs (tracked explicitly in figures/ only)
*.pdf
!manuscript/figures/*.pdf

# Python
__pycache__/
*.pyc

# macOS
.DS_Store

# Personal notes with sensitive paths
notes/git_overleaf_pipeline.md
```

---

## Part 9 — Going public and getting a DOI

### When to make the repo public

Make it public **at the same moment** you post the preprint. The two should appear together so anyone reading the preprint can immediately access the code.

```
Preprint posted → Repo made public → Zenodo DOI minted → DOI in paper
```

### Making the repo public

GitHub → Settings → Danger Zone → **Change visibility** → Public

### Getting a permanent DOI via Zenodo

1. Go to `https://zenodo.org`
2. Log in with your GitHub account
3. Go to **GitHub** tab → enable your repo
4. Back in GitHub, create a release:
   ```bash
   git tag -a v1.0 -m "bioRxiv submission"
   git push origin v1.0
   ```
5. GitHub → Releases → **Create release** from tag v1.0
6. Zenodo automatically mints a DOI

Replace the placeholder in your manuscript:
```
[repository DOI to be inserted before submission]
↓
https://doi.org/10.5281/zenodo.XXXXXXX
```

---

## Part 10 — Troubleshooting

### "Updates were rejected because the remote contains work you do not have"

```bash
git pull origin main --rebase
git push origin main
```

### "Nothing to commit, working tree clean" but you expect changes

```bash
git status         # confirm which directory you are in
ls                 # check the files are actually there
```

### You committed something you should not have

```bash
# Undo last commit, keep changes in working directory
git reset --soft HEAD~1

# Or if you already pushed — create a new commit that reverses it
git revert HEAD
git push origin main
```

### You want to go back to a previous version of one file

```bash
git log --oneline manuscript/main.tex    # find the commit hash
git checkout abc1234 -- manuscript/main.tex
```

### Merge conflict you cannot resolve

```bash
# Accept entirely one side or the other
git checkout --theirs manuscript/main.tex   # take their version
git checkout --ours manuscript/main.tex     # keep your version
git add manuscript/main.tex
git commit -m "Resolve conflict"
```

### "Your token has expired" (Overleaf)

Go to Overleaf → Account → Git Integration → **Add another token**. Use the new token as your password.

---

## Appendix — Repository structure for a PNAS paper

```
additive-channels/
├── README.md                          Project description and instructions
├── LICENSE                            MIT for code; manuscript © authors
├── .gitignore                         What git should ignore
│
├── manuscript/
│   ├── main.tex                       Main text LaTeX source
│   ├── si.tex                         Supplementary Information
│   ├── references.bib                 BibTeX bibliography
│   ├── coverLetter.docx               Cover letter
│   │
│   ├── figures/
│   │   ├── figures_script.py          Single script generating all figures
│   │   ├── figure1_*.pdf              Figure PDFs
│   │   └── figS1_*.pdf                SI figure PDFs
│   │
│   └── mathematica/
│       ├── notebook.nb                Mathematica notebook
│       └── notebook.wl                Wolfram Language script
│
└── notes/                             Local notes (gitignored if sensitive)
```

---

## Remotes summary

| Remote | URL | Purpose |
|--------|-----|---------|
| `origin` | `https://github.com/USERNAME/REPO.git` | GitHub archive |
| `overleaf` | `https://git@git.overleaf.com/PROJECT_ID` | Overleaf writing env |

---

*This guide uses the `additive-channels` repository as its running example.*
*Repository: `https://github.com/dortizbarrientos/additive-channels`*
