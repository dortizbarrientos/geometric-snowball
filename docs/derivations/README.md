# Derivation documents

`full_derivations.tex` is the top-level source of the self-contained teaching monograph; its chapter bodies are maintained under `chapters/` so mathematical and biological revisions can be reviewed in focused diffs. It develops the linear algebra, perturbation theory, hybrid-fitness bridge, standing-variation dynamics, and early-to-late speciation argument, with verbal, geometric, and biological interpretations following the mathematical steps.

`technical_note.tex` is a compact mathematical development useful for orienting the software implementation.

The generated PDFs and supporting figures are intentionally not tracked in the clean v2 source tree. Rebuild the full set with:

```bash
make docs
```

For a direct LaTeX-only rebuild after the figures and audit have been generated:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error docs/derivations/full_derivations.tex
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error docs/derivations/technical_note.tex
```

The monograph is the conceptual source of truth. The tested package is the computational source of truth. When they disagree, open an issue and add a numerical identity test before changing either one.
