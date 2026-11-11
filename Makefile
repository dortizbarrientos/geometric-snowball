.PHONY: install test audit simulations figures senecio-figure docs check clean

install:
	python -m pip install -e ".[all]"

test:
	pytest

audit:
	python scripts/run_derivation_audit.py --outdir build/audit

simulations:
	geometric-snowball reproduce --outdir build/simulations

figures:
	python scripts/reproduce_figures.py --outdir build/figures

senecio-figure:
	python analysis/senecio/senecio_figure.py --outfile build/senecio_pair_descriptive.pdf

docs:
	python scripts/run_derivation_audit.py --outdir docs/derivations/checks
	python scripts/reproduce_figures.py --outdir docs/derivations/figures
	python scripts/make_pedagogical_figures.py --outdir docs/derivations/figures
	python analysis/senecio/senecio_figure.py --outfile docs/derivations/figures/senecio_pair_descriptive.pdf
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error docs/derivations/full_derivations.tex
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error docs/derivations/technical_note.tex

check: test audit simulations figures senecio-figure

clean:
	rm -rf build dist .pytest_cache .ruff_cache src/*.egg-info src/geometric_snowball/*.egg-info
	rm -rf docs/derivations/figures docs/derivations/checks
	rm -f docs/derivations/full_derivations.pdf docs/derivations/technical_note.pdf
	latexmk -C -cd docs/derivations/full_derivations.tex 2>/dev/null || true
	latexmk -C -cd docs/derivations/technical_note.tex 2>/dev/null || true
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
