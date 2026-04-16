#!/usr/bin/env python3
"""
03b_table_celegans_summary.py
==============================

PURPOSE
-------
Build a publication-ready Word table (and LaTeX version) summarising the
posterior analysis of the C. elegans populations. Reports posterior
median and 95% credible interval for the channel-anisotropy ratio,
spectral gap, and channel rotation, for all four populations in the
high-salt environment.

Also computes the spectral-gap erosion (percent loss relative to the
A6140 ancestor) as a direct test of the prediction that sustained
directional selection collapses the spectral gap.

INPUT
-----
    results/tables/celegans_posterior_summary.csv

OUTPUTS
-------
    results/tables/table_celegans_summary.docx
    results/tables/table_celegans_summary.tex
"""

from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# LOAD AND RESHAPE THE SUMMARY
# ----------------------------------------------------------------------
df = pd.read_csv(TABLES / "celegans_posterior_summary.csv")

# Focus on high-salt (the selection environment), the four populations
df = df[df.environment == "NaCl"].copy()

# Pivot so that each population has one row
wanted_qs = ["ratio", "delta", "rotation_deg"]
rows = []
for pop in ["A6140", "GA150", "GA250", "GA450"]:
    sub = df[df.population == pop]
    row = {"population": pop,
           "n_measurement": int(sub.n_measurement.iloc[0])}
    for q in wanted_qs:
        r = sub[sub.quantity == q].iloc[0]
        row[q + "_median"] = r["median"]
        row[q + "_lo95"] = r.lo95
        row[q + "_hi95"] = r.hi95
    rows.append(row)

summary = pd.DataFrame(rows)

# Spectral-gap erosion vs ancestor
anc_delta = summary.loc[summary.population == "A6140", "delta_median"].iloc[0]
summary["gap_loss_pct"] = 100 * (1 - summary["delta_median"] / anc_delta)

print("Final summary table:")
print(summary.to_string(index=False))


# ----------------------------------------------------------------------
# WORD TABLE
# ----------------------------------------------------------------------
# Use python-docx directly for precise control
from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


def set_cell_bg(cell, hex_color):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_color)
    tc_pr.append(shd)


def fmt_ci(med, lo, hi, digits=2):
    return f"{med:.{digits}f} [{lo:.{digits}f}, {hi:.{digits}f}]"


doc = Document()

# Title line
p = doc.add_paragraph()
r = p.add_run("Table 1 | Posterior summary of the channel-operator geometry in ")
r.bold = True
r = p.add_run("Caenorhabditis elegans"); r.bold = True; r.italic = True
r = p.add_run(" under selection for high-salt adaptation."); r.bold = True

p = doc.add_paragraph()
p.add_run("Posterior median and 95% credible interval computed from 10,000 "
          "MCMCglmm draws per population. Spectral-gap erosion is reported "
          "relative to the A6140 ancestor. Rotation is measured between "
          "g_max (leading eigenvector of G) and the channel direction "
          "derived from the rank-1 surrogate H_tilde = beta beta^T / ||beta||^2 + eps I, "
          "where beta is the standardized selection gradient from Mallard et al.").italic = True

# Table
headers = ["Population", "n", "lambda1/lambda2", "Spectral gap delta",
           "Gap erosion (%)", "Rotation (deg)"]
t = doc.add_table(rows=1 + len(summary), cols=len(headers))
t.alignment = WD_TABLE_ALIGNMENT.CENTER

# Header row
for j, h in enumerate(headers):
    c = t.rows[0].cells[j]
    c.text = ""
    p = c.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(h)
    run.bold = True
    run.font.size = Pt(10)
    run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    set_cell_bg(c, "2E5090")
    c.vertical_alignment = WD_ALIGN_VERTICAL.CENTER

# Data rows
for i, row in summary.iterrows():
    cells = t.rows[1 + i].cells

    def fill(idx, text, align="center"):
        c = cells[idx]
        c.text = ""
        p = c.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if align == "center" else WD_ALIGN_PARAGRAPH.LEFT
        run = p.add_run(text)
        run.font.size = Pt(9)
        if i % 2 == 0:
            set_cell_bg(c, "F2F6FA")
        c.vertical_alignment = WD_ALIGN_VERTICAL.CENTER

    fill(0, row["population"])
    fill(1, str(row["n_measurement"]))
    fill(2, fmt_ci(row["ratio_median"], row["ratio_lo95"], row["ratio_hi95"], 2))
    fill(3, fmt_ci(row["delta_median"], row["delta_lo95"], row["delta_hi95"], 3))
    gap = row["gap_loss_pct"]
    fill(4, "—" if row["population"] == "A6140" else f"{gap:.0f}")
    fill(5, fmt_ci(row["rotation_deg_median"], row["rotation_deg_lo95"], row["rotation_deg_hi95"], 1))

# Set column widths
widths = [Cm(2.0), Cm(1.2), Cm(3.2), Cm(3.5), Cm(2.0), Cm(3.5)]
for row in t.rows:
    for j, c in enumerate(row.cells):
        c.width = widths[j]

out_docx = TABLES / "table_celegans_summary.docx"
doc.save(out_docx)
print(f"\nWrote {out_docx}")

# ----------------------------------------------------------------------
# LATEX VERSION (for supplementary or Overleaf)
# ----------------------------------------------------------------------
tex = r"""\begin{table}[ht]
\centering
\small
\begin{tabular}{lrllrr}
\toprule
\textbf{Population} & $n$ & $\boldsymbol{\lambda_1/\lambda_2}$ & \textbf{Gap $\boldsymbol{\delta}$} & \textbf{Erosion (\%)} & \textbf{Rotation (deg)} \\
\midrule
"""
for _, row in summary.iterrows():
    tex += (
        f"{row['population']} & {row['n_measurement']} & "
        f"{fmt_ci(row['ratio_median'], row['ratio_lo95'], row['ratio_hi95'], 2)} & "
        f"{fmt_ci(row['delta_median'], row['delta_lo95'], row['delta_hi95'], 3)} & "
    )
    tex += ("---" if row["population"] == "A6140" else f"{row['gap_loss_pct']:.0f}")
    tex += (f" & {fmt_ci(row['rotation_deg_median'], row['rotation_deg_lo95'], row['rotation_deg_hi95'], 1)} \\\\\n")
tex += r"""\bottomrule
\end{tabular}
\caption{\textbf{Posterior summary of the channel-operator geometry in
\textit{Caenorhabditis elegans}.} Posterior median and 95\% credible
interval from 10,000 MCMCglmm draws. Erosion reports percent loss of
the spectral gap relative to the A6140 ancestor.}
\label{tab:celegans_summary}
\end{table}
"""
out_tex = TABLES / "table_celegans_summary.tex"
out_tex.write_text(tex)
print(f"Wrote {out_tex}")
