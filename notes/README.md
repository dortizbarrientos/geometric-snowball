# Notes

Companion learning materials for mathematical work related to *Additive Channels in curved spaces*.

> **Additive Channels in Curved Fitness Landscapes**
> Ortiz-Barrientos & Cooper · submitted to *PNAS* · The University of Queensland

> **Additive Channels in Curved Fitness Landscapes**
> Ortiz-Barrientos & Cooper · submitted to *PNAS* · The University of Queensland

---

## Structure

## Structure
```
notes/
├── prereqs/        Foundational mathematics (before the paper)
├── apps/           Interactive topic apps (one per key result)
│   └── METexample/ Introductory apps — what is a Multi-Environment Trial?
├── github/         Git + Overleaf workflow guide
└── README.md       This file
```

---

## prereqs/

Four tutorials covering the mathematical background assumed by the paper.
Work through these before the apps.

| File | Topic |
|---|---|
| `math_tutorial.html` | Overview of all prerequisite concepts |
| `tutorial1_gradient.html` | What is a gradient? |
| `tutorial2_matrix.html` | What is the G matrix? |
| `tutorial3_5.html` | Variance, Taylor expansion, eigenvalues |

---

## apps/

Interactive learning apps — one per key result in the paper.
Each app has four learning lanes: Plain English · Geometry · Biology · Mathematics.
Live sliders connect directly to the equations.

| File | Topic | Key result |
|---|---|---|
| `learning_apps.html` | All five topics, single tabbed file | — |
| `app1_additivity_index.html` | The additivity index | Ag = V_lin / (V_lin + V_quad) |
| `app2_rate_diagnostic.html` | The rate diagnostic | dAg/dt = Ag(1 - Ag) r(t) |
| `app3_moving_optimum.html` | Moving optimum closed form | Ag = d^2 / (d^2 + G/2) |
| `app4_gradient_rotation.html` | Gradient rotation — Act III exit | V_lin(θ) = ‖b‖² (λ₁ cos²θ + λ₂ sin²θ) |
| `app5_introgression.html` | Introgression criterion | δAg > 0  iff  (b·v)² / V_lin  >  v'HGHv / V_quad |

All apps point to the companion Mathematica notebook:
`manuscript/mathematica/AdditiveChannels_Expanded_Annotated.nb`

### apps/METexample/

Introductory series for readers unfamiliar with Multi-Environment Trials.
Work through these before the main apps if the MET context is new.

| File | Topic |
|---|---|
| `app0_what_is_a_MET.html` | What is a Multi-Environment Trial? |
| `app1_genetic_building_blocks_MET.html` | Genetic building blocks |
| `app2_reaction_norm_model_MET.html` | The reaction norm model |
| `app3_cross_env_covariance_MET.html` | Cross-environment genetic covariance |
| `app4_linear_vs_quad_GEI_MET.html` | Linear vs quadratic GEI |
| `app5_reading_diagnostics_MET.html` | Reading the diagnostics |

---

## github/

Workflow documentation for the Git + Overleaf pipeline used during manuscript preparation.

| File | Format | Contents |
|---|---|---|
| `github_tutorial_academics.html` | Interactive app | Step-by-step GitHub guide for academics |
| `github_tutorial_academics.md` | Markdown | Plain-text version of the same guide |
| `github_tutorial_academics.pptx` | Slides | Presentation version (10 slides) |

The file `git_overleaf_pipeline.md` in the root of `notes/` is gitignored —
it contains local paths specific to this machine and is not tracked.

---

## How to use

Open any HTML file directly in a browser — no server required.
Start with `prereqs/` if the mathematics is unfamiliar.
If the MET context is new, start with `apps/METexample/` before the main apps.
Move to `apps/` for the key results from the paper.
Equations render via MathJax (internet connection required on first load).
