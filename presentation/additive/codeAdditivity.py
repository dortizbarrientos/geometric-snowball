import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib import patheffects

# --- Figure setup (square, clean, "publication-ish") ---
fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.set_aspect('equal')
ax.axis('off')

# --- Node layout: 8 contributors around, 1 outcome in centre ---
labels = list("ABCDEFGH")
n = len(labels)
angles = np.linspace(np.deg2rad(150), np.deg2rad(150) - 2*np.pi, n, endpoint=False)
R = 2.55
pos = {lab: (R*np.cos(a), R*np.sin(a)) for lab, a in zip(labels, angles)}
pos["Y"] = (0.0, 0.0)  # outcome

# --- Styling choices (kept simple, cool-toned, strong contrast for labels) ---
# Use matplotlib's default colour cycle, but "cool it down" by blending towards blue.
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
def coolen(hexcol, w=0.35):
    # blend towards a dark blue
    import matplotlib.colors as mcolors
    c = np.array(mcolors.to_rgb(hexcol))
    target = np.array(mcolors.to_rgb("#1f4e79"))  # dark blue-ish
    return tuple((1-w)*c + w*target)

node_face = {lab: coolen(cycle[i % len(cycle)], w=0.55) for i, lab in enumerate(labels)}
node_face["Y"] = coolen("#9467bd", w=0.35)  # purple-ish for outcome

edge_col = "#d7866b"  # muted orange (like the example), used sparingly

# --- Draw edges: each contributor points to outcome (no contributor–contributor edges) ---
for i, lab in enumerate(labels):
    x0, y0 = pos[lab]
    x1, y1 = pos["Y"]
    # small alternating curvature to avoid perfect overlap
    rad = 0.18 if i % 2 == 0 else -0.18
    arr = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle='-',
        mutation_scale=12,
        lw=2.2,
        color=edge_col,
        connectionstyle=f"arc3,rad={rad}",
        alpha=0.95,
        capstyle='round',
        joinstyle='round'
    )
    ax.add_patch(arr)

# --- Draw nodes ---
def draw_node(x, y, r, face, text, text_size=22):
    circ = Circle((x, y), r, facecolor=face, edgecolor="#243447", lw=3.2)
    ax.add_patch(circ)
    t = ax.text(x, y, text, ha='center', va='center',
                color="white", fontsize=text_size, fontweight='bold', family='sans-serif')
    t.set_path_effects([patheffects.withStroke(linewidth=1.5, foreground=(0,0,0,0.15))])

for lab in labels:
    draw_node(*pos[lab], r=0.52, face=node_face[lab], text=lab, text_size=22)

draw_node(*pos["Y"], r=0.60, face=node_face["Y"], text="Trait\n(Y)", text_size=16)

# --- Add simple "additivity" cue: sum formula + no-interaction statement ---
ax.text(0, -3.15, r"$Y \;=\; A + B + C + D + E + F + G + H$",
        ha='center', va='center', fontsize=18, color="#36454f")

ax.text(0, -3.55, "0 interactions (each part adds independently)",
        ha='center', va='center', fontsize=16, color=edge_col, style='italic')

# --- Tight limits with a bit of padding ---
ax.set_xlim(-3.6, 3.6)
ax.set_ylim(-3.9, 3.6)

# Save outputs
png_path = "/mnt/data/additivity_diagram.png"
pdf_path = "/mnt/data/additivity_diagram.pdf"
plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
png_path, pdf_path
