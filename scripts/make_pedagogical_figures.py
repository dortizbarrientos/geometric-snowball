from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Arc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'

# A restrained, colour-blind-safe palette.
INK = '#222222'
GREY = '#7A7A7A'
LIGHT = '#E9ECEF'
BLUE = '#0072B2'
SKY = '#56B4E9'
ORANGE = '#E69F00'
VERM = '#D55E00'
GREEN = '#009E73'
PURPLE = '#CC79A7'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.linewidth': 0.7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


def clean(ax):
    ax.set_aspect('equal', adjustable='box')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=0.6, length=3)


def ellipse_from_cov(ax, cov, center=(0, 0), nstd=1.0, **kwargs):
    vals, vecs = np.linalg.eigh(0.5 * (cov + cov.T))
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = math.degrees(math.atan2(vecs[1, 0], vecs[0, 0]))
    e = Ellipse(center, 2*nstd*math.sqrt(vals[0]), 2*nstd*math.sqrt(vals[1]), angle=angle, **kwargs)
    ax.add_patch(e)
    return vals, vecs


def fig_G_H_S():
    # Three aligned panels in one pedagogical plate.
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.55))

    G = np.array([[2.2, 1.0], [1.0, 1.0]])
    H = np.array([[1.8, -0.75], [-0.75, 1.25]])
    valsG, vecsG = np.linalg.eigh(G)
    Gh = vecsG @ np.diag(np.sqrt(valsG)) @ vecsG.T
    S = Gh @ H @ Gh

    # G panel
    ax = axes[0]
    clean(ax)
    vals, vecs = ellipse_from_cov(ax, G, fc=SKY, ec=BLUE, alpha=0.30, lw=1.5)
    for j, lab in [(0, r'$\mathbf{g}_{\max}$'), (1, r'$\mathbf{g}_2$')]:
        u = vecs[:, j]
        scale = 1.35 * math.sqrt(vals[j])
        ax.add_patch(FancyArrowPatch((0, 0), tuple(scale*u), arrowstyle='-|>', mutation_scale=10,
                                     lw=1.3, color=BLUE if j == 0 else GREY))
        ax.text(*(1.08*scale*u), lab, color=BLUE if j == 0 else GREY, ha='center', va='center')
    ax.set_xlim(-2.4, 2.4); ax.set_ylim(-2.2, 2.2)
    ax.set_xlabel('trait 1'); ax.set_ylabel('trait 2')
    ax.set_title(r'Genetic covariance $\mathbf{G}$')
    ax.text(-2.25, 1.95, 'Variance ellipse:\nlong axes contain more\nadditive variation', va='top', color=INK)

    # H panel: draw equal fitness-loss contours d'Hd = c
    ax = axes[1]
    clean(ax)
    Hinv = np.linalg.inv(H)
    for nstd, alpha in [(0.75, 0.18), (1.15, 0.13), (1.55, 0.09)]:
        ellipse_from_cov(ax, Hinv, nstd=nstd, fc=ORANGE, ec=VERM, alpha=alpha, lw=1.2)
    d = np.array([1.0, 0.85])
    grad = -H @ d
    ax.plot(*d, 'o', color=INK, ms=4)
    ax.add_patch(FancyArrowPatch(tuple(d), tuple(d + 0.85*grad), arrowstyle='-|>', mutation_scale=11,
                                 lw=1.6, color=VERM))
    ax.text(*(d + 0.47*grad + np.array([0.12, 0.08])), r'$\boldsymbol{\beta}=-\mathbf{H}\mathbf{d}$', color=VERM)
    ax.plot(0, 0, marker='*', color=GREEN, ms=9)
    ax.text(0.08, -0.18, 'optimum', color=GREEN)
    ax.set_xlim(-2.4, 2.4); ax.set_ylim(-2.2, 2.2)
    ax.set_xlabel('trait 1'); ax.set_ylabel('trait 2')
    ax.set_title(r'Fitness curvature $\mathbf{H}$')
    ax.text(-2.25, 1.95, 'Fitness-loss contours:\ntight directions are\nmore strongly penalised', va='top', color=INK)

    # S panel in whitened coordinates
    ax = axes[2]
    clean(ax)
    vals, vecs = np.linalg.eigh(S)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    # unit circle = whitened G
    ax.add_patch(Ellipse((0, 0), 2.15, 2.15, fc=SKY, ec=BLUE, alpha=0.10, lw=1.0, ls='--'))
    # represent S as rate/curvature ellipse using covariance-like scaling for display
    display_cov = vecs @ np.diag(vals / vals.max() * np.array([2.1, 0.75])) @ vecs.T
    ellipse_from_cov(ax, display_cov, fc=PURPLE, ec=PURPLE, alpha=0.22, lw=1.5)
    for j in range(2):
        u = vecs[:, j]
        scale = 1.45 if j == 0 else 1.05
        ax.add_patch(FancyArrowPatch((0, 0), tuple(scale*u), arrowstyle='-|>', mutation_scale=11,
                                     lw=1.6, color=PURPLE if j == 0 else GREY))
    ax.text(*(1.55*vecs[:, 0]), 'leading channel', color=PURPLE, ha='center', va='center')
    ax.set_xlim(-2.4, 2.4); ax.set_ylim(-2.2, 2.2)
    ax.set_xlabel('genetically whitened trait 1'); ax.set_ylabel('genetically whitened trait 2')
    ax.set_title(r'Joint operator $\mathbf{S}=\mathbf{G}^{1/2}\mathbf{H}\mathbf{G}^{1/2}$')
    ax.text(-2.25, 1.95, 'Whitening makes genetic\nvariance spherical; $\\mathbf{S}$ then\nreveals joint response modes', va='top', color=INK)

    fig.tight_layout(w_pad=2.2)
    path = OUT / 'pedagogy_G_H_S.pdf'
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def fig_principal_angles():
    fig = plt.figure(figsize=(9.8, 3.7))
    ax1 = fig.add_axes([0.05, 0.14, 0.27, 0.75])
    ax2 = fig.add_axes([0.37, 0.14, 0.27, 0.75])
    ax3 = fig.add_axes([0.69, 0.14, 0.27, 0.75])

    # 1D angle
    for ax in (ax1, ax2, ax3):
        ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.15, 1.25)
    th = math.radians(38)
    u1 = np.array([1., 0.]); u2 = np.array([math.cos(th), math.sin(th)])
    ax1.plot([-u1[0], u1[0]], [-u1[1], u1[1]], color=BLUE, lw=2)
    ax1.plot([-u2[0], u2[0]], [-u2[1], u2[1]], color=ORANGE, lw=2)
    ax1.add_patch(Arc((0,0), 0.8, 0.8, theta1=0, theta2=38, color=INK, lw=1.2))
    ax1.text(0.47, 0.14, r'$\psi$', color=INK)
    ax1.text(0, 1.12, 'Principal angle', ha='center', weight='bold')
    ax1.text(0, -1.06, r'$\cos\psi=|\mathbf{u}_1^\top\mathbf{u}_2|$', ha='center')

    # Projection length
    p = np.dot(u1, u2)*u1
    ax2.add_patch(FancyArrowPatch((0,0), tuple(u2), arrowstyle='-|>', mutation_scale=12, lw=2, color=ORANGE))
    ax2.add_patch(FancyArrowPatch((0,0), tuple(p), arrowstyle='-|>', mutation_scale=12, lw=2, color=BLUE))
    ax2.plot([u2[0], p[0]], [u2[1], p[1]], color=GREY, ls='--', lw=1.2)
    ax2.plot([-1.1,1.1],[0,0], color=BLUE, lw=1.2, alpha=0.5)
    ax2.text(p[0]/2, -0.13, r'length $=\cos\psi$', ha='center', color=BLUE)
    ax2.text(0, 1.12, 'Projection interpretation', ha='center', weight='bold')
    ax2.text(0, -1.06, 'Alignment is the fraction of a unit\ndirection retained after projection', ha='center')

    # k-volume interpretation using a slanted parallelogram projected to baseline
    base = np.array([[-0.85,-0.40], [0.25,-0.40], [0.78,0.42], [-0.32,0.42]])
    proj = base.copy(); proj[:,1] *= 0.42
    ax3.fill(base[:,0], base[:,1], fc=ORANGE, ec=VERM, alpha=0.23, lw=1.5)
    ax3.fill(proj[:,0], proj[:,1]-0.55, fc=SKY, ec=BLUE, alpha=0.30, lw=1.5)
    for i in range(4):
        ax3.plot([base[i,0], proj[i,0]], [base[i,1], proj[i,1]-0.55], color=GREY, ls=':', lw=0.8)
    ax3.text(0, 1.12, '$k$-volume interpretation', ha='center', weight='bold')
    ax3.text(0, -1.06, r'volume ratio $=\prod_j\cos\psi_j$'+'\n'+r'deficit $=1-\prod_j\cos^2\psi_j$', ha='center')

    path = OUT / 'pedagogy_principal_angles.pdf'
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def fig_hybrid_load():
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.7))
    theta = math.radians(52)
    v1 = np.array([math.cos(theta/2), math.sin(theta/2)])
    v2 = np.array([math.cos(theta/2), -math.sin(theta/2)])
    def projector(v): return np.outer(v, v)
    dg, dh = 2.1, 1.5
    gp, hp = 0.35, 2.0
    G1 = gp*np.eye(2) + dg*projector(v1)
    G2 = gp*np.eye(2) + dg*projector(v2)
    H1 = hp*np.eye(2) - dh*projector(v1)
    H2 = hp*np.eye(2) - dh*projector(v2)
    PH = 0.5*(G1+G2)

    objects = [(G1,H1,v1,'Population 1'), (G2,H2,v2,'Population 2'), (PH,0.5*(H1+H2),np.array([1.,0.]),'Recombinant hybrid')]
    for ax,(G,H,v,title) in zip(axes,objects):
        clean(ax)
        # fitness contours from H^{-1}
        Hinv=np.linalg.inv(H)
        for s,a in [(0.8,0.10),(1.25,0.08),(1.7,0.06)]:
            ellipse_from_cov(ax,Hinv,nstd=s,fc=ORANGE,ec=VERM,alpha=a,lw=1.0)
        ellipse_from_cov(ax,G,nstd=0.75,fc=SKY if 'Hybrid' not in title else PURPLE,
                         ec=BLUE if 'Hybrid' not in title else PURPLE,alpha=0.32,lw=1.6)
        ax.add_patch(FancyArrowPatch(tuple(-1.3*v),tuple(1.3*v),arrowstyle='-|>',mutation_scale=11,
                                     lw=1.5,color=GREEN if 'Hybrid' not in title else GREY))
        ax.set_xlim(-2.2,2.2); ax.set_ylim(-2.0,2.0)
        ax.set_xlabel('trait 1'); ax.set_ylabel('trait 2')
        ax.set_title(title)
    axes[0].text(-2.05,1.78,'Genetic variance lies along\na shallow local ridge',va='top')
    axes[1].text(-2.05,1.78,'The second population has\na differently oriented ridge',va='top')
    axes[2].text(-2.05,1.78,'Averaging the two covariances\nplaces variance into directions\nthat are steep in parental habitats',va='top')
    fig.suptitle('Why channel rotation can become hybrid load only under covariance-curvature coadaptation', y=1.01, weight='bold')
    fig.tight_layout(w_pad=2.2)
    path=OUT/'pedagogy_hybrid_load.pdf'
    fig.savefig(path,bbox_inches='tight')
    plt.close(fig)
    return path


def fig_handoff():
    fig, ax = plt.subplots(figsize=(10.5, 3.15))
    ax.axis('off'); ax.set_xlim(0, 10.5); ax.set_ylim(0, 3.2)
    labels = [
        ('Divergent\nselection', 0.25, BLUE),
        ('Channel\nrotation', 2.15, PURPLE),
        ('Curvature-weighted\nhybrid load', 4.05, VERM),
        ('Lower effective\ngene flow', 6.10, ORANGE),
        ('Independent\nsubstitutions', 8.05, GREEN),
        ('DMI\nsnowball', 9.65, INK),
    ]
    widths=[1.45,1.45,1.75,1.55,1.45,0.95]
    for (txt,x,col),w in zip(labels,widths):
        box=FancyBboxPatch((x,1.16),w,0.86,boxstyle='round,pad=0.08,rounding_size=0.08',
                           fc=col,ec=col,alpha=0.16,lw=1.4)
        ax.add_patch(box); ax.text(x+w/2,1.59,txt,ha='center',va='center',color=INK,weight='bold')
    for i in range(len(labels)-1):
        x=labels[i][1]+widths[i]+0.05; x2=labels[i+1][1]-0.05
        ax.add_patch(FancyArrowPatch((x,1.59),(x2,1.59),arrowstyle='-|>',mutation_scale=12,lw=1.3,color=GREY))
    ax.text(0.3,2.65,'Early, ecological, potentially reversible',color=BLUE,weight='bold')
    ax.plot([0.3,7.45],[2.48,2.48],color=BLUE,lw=1.2)
    ax.text(7.95,2.65,'Later, intrinsic, increasingly persistent',color=GREEN,weight='bold')
    ax.plot([7.95,10.35],[2.48,2.48],color=GREEN,lw=1.2)
    ax.text(5.2,0.46,r'$m_{\rm eff}=m_0e^{-L_E}$',ha='center',fontsize=11,color=ORANGE)
    ax.text(8.95,0.46,r'$D=qK_1K_2,\qquad I_D=1-e^{-\eta D}$',ha='center',fontsize=11,color=GREEN)
    ax.text(0.3,0.42,r'$d_C^2(t)=t^2\|\mathcal{K}\|_F^2+O(t^3)$',fontsize=11,color=PURPLE)
    path=OUT/'pedagogy_handoff.pdf'
    fig.savefig(path,bbox_inches='tight')
    plt.close(fig)
    return path


def fig_tensor_bridge():
    fig, ax = plt.subplots(figsize=(10.6, 4.0))
    ax.axis('off'); ax.set_xlim(0,10.6); ax.set_ylim(0,4.0)
    # Draw three covariance matrices as ellipses
    centers=[(0.8,2.7),(0.8,1.85),(0.8,1.0)]
    angles=[20,-5,-30]
    cols=[BLUE,GREEN,ORANGE]
    for idx, ((cx,cy),ang,col) in enumerate(zip(centers,angles,cols), start=1):
        ax.add_patch(Ellipse((cx,cy),1.25,0.43,angle=ang,fc=col,ec=col,alpha=0.24,lw=1.4))
        ax.text(cx-0.04,cy,r'$\mathbf{G}_{%d}$' % idx,ha='center',va='center')
    ax.text(0.8,3.55,'Symmetric matrices\nin trait space',ha='center',weight='bold')
    # arrow to svec
    ax.add_patch(FancyArrowPatch((1.6,2.0),(2.45,2.0),arrowstyle='-|>',mutation_scale=13,lw=1.4,color=GREY))
    ax.text(2.02,2.26,r'$\operatorname{svec}$',ha='center',color=INK)
    # matrix-space cloud and eigentensor axis
    pts=np.array([[3.2,1.25],[3.75,2.45],[4.35,2.85]])
    for p,col in zip(pts,cols): ax.plot(*p,'o',ms=7,color=col)
    ax.plot([2.85,4.75],[0.95,3.15],color=PURPLE,lw=2)
    ax.text(3.75,3.55,'Matrix space',ha='center',weight='bold')
    ax.text(4.55,3.05,r'leading eigentensor $\mathbf{E}_1$',rotation=48,color=PURPLE,ha='center')
    ax.text(3.75,0.52,r'$c_{i1}=\langle\mathbf{G}_i,\mathbf{E}_1\rangle_F$',ha='center')
    # arrow to response operator
    ax.add_patch(FancyArrowPatch((4.9,2.0),(5.75,2.0),arrowstyle='-|>',mutation_scale=13,lw=1.4,color=GREY))
    ax.text(5.33,2.30,r'add $\mathbf{H}_i$',ha='center')
    # S matrices as ellipses
    centers2=[(6.45,2.7),(6.45,1.85),(6.45,1.0)]
    angles2=[42,10,-42]
    for idx, ((cx,cy),ang,col) in enumerate(zip(centers2,angles2,cols), start=1):
        ax.add_patch(Ellipse((cx,cy),1.25,0.40,angle=ang,fc=col,ec=col,alpha=0.22,lw=1.4))
        ax.text(cx,cy,r'$\mathbf{S}_{%d}$' % idx,ha='center',va='center')
    ax.text(6.45,3.55,r'Response operators $\mathbf{S}_i$',ha='center',weight='bold')
    # arrow to channel geometry
    ax.add_patch(FancyArrowPatch((7.25,2.0),(8.1,2.0),arrowstyle='-|>',mutation_scale=13,lw=1.4,color=GREY))
    ax.text(7.68,2.28,'eigendecompose',ha='center')
    # channel lines and angle
    origin=np.array([9.25,1.95])
    th1=math.radians(25); th2=math.radians(-38)
    for th,col in [(th1,BLUE),(th2,ORANGE)]:
        u=np.array([math.cos(th),math.sin(th)])
        ax.plot([origin[0]-0.85*u[0],origin[0]+0.85*u[0]],
                [origin[1]-0.85*u[1],origin[1]+0.85*u[1]],color=col,lw=2.2)
    ax.add_patch(Arc(origin,0.85,0.85,theta1=-38,theta2=25,color=INK,lw=1.2))
    ax.text(9.62,2.02,r'$\psi$',ha='center')
    ax.text(9.25,3.55,'Channel subspaces\nin trait geometry',ha='center',weight='bold')
    ax.text(9.25,0.52,r'$d_C^2=\sum_j\sin^2\psi_j$',ha='center')
    path=OUT/'pedagogy_tensor_bridge.pdf'
    fig.savefig(path,bbox_inches='tight')
    plt.close(fig)
    return path


def main():
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, default=ROOT / 'figures')
    args = parser.parse_args()
    OUT = args.outdir
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [fig_G_H_S(), fig_principal_angles(), fig_hybrid_load(), fig_handoff(), fig_tensor_bridge()]
    for path in paths:
        print(path)


if __name__ == '__main__':
    main()
