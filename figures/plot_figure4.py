import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from Bio import Phylo

# ── figstyle module (inline) ────────────────────────────────────────────────

MM = 1.0 / 25.4
COLW = 183.0

C = {
    "ink":      "#242424",
    "muted":    "#6E6E6E",
    "grid":     "#D7D7D7",
    "branch":   "#B7B7B7",
    "neutral":  "#E8E8E8",
    "absent":   "#ECECEC",
    "leader":   "#9A9A9A",
    "tmi":       "#3E6FA3",
    "tmi_light": "#B9CCE2",
    "mp":        "#B85C55",
    "mp_light":  "#E2BBB7",
    "avium":       "#D19A32",
    "colombiense": "#8D6A9F",
    "timonense":   "#D19A32",
    "other_pub":   "#A9A9A9",
    "direct":    "#9AA0A6",
    "recovered": "#2F7F7A",
    "excluded":  "#B95A55",
    "primary":   "#3D8875",
    "auxiliary": "#D9A441",
    "negative":  "#E8E8E8",
    "pairA": "#3F6E9A", "pairB": "#5B8C72", "pairC": "#C28A32", "pairD": "#9A6887",
    "alarm": "#B33A3A",
    "warning": "#B68A2D",
    "intervening": "#D9D9D9",
}


def setup(base=8, mid=7, small=6):
    preferred_font = os.environ.get("NTM_FIGURE_FONT", "Arial")
    font_stack = list(dict.fromkeys([preferred_font, "Arial", "Helvetica", "DejaVu Sans"]))
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": font_stack,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.6,
        "axes.edgecolor": C["ink"], "axes.labelcolor": C["ink"],
        "text.color": C["ink"], "xtick.color": C["ink"], "ytick.color": C["ink"],
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.6, "ytick.major.size": 2.6,
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "xtick.labelsize": small, "ytick.labelsize": small,
        "legend.fontsize": small, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": None,
        "axes.titlepad": 6,
    })
    return dict(base=base, mid=mid, small=small)


def set_frame(ax, keep="lb"):
    m = {"l": "left", "r": "right", "t": "top", "b": "bottom"}
    for k, name in m.items():
        ax.spines[name].set_visible(k in keep)
    ax.tick_params(top=False, right=False)
    if "l" not in keep: ax.tick_params(left=False)
    if "b" not in keep: ax.tick_params(bottom=False)


def panel_letter(ax, letter, dx=-0.02, dy=1.02, size=11):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right", color=C["ink"])


def load_prune_ladderize(tree_path, keep_ids=None, anchor_id=None):
    tree = Phylo.read(tree_path, "newick")
    if anchor_id is not None:
        anchor = next((t for t in tree.get_terminals() if t.name == anchor_id), None)
        if anchor is not None:
            tree.root_with_outgroup(anchor)
    if keep_ids is not None:
        keep = set(keep_ids)
        for tip in list(tree.get_terminals()):
            if tip.name not in keep:
                tree.prune(tip)
    tree.ladderize(reverse=True)
    return tree


def clade_xy(tree):
    x = tree.depths()
    if not max(x.values()):
        x = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    y = {t: i for i, t in enumerate(reversed(terminals))}
    def assign(cl):
        for ch in cl.clades:
            assign(ch)
        if cl not in y:
            y[cl] = float(np.mean([y[ch] for ch in cl.clades]))
    assign(tree.root)
    return x, y


def draw_tree(ax, tree, lw=0.5, color=None):
    color = color or C["branch"]
    x, y = clade_xy(tree)
    parent_of = {}
    for cl in tree.find_clades(order="preorder"):
        for ch in cl.clades:
            parent_of[ch] = cl
    for cl in tree.find_clades(order="preorder"):
        if cl.clades:
            cy = [y[ch] for ch in cl.clades]
            ax.plot([x[cl], x[cl]], [min(cy), max(cy)], color=color, lw=lw, zorder=1,
                    solid_capstyle="round")
        if cl in parent_of:
            p = parent_of[cl]
            ax.plot([x[p], x[cl]], [y[cl], y[cl]], color=color, lw=lw, zorder=1,
                    solid_capstyle="round")
    order = sorted(tree.get_terminals(), key=lambda t: y[t])
    set_frame(ax, "b")
    ax.set_yticks([])
    ax.set_ylim(-1, len(tree.get_terminals()))
    return x, y, order


def presence_strip(ax, order, y, presence_by_tip, blocks, block_colors,
                   x0=0.0, cell_w=1.0, cell_h=0.82, gap=0.12, absent=None):
    absent = absent or C["absent"]
    for ci, block in enumerate(blocks):
        cx = x0 + ci * (cell_w + gap)
        for tip in order:
            present = bool(presence_by_tip.get(tip.name, {}).get(block, False))
            face = block_colors[block] if present else absent
            ax.add_patch(Rectangle((cx, y[tip] - cell_h/2), cell_w, cell_h,
                                    facecolor=face, edgecolor="white", linewidth=0.3, zorder=2))
    ax.set_xlim(x0 - gap, x0 + len(blocks) * (cell_w + gap))
    ax.set_ylim(-1, len(order))
    ax.axis("off")
    return [x0 + ci*(cell_w+gap) + cell_w/2 for ci in range(len(blocks))]


def dumbbell(ax, y, left_x, right_x, left_c, right_c, stem=None, s=34, lw=1.0):
    stem = stem or C["leader"]
    ax.plot([left_x, right_x], [y, y], color=stem, lw=lw, zorder=1, solid_capstyle="round")
    ax.scatter([left_x], [y], s=s, color=left_c, zorder=3, edgecolor="white", linewidth=0.4)
    ax.scatter([right_x], [y], s=s, color=right_c, zorder=3, edgecolor="white", linewidth=0.4)


def _one_gene_arrow(ax, x, y, w, h, strand, color):
    head = min(0.42 * w, 0.34)
    if strand >= 0:
        verts = [(x, y-h/2), (x+w-head, y-h/2), (x+w-head, y-h*0.8),
                 (x+w, y), (x+w-head, y+h*0.8), (x+w-head, y+h/2), (x, y+h/2), (x, y-h/2)]
    else:
        verts = [(x+w, y-h/2), (x+head, y-h/2), (x+head, y-h*0.8),
                 (x, y), (x+head, y+h*0.8), (x+head, y+h/2), (x+w, y+h/2), (x+w, y-h/2)]
    codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="white",
                           linewidth=0.3, zorder=3))


def gene_arrows(ax, y, strands, color, x0=0.0, glyph_w=1.0, gap=0.28, h=0.5,
                grey_mask=None, grey=None):
    grey = grey or C["intervening"]
    x = x0
    for i, s in enumerate(strands):
        col = grey if (grey_mask and grey_mask[i]) else color
        _one_gene_arrow(ax, x, y, glyph_w, h, s, col)
        x += glyph_w + gap
    return x


def save(fig, stem, dpi=300):
    output_dir = os.environ.get("NTM_FIGURE_OUTPUT_DIR", ".")
    os.makedirs(output_dir, exist_ok=True)
    pdf = os.path.join(output_dir, f"{stem}.pdf")
    png = os.path.join(output_dir, f"{stem}.png")
    fig.savefig(pdf, facecolor="white", bbox_inches="tight")
    fig.savefig(png, dpi=dpi, facecolor="white", bbox_inches="tight")
    return pdf, png


# ── load data ───────────────────────────────────────────────────────────────

PHY_TREE = "phylogeny/intracellulare_complex_94.treefile"

f3a = pd.read_csv("source_data/Figure4a_tree_aligned_interval_presence.tsv", sep="\t")
f3b = pd.read_csv("source_data/Figure4b_full_atlas_interval_prevalence.tsv", sep="\t")
f3c = pd.read_csv("source_data/Figure4c_bioproject_stratified_statistics.tsv", sep="\t")
f3d = pd.read_csv("source_data/Figure4d_representative_gene_order.tsv", sep="\t")

BLOCK_ORDER = [
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
]
SHORT = {
    "TMI_aromatic_catabolism_associated_block": "Aromatic",
    "TMI_methyltransferase_hydrolase_cupin_block": "Methyl./cupin",
    "MP_MIP_nitrogen_redox_associated_block": "N-redox",
    "MP_MIP_oxidoreductase_associated_block": "Oxidoreductase",
}

ids75 = set(f3a.tree_id)
lin = dict(zip(f3a.tree_id, f3a.accessory_lineage))
src = dict(zip(f3a.tree_id, f3a.source))
samp = dict(zip(f3a.tree_id, f3a.sample_id))
pres = {r.tree_id: {b: int(r[b]) for b in BLOCK_ORDER} for _, r in f3a.iterrows()}

# ── build figure ────────────────────────────────────────────────────────────

setup()
block_color = {b: (C["tmi"] if b.startswith("TMI_") else C["mp"]) for b in BLOCK_ORDER}


def lin_color(tag):
    return C["tmi"] if tag == "MI_TMI_lineage" else C["mp"]


def build_fig3_v3():
    tree = load_prune_ladderize(PHY_TREE, keep_ids=ids75, anchor_id="ANCHOR_GCF_002219285.1")
    fig = plt.figure(figsize=(COLW*MM, 220*MM))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 0.92], height_ratios=[1.06, 0.9, 1.06],
                          left=0.035, right=0.965, top=0.92, bottom=0.055, wspace=0.30, hspace=0.5)
    # Panel a: tree | label-col | heat
    gsa = gs[:, 0].subgridspec(1, 3, width_ratios=[1.0, 0.14, 0.30], wspace=0.0)
    ax_tree = fig.add_subplot(gsa[0, 0])
    ax_lab = fig.add_subplot(gsa[0, 1], sharey=ax_tree)
    ax_heat = fig.add_subplot(gsa[0, 2], sharey=ax_tree)
    x, y, order = draw_tree(ax_tree, tree, lw=0.5)
    n = len(order)
    for tip in order:
        nm = tip.name
        col = lin_color(lin.get(nm, ""))
        if src.get(nm) == "local":
            ax_tree.scatter(x[tip], y[tip], marker="*", s=34, facecolor=col, edgecolor="white", linewidth=0.4, zorder=4)
            ax_lab.text(0.08, y[tip], samp[nm], fontsize=6.0, fontweight="bold", ha="left", va="center", color=col)
        else:
            ax_tree.scatter(x[tip], y[tip], marker="o", s=6, facecolor=col, edgecolor="none", zorder=3)
    ax_lab.set_xlim(0, 1)
    ax_lab.axis("off")
    ax_lab.set_ylim(-1, n+7)
    ax_tree.set_xlabel("Substitutions per core SNP site")
    ax_tree.set_xlim(0, max(x.values())*1.02)
    ax_tree.set_title("TMI/MP-MIP complete-core SNP phylogeny (n = 75)", loc="left", pad=6, x=0.02)
    ax_tree.set_ylim(-1, n+7)
    centers = presence_strip(ax_heat, order, y, pres, BLOCK_ORDER, block_color,
                             x0=0.0, cell_w=1.0, cell_h=0.82, gap=0.18, absent=C["absent"])
    ax_heat.set_ylim(-1, n+7)
    for cx, b in zip(centers, BLOCK_ORDER):
        ax_heat.text(cx, n+0.4, SHORT[b], rotation=90, ha="center", va="bottom",
                     fontsize=6.0, color=block_color[b])
    leg_a = [
        Line2D([], [], marker="o", linestyle="", markerfacecolor=C["tmi"], markeredgecolor="none", markersize=5, label="TMI"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor=C["mp"], markeredgecolor="none", markersize=5, label="MP-MIP"),
        Line2D([], [], marker="*", linestyle="", markerfacecolor="#555", markeredgecolor="white", markersize=7, label="Local genome"),
    ]
    ax_tree.legend(handles=leg_a, loc="upper left", bbox_to_anchor=(0.02, 0.99), frameon=False,
                   handletextpad=0.4, labelspacing=0.3, fontsize=6)
    panel_letter(ax_tree, "a", dx=-0.055, dy=1.005)
    ys = list(range(len(BLOCK_ORDER)))[::-1]

    # Panel b
    axb = fig.add_subplot(gs[0, 1])
    for yy, b in zip(ys, BLOCK_ORDER):
        r = f3b[f3b.block_id == b].iloc[0]
        tmi_pct = 100*r.full_public_tmi_present/r.full_public_tmi_total
        mp_pct  = 100*r.full_public_mp_mip_present/r.full_public_mp_mip_total
        dumbbell(axb, yy, mp_pct, tmi_pct, C["mp"], C["tmi"], s=40, lw=1.4)
        hi, lo = (tmi_pct, mp_pct) if tmi_pct >= mp_pct else (mp_pct, tmi_pct)
        hi_c = C["tmi"] if tmi_pct >= mp_pct else C["mp"]
        lo_c = C["mp"] if tmi_pct >= mp_pct else C["tmi"]
        axb.text(hi+3, yy, f"{hi:.0f}", va="center", ha="left", fontsize=6, color=hi_c, fontweight="bold")
        axb.text(lo-3, yy, f"{lo:.0f}", va="center", ha="right", fontsize=6, color=lo_c, fontweight="bold")
        axb.text(0, yy+0.34, SHORT[b], va="bottom", ha="left", fontsize=6, color=C["ink"])
    axb.set_yticks([])
    axb.set_xlim(-12, 112)
    axb.set_xticks([0, 25, 50, 75, 100])
    axb.set_xlabel("Public genomes with candidate interval (%)")
    set_frame(axb, "b")
    axb.set_ylim(-0.6, len(BLOCK_ORDER)-0.05)
    axb.set_title("Full-atlas prevalence", loc="left", pad=16)
    axb.text(0.0, 1.05, "TMI (n=64)", transform=axb.transAxes, color=C["tmi"], fontsize=6, fontweight="bold", ha="left")
    axb.text(0.44, 1.05, "MP-MIP (n=43)", transform=axb.transAxes, color=C["mp"], fontsize=6, fontweight="bold", ha="left")
    panel_letter(axb, "b", dx=-0.02, dy=1.13)

    # Panel c
    axc = fig.add_subplot(gs[1, 1])
    for yy, b in zip(ys, BLOCK_ORDER):
        r = f3c[f3c.block_id == b].iloc[0]
        col = block_color[b]
        axc.plot([r.mantel_haenszel_odds_ratio_95ci_lower, r.mantel_haenszel_odds_ratio_95ci_upper],
                 [yy, yy], color=col, lw=1.4, solid_capstyle="round", zorder=2)
        axc.scatter([r.mantel_haenszel_common_odds_ratio_expected_direction], [yy], s=34, color=col,
                    edgecolor="white", linewidth=0.4, zorder=3)
        axc.text(1.18, yy+0.30, SHORT[b], va="bottom", ha="left", fontsize=6, color=C["ink"])
    axc.axvline(1, color="#999", linestyle="--", linewidth=0.8, zorder=1)
    axc.set_xscale("log")
    axc.set_xlim(0.8, 260)
    axc.set_xticks([1, 10, 100])
    axc.set_xticklabels(["1", "10", "100"])
    axc.set_yticks([])
    set_frame(axc, "b")
    axc.set_ylim(-0.6, len(BLOCK_ORDER)-0.05)
    axc.set_xlabel("Mantel–Haenszel odds ratio (95% CI)")
    axc.set_title("Association within mixed BioProjects", loc="left", pad=6)
    panel_letter(axc, "c", dx=-0.02, dy=1.05)

    # Panel d
    axd = fig.add_subplot(gs[2, 1])
    for yy, b in zip(ys, BLOCK_ORDER):
        sub = f3d[f3d.block_id == b].sort_values("start")
        strands = [1 if s == "+" else -1 for s in sub.strand]
        grey_mask = [not bool(v) for v in sub.stable_family]
        n_stable = int(sub.stable_family.sum())
        gene_arrows(axd, yy, strands, block_color[b], x0=0.0, glyph_w=1.0, gap=0.30, h=0.52,
                    grey_mask=grey_mask, grey=C["intervening"])
        axd.text(0, yy+0.42, f"{SHORT[b]} ({n_stable})", va="bottom", ha="left", fontsize=6, color=block_color[b])
    axd.set_xlim(-0.5, max(len(f3d[f3d.block_id == b]) for b in BLOCK_ORDER)*1.34)
    axd.set_ylim(-0.8, len(BLOCK_ORDER)-0.0)
    axd.set_xticks([])
    axd.set_yticks([])
    axd.axis("off")
    axd.set_title("Representative local gene order", loc="left", pad=6)
    axd.text(0.0, -0.06, "coloured, stable family;  grey, intervening ORF;  row scales independent",
             transform=axd.transAxes, fontsize=5.8, color=C["muted"], ha="left", va="top")
    panel_letter(axd, "d", dx=-0.02, dy=1.05)
    return fig


fig = build_fig3_v3()
for ax in fig.axes:
    if not ax.axison:
        ax.set_xticks([])
        ax.set_yticks([])

pdf3, png3 = save(fig, "Figure_4")
print("saved:", pdf3, png3)
