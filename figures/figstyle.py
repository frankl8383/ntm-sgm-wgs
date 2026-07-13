"""
figstyle.py — Canonical publication style + drawing primitives for the MAC/MGen manuscript.

Self-contained (no external skill dependency) so it is portable into any kernel via the
shared figure bundle. Load the `figure-style` skill separately for the §9 render-then-verify
QA checklist; use THIS module for the house look and the shared panel primitives, so every
figure in the paper shares one palette, one font ladder, and identical tree/heatmap/forest glyphs.

Usage:
    import figstyle as fs
    fs.setup()                       # sets Arial + rcParams; call once before plotting
    C = fs.C                         # canonical color dict
    fig, axes = plt.subplots(...)
    fs.panel_letter(ax, 'a')         # bold panel tag, outside axes
    fs.set_frame(ax, 'lb')           # keep only left+bottom spines
"""
from __future__ import annotations
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D

MM = 1.0 / 25.4          # millimetres -> inches
COLW = 183.0             # MGen double-column width in mm

# ----------------------------------------------------------------------------- palette
# Canonical, reconciled across all figures. Color IS the cross-reference (figure-style §4.1):
# an entity keeps ONE hue in every figure, line/fill/marker/text/heatmap alike.
C = {
    # neutrals / ink
    "ink":      "#242424",
    "muted":    "#6E6E6E",
    "grid":     "#D7D7D7",
    "branch":   "#B7B7B7",
    "neutral":  "#E8E8E8",
    "absent":   "#ECECEC",
    "leader":   "#9A9A9A",

    # ---- core scientific entities: the two accessory lineages ----
    "tmi":       "#3E6FA3",  # typical M. intracellulare (TMI)  -- blue
    "tmi_light": "#B9CCE2",
    "mp":        "#B85C55",  # M. paraintracellulare / MIP (MP-MIP) -- terracotta
    "mp_light":  "#E2BBB7",

    # ---- species / public-context groups ----
    "avium":       "#D19A32",  # M. avium public context -- amber
    "colombiense": "#8D6A9F",  # M. colombiense -- purple
    "timonense":   "#C58A3E",  # M. timonense boundary -- distinct warm amber
    "other_pub":   "#A9A9A9",  # other public label -- grey
    "chimaera":    "#8A78A8",  # chimaera-adjacent (S1 context) -- lavender
    "yongonense":  "#6FA1B0",  # yongonense-boundary (S1 context) -- teal-blue
    "marseillense":"#9C9C9C",  # marseillense (distant; omitted from S1 display) -- grey

    # ---- cohort / recovery status ----
    "direct":    "#9AA0A6",  # directly retained -- neutral grey
    "recovered": "#2F7F7A",  # rescued / recovered -- teal
    "excluded":  "#B95A55",  # residual within-MAC mixture / failure control -- muted red

    # ---- two-route reconstruction outcomes ----
    "primary":   "#3D8875",  # primary rule positive (20-80%) -- green
    "auxiliary": "#D9A441",  # auxiliary-only positive -- gold
    "negative":  "#E8E8E8",  # neither rule positive -- pale grey

    # ---- challenge-pair identities (Fig 4) ----
    "pairA": "#3F6E9A", "pairB": "#5B8C72", "pairC": "#C28A32", "pairD": "#9A6887",

    # ---- reserved alarm hue (never a data-series color) ----
    "alarm": "#B33A3A",
    "warning": "#B68A2D",
}

# lineage light-fill companions keyed by lineage color
LIGHT = {C["tmi"]: C["tmi_light"], C["mp"]: C["mp_light"]}

# diverging map for threshold-ratio matrices, centered at the semantic 1.0
DIVERGING = mpl.colors.LinearSegmentedColormap.from_list(
    "burden_ratio", ["#496E88", "#F7F7F7", "#B96346"]
)

def setup(base=8, mid=7, small=6):
    """Set Arial + a 3-size role ladder + outward ticks + vector-friendly output."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42,          # editable text in vector output
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

# ----------------------------------------------------------------------------- primitives
def set_frame(ax, keep="lb"):
    """Keep only the named spines. keep chars: l,r,t,b ; '' removes all."""
    m = {"l": "left", "r": "right", "t": "top", "b": "bottom"}
    for k, name in m.items():
        ax.spines[name].set_visible(k in keep)
    ax.tick_params(top=False, right=False)
    if "l" not in keep: ax.tick_params(left=False)
    if "b" not in keep: ax.tick_params(bottom=False)

def panel_letter(ax, letter, dx=-0.02, dy=1.02, size=11):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right", color=C["ink"])

def title_left(ax, text, pad=6, size=None):
    ax.set_title(text, loc="left", pad=pad, fontsize=size)

# ---- phylogeny -------------------------------------------------------------
def load_prune_ladderize(tree_path, keep_ids=None, anchor_id=None):
    """Read newick, optionally root on anchor, prune to keep_ids, ladderize(reverse)."""
    from Bio import Phylo
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
    """x = root-to-tip depth (branch lengths), y = evenly spaced terminals + parent means."""
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
    """Rectangular phylogram. Returns (x, y, ordered_terminals)."""
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
    """Aligned presence matrix beside a tree. presence_by_tip[tip.name][block] -> bool.
    block_colors: dict block-> color used when present. Absent cells get `absent` fill."""
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

# ---- dumbbell / lollipop ---------------------------------------------------
def dumbbell(ax, y, left_x, right_x, left_c, right_c, stem=None, s=34, lw=1.0):
    stem = stem or C["leader"]
    ax.plot([left_x, right_x], [y, y], color=stem, lw=lw, zorder=1, solid_capstyle="round")
    ax.scatter([left_x], [y], s=s, color=left_c, zorder=3, edgecolor="white", linewidth=0.4)
    ax.scatter([right_x], [y], s=s, color=right_c, zorder=3, edgecolor="white", linewidth=0.4)

def lollipop(ax, y, x, color, base=0.0, s=34, lw=0.9, orient="h"):
    if orient == "h":
        ax.plot([base, x], [y, y], color=C["leader"], lw=lw, zorder=1)
        ax.scatter([x], [y], s=s, color=color, zorder=3, edgecolor="white", linewidth=0.4)
    else:
        ax.plot([y, y], [base, x], color=C["leader"], lw=lw, zorder=1)
        ax.scatter([y], [x], s=s, color=color, zorder=3, edgecolor="white", linewidth=0.4)

# ---- forest plot (odds ratios on log axis) ---------------------------------
def forest(ax, ys, centers, los, his, colors, s=32, lw=1.2, cap=0.0):
    for y, c, lo, hi, col in zip(ys, centers, los, his, colors):
        ax.plot([lo, hi], [y, y], color=col, lw=lw, zorder=2, solid_capstyle="round")
        ax.scatter([c], [y], s=s, color=col, zorder=3, edgecolor="white", linewidth=0.4)

# ---- flow / schematic boxes ------------------------------------------------
def flow_box(ax, x, y, w, h, text, n=None, fc=None, ec=None, fs_txt=7.2, fs_n=8.5):
    fc = fc or C["neutral"]; ec = ec or "none"
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                         facecolor=fc, edgecolor=ec, linewidth=0.7, zorder=2)
    ax.add_patch(box)
    label = text if n is None else f"{text}"
    ax.text(x + w/2, y + h*0.60, text, ha="center", va="center", fontsize=fs_txt,
            color=C["ink"], zorder=3, linespacing=1.15)
    if n is not None:
        ax.text(x + w/2, y + h*0.26, f"n = {n}", ha="center", va="center",
                fontsize=fs_n, fontweight="bold", color=C["ink"], zorder=3)
    return (x, y, w, h)

def flow_arrow(ax, xy_from, xy_to, color=None, lw=0.9, rad=0.0):
    color = color or C["muted"]
    ar = FancyArrowPatch(xy_from, xy_to, arrowstyle="-|>", mutation_scale=8,
                         color=color, lw=lw, shrinkA=2, shrinkB=2,
                         connectionstyle=f"arc3,rad={rad}", zorder=1)
    ax.add_patch(ar)

# ---- gene-arrow track ------------------------------------------------------
def gene_arrows(ax, y, strands, color, x0=0.0, glyph_w=1.0, gap=0.28, h=0.5,
                grey_mask=None, grey=None):
    """Draw a row of gene arrows. strands: iterable of +1/-1. grey_mask: bools -> intervening ORF."""
    grey = grey or C["intervening"] if "intervening" in C else "#D9D9D9"
    x = x0
    for i, s in enumerate(strands):
        col = grey if (grey_mask and grey_mask[i]) else color
        _one_gene_arrow(ax, x, y, glyph_w, h, s, col)
        x += glyph_w + gap
    return x

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

# ---- categorical / value matrix -------------------------------------------
def value_matrix(ax, M, cmap, norm=None, row_labels=None, col_labels=None,
                 annot=None, annot_fmt="{:.1f}", text_thresh=None, grid=True):
    """Heatmap with per-cell annotation. M is 2D array; annot overrides text if given."""
    im = ax.imshow(M, aspect="auto", cmap=cmap, norm=norm, interpolation="none")
    nr, nc = M.shape
    if grid:
        for e in range(nc+1): ax.axvline(e-0.5, color="white", lw=0.8)
        for e in range(nr+1): ax.axhline(e-0.5, color="white", lw=0.8)
    if annot is not None or text_thresh is not None:
        for i in range(nr):
            for j in range(nc):
                v = M[i, j]
                txt = annot[i][j] if annot is not None else annot_fmt.format(v)
                if txt == "" or txt is None: continue
                dark = (text_thresh is not None and abs(v) >= text_thresh)
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.2,
                        color="white" if dark else C["ink"], zorder=3)
    ax.set_xticks(range(nc)); ax.set_yticks(range(nr))
    if col_labels is not None: ax.set_xticklabels(col_labels)
    if row_labels is not None: ax.set_yticklabels(row_labels)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    return im

def legend_swatches(ax, entries, loc="upper right", ncol=1, title=None, **kw):
    """entries: list of (color, label) or (marker, color, label)."""
    handles = []
    for e in entries:
        if len(e) == 2:
            c, lab = e
            handles.append(Rectangle((0,0),1,1, facecolor=c, edgecolor="none", label=lab))
        else:
            mk, c, lab = e
            handles.append(Line2D([],[], marker=mk, linestyle="", markerfacecolor=c,
                                  markeredgecolor="white", markersize=6, label=lab))
    leg = ax.legend(handles=handles, loc=loc, ncol=ncol, frameon=False,
                    handletextpad=0.5, columnspacing=1.0, borderaxespad=0.3,
                    title=title, **kw)
    return leg

def qa_report(fig):
    """Strong render-then-verify check. Catches THREE collision classes that a
    naive text-vs-text scan misses:
      (1) any visible text extending past the figure bounds (out-of-bounds),
      (2) a label straying over a DIFFERENT axes' data region (cross-panel),
      (3) text-vs-text overlaps within the same axes.
    Returns a list of (kind, detail) tuples; empty list == clean.
    Run this before save() on every deliverable figure."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    figbb = fig.get_window_extent(r)
    axbb = {ax: ax.get_window_extent(r) for ax in fig.axes}
    txts = [t for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    probs = []
    for t in txts:
        bb = t.get_window_extent(r)
        if bb.x1 > figbb.x1+1 or bb.x0 < figbb.x0-1 or bb.y1 > figbb.y1+1 or bb.y0 < figbb.y0-1:
            probs.append(("OOB", t.get_text()[:24]))
    for t in txts:
        home = t.axes; bb = t.get_window_extent(r)
        if home is None: continue
        for ax, abb in axbb.items():
            if ax is home: continue
            if bb.overlaps(abb):
                ix = max(0, min(bb.x1,abb.x1)-max(bb.x0,abb.x0))
                iy = max(0, min(bb.y1,abb.y1)-max(bb.y0,abb.y0))
                if ix*iy > 15:
                    probs.append(("CROSS_PANEL", t.get_text()[:24])); break
    tb = [(t, t.get_window_extent(r)) for t in txts]
    seen = set()
    for i,(a,ba) in enumerate(tb):
        for b,bb in tb[i+1:]:
            if ba.overlaps(bb) and a.axes is b.axes:
                k = tuple(sorted([a.get_text()[:15], b.get_text()[:15]]))
                if k not in seen: seen.add(k); probs.append(("TXT_TXT", k))
    # (4) text crossing a reference Line2D drawn in the same axes (axvline/axhline/
    #     threshold/leader). Straight ref lines are thin; a label sitting on one reads
    #     as a collision even though it is neither text nor an axes box.
    for t in txts:
        home = t.axes
        if home is None: continue
        bb = t.get_window_extent(r)
        for ln in home.get_lines():
            if not ln.get_visible(): continue
            try:
                lbb = ln.get_window_extent(r)
            except Exception:
                continue
            # only flag near-vertical/near-horizontal reference lines (thin bbox on one axis)
            thin = (lbb.width <= 3 or lbb.height <= 3)
            if not thin: continue
            # a zero-thickness ref line (axvline/axhline) has width or height == 0, for which
            # Bbox.overlaps() returns False; pad the thin dimension so it registers.
            lx0, lx1 = lbb.x0, lbb.x1; ly0, ly1 = lbb.y0, lbb.y1
            if lx1 - lx0 < 1: lx0 -= 1.5; lx1 += 1.5
            if ly1 - ly0 < 1: ly0 -= 1.5; ly1 += 1.5
            ix = max(0, min(bb.x1, lx1) - max(bb.x0, lx0))
            iy = max(0, min(bb.y1, ly1) - max(bb.y0, ly0))
            if ix*iy > 4:
                probs.append(("TXT_ON_LINE", t.get_text()[:24])); break
    return probs

def save(fig, stem, dpi=300):
    """Save PDF (vector) + PNG (raster) with the house defaults."""
    fig.savefig(f"{stem}.pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, facecolor="white", bbox_inches="tight")
    return f"{stem}.pdf", f"{stem}.png"

C["intervening"] = "#D9D9D9"
