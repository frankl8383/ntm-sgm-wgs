import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import sys, os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle, Patch
import importlib, types

# Load figstyle
import importlib.util
spec = importlib.util.spec_from_file_location("figstyle", "figstyle.py")
fs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)
fs.setup()

C, MM, COLW = fs.C, fs.MM, fs.COLW

BLOCK_ORDER = ["TMI_aromatic_catabolism_associated_block", "TMI_methyltransferase_hydrolase_cupin_block",
               "MP_MIP_nitrogen_redox_associated_block", "MP_MIP_oxidoreductase_associated_block"]
SHORT = {"TMI_aromatic_catabolism_associated_block": "Aromatic",
         "TMI_methyltransferase_hydrolase_cupin_block": "Methyl./cupin",
         "MP_MIP_nitrogen_redox_associated_block": "N-redox",
         "MP_MIP_oxidoreductase_associated_block": "Oxidoreductase"}
HDR = {"TMI_aromatic_catabolism_associated_block": "Aromatic",
       "TMI_methyltransferase_hydrolase_cupin_block": "Methyl./\ncupin",
       "MP_MIP_nitrogen_redox_associated_block": "N-redox",
       "MP_MIP_oxidoreductase_associated_block": "Oxido-\nreductase"}
BLOCK_LINE = {b: ("tmi" if b.startswith("TMI_") else "mp") for b in BLOCK_ORDER}
DATASETS = ["Separate public set", "Complete genomes", "One per project", "Near-clone reduced"]
LINES_B = ["TMI", "MP-MIP"]
CTX = ["block_and_both_flanks_supported", "both_flanks_supported_block_absent", "partial_or_divergent_context"]


def _tf(v):
    return v is True or str(v).strip().lower() == "true"


def load_fig5_data():
    ann = pd.read_csv("source_data/Figure5a_external_annotations.tsv", sep="\t")
    qc = pd.read_csv("source_data/Figure5a_external_qc.tsv", sep="\t")
    calls = pd.read_csv("source_data/Figure5a_external_interval_calls.tsv", sep="\t")
    reads = pd.read_csv("source_data/Figure5a_external_family_read_support.tsv", sep="\t")
    b5 = pd.read_csv("source_data/Figure5b_public_sensitivities.tsv", sep="\t")
    c5 = pd.read_csv("source_data/Figure5c_flanking_context.tsv", sep="\t")

    A = ann.merge(qc[["sample_id", "checkm2_contamination", "contigs_ge_500_bp"]], on="sample_id")
    A["lin_rank"] = A["frozen_expected_lineage"].map({"TMI": 0, "MP_MIP": 1})
    A["qual"] = A["principal_lineage_qualified"].apply(_tf)
    A["qual_rank"] = (~A["qual"]).astype(int)
    A = A.sort_values(["lin_rank", "qual_rank", "sample_id"]).reset_index(drop=True)

    pres = np.zeros((len(A), 4), dtype=bool)
    ci = calls.set_index(["sample_id", "block_id"])
    for i, sid in enumerate(A.sample_id):
        for j, b in enumerate(BLOCK_ORDER):
            try:
                pres[i, j] = _tf(ci.loc[(sid, b), "syntenic_block_present"])
            except KeyError:
                pres[i, j] = False
    reads = reads.assign(agree=reads["read_assembly_content_agreement"].apply(_tf))
    dis = [(r.sample_id, r.block_id) for r in reads[~reads.agree].itertuples()]
    agree_str = f"{int(reads.agree.sum())}/{len(reads)}"

    B = {ds: {b: {} for b in BLOCK_ORDER} for ds in DATASETS}
    for ds in DATASETS:
        for b in BLOCK_ORDER:
            for lin in LINES_B:
                r = b5[(b5.dataset == ds) & (b5.block_id == b) & (b5.lineage == lin)].iloc[0]
                B[ds][b][lin] = (int(r.present), int(r.total), float(r.prevalence))

    def ctx_count(block, lin, cls):
        return int(c5[(c5.block_id == block) & (c5.accessory_lineage == lin)
                      & (c5.structural_context_class == cls)].genomes.sum())
    return dict(A=A, pres=pres, dis=dis, agree=agree_str, B=B, c5=c5, ctx_count=ctx_count)


def _panel_a(ax, A, pres, dis, agree, C):
    n = len(A)
    n_tmi = int((A.frozen_expected_lineage == "TMI").sum())
    n_qual = int(A.qual.sum())
    ax.set_xlim(-3.5, 3.6); ax.set_ylim(n + 0.05, -1.55)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("External evaluation of 9 PRJNA983112 genomes", loc="left", pad=6, fontsize=8)

    # anchor-qualified group band (rows 0..n_qual-1)
    ax.add_patch(Rectangle((-3.42, -0.44), 0.13, n_qual - 1 + 0.88, facecolor=C["recovered"],
                           edgecolor="none", alpha=0.85, zorder=1))
    ax.text(-3.28, (n_qual - 1) / 2.0, "anchor-qualified (\u2265 98.60% ANI)", rotation=90,
            ha="center", va="center", fontsize=5.6, color=C["recovered"], fontweight="bold")

    # column headers
    ax.text(-3.15, -0.72, "Isolate", ha="left", va="center", fontsize=5.6, color=C["muted"])
    ax.text(-1.15, -0.72, "ANI", ha="right", va="center", fontsize=5.6, color=C["muted"])
    ax.text(-0.62, -0.72, "QC", ha="center", va="center", fontsize=5.6, color=C["muted"])
    for j, b in enumerate(BLOCK_ORDER):
        ax.text(j, -0.80, HDR[b], ha="center", va="center", fontsize=5.7, fontweight="bold",
                color=C[BLOCK_LINE[b]], linespacing=1.0)

    dis_set = set(dis)
    for i, row in A.iterrows():
        sid = row.sample_id
        lin = "tmi" if row.frozen_expected_lineage == "TMI" else "mp"
        # lineage chip
        ax.add_patch(Rectangle((-3.18, i - 0.40), 0.10, 0.80, facecolor=C[lin], edgecolor="none"))
        ax.text(-3.02, i, sid, ha="left", va="center", fontsize=6.2,
                fontweight="bold" if row.qual else "normal",
                color=C["ink"] if row.qual else C["muted"])
        ax.text(-1.15, i, f"{row.principal_lineage_anchor_max_ani:.1f}", ha="right", va="center",
                fontsize=5.8, color=C["ink"] if row.qual else C["muted"])
        # QC triangle for external_qc_pass == False
        if not _tf(row.external_qc_pass):
            ax.scatter(-0.62, i, marker="^", s=26, facecolor=C["warning"], edgecolor="white",
                       linewidth=0.4, zorder=4)
        # presence cells
        for j, b in enumerate(BLOCK_ORDER):
            present = pres[i, j]
            fc = C[BLOCK_LINE[b]] if present else C["absent"]
            ax.add_patch(Rectangle((j - 0.46, i - 0.44), 0.92, 0.88, facecolor=fc,
                                    edgecolor="white", linewidth=0.6, zorder=2))
            if (sid, b) in dis_set:
                ax.add_patch(Rectangle((j - 0.46, i - 0.44), 0.92, 0.88, facecolor="none",
                                       edgecolor=C["alarm"], linewidth=1.4, zorder=5))
    # lineage divider between TMI and MP-MIP (starts right of the rotated anchor label)
    ax.plot([-3.05, 3.5], [n_tmi - 0.5, n_tmi - 0.5], color=C["muted"], lw=0.6, zorder=3)
    ax.text(3.55, (n_tmi - 1) / 2.0, "TMI", rotation=90, ha="center", va="center",
            fontsize=5.6, color=C["tmi"], fontweight="bold")
    ax.text(3.55, (n_tmi + n - 1) / 2.0, "MP-MIP", rotation=90, ha="center", va="center",
            fontsize=5.6, color=C["mp"], fontweight="bold")
    # caption
    ax.text(-3.42, n + 0.55,
            f"Filled = one-contig interval call.  {agree}/family-content agreement.\n"
            "\u25b2 QC contig-count warning (mv17).   \u25a1 red = read/assembly disagreement (mv8).",
            ha="left", va="top", fontsize=5.3, color=C["muted"], linespacing=1.3)


def _panel_b(ax, B, C):
    ax.set_xlim(-4.7, 7.7); ax.set_ylim(len(DATASETS) + 0.55, -1.65)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Interval prevalence across 4 public dataset definitions", loc="left", pad=6,
                 fontsize=8, x=-0.145)
    # block headers spanning T,M sub-cols
    for g, b in enumerate(BLOCK_ORDER):
        xc = 2 * g + 0.5
        ax.text(xc, -1.28, HDR[b], ha="center", va="center", fontsize=5.7, fontweight="bold",
                color=C[BLOCK_LINE[b]], linespacing=1.0)
        ax.text(2 * g, -0.60, "T", ha="center", va="center", fontsize=5.6, color=C["tmi"], fontweight="bold")
        ax.text(2 * g + 1, -0.60, "M", ha="center", va="center", fontsize=5.6, color=C["mp"], fontweight="bold")
    for r, ds in enumerate(DATASETS):
        ax.text(-0.7, r, ds, ha="right", va="center", fontsize=5.8, color=C["ink"])
        for g, b in enumerate(BLOCK_ORDER):
            for k, lin in enumerate(LINES_B):
                present, total, prev = B[ds][b][lin]
                base = C["tmi"] if lin == "TMI" else C["mp"]
                fc = mpl.colors.to_rgba(base, 0.12 + 0.88 * prev)
                x = 2 * g + k
                ax.add_patch(Rectangle((x - 0.47, r - 0.45), 0.94, 0.90, facecolor=fc,
                                       edgecolor="white", linewidth=0.6, zorder=2))
                ax.text(x, r, f"{present}/{total}", ha="center", va="center", fontsize=4.9,
                        color="white" if prev >= 0.6 else C["ink"], zorder=4)
    # group separators
    for xb in (1.5, 3.5, 5.5):
        ax.plot([xb, xb], [-0.5, len(DATASETS) - 0.5], color="white", lw=1.8, zorder=3)
    ax.text(-4.6, len(DATASETS) + 0.05,
            "T = TMI, M = MP-MIP.  Cell = genomes with interval / tested;\n"
            "fill depth tracks prevalence. Lineage-specific pattern holds\nacross all four dataset definitions.",
            ha="left", va="top", fontsize=5.3, color=C["muted"], linespacing=1.3)


def _panel_c(ax, ctx_count, C):
    oxi = "MP_MIP_oxidoreductase_associated_block"
    met = "TMI_methyltransferase_hydrolase_cupin_block"
    # (label, block, accessory_lineage, present_colour_key, ytick_colour_key)
    rows = [
        ("MP-MIP \u00b7 own",      oxi, "MI_MP_MIP_lineage", "mp",  "mp",  3.5),
        ("TMI \u00b7 opposing",    oxi, "MI_TMI_lineage",    "mp",  "tmi", 2.7),
        ("TMI \u00b7 own",         met, "MI_TMI_lineage",    "tmi", "tmi", 1.3),
        ("MP-MIP \u00b7 opposing", met, "MI_MP_MIP_lineage", "tmi", "mp",  0.5),
    ]
    flanks_absent_fc, partial_fc = C["intervening"], C["absent"]
    ax.set_xlim(-14.0, 30.5); ax.set_ylim(-0.75, 4.65)
    ax.set_yticks([])
    ax.set_xticks([0, 10, 20]); ax.set_xlabel("Complete genomes", fontsize=7)
    fs.set_frame(ax, "b")
    ax.spines["bottom"].set_bounds(0, 26)
    ax.set_title("Flanking context at two candidate loci (complete assemblies)", loc="left",
                 pad=6, fontsize=8)
    hh = 0.30
    for lab, blk, lin, pcol, tcol, y in rows:
        present = ctx_count(blk, lin, CTX[0])
        fabs = ctx_count(blk, lin, CTX[1])
        part = ctx_count(blk, lin, CTX[2])
        ax.text(-6.3, y, lab, ha="left", va="center", fontsize=6.0, color=C[tcol])
        x0 = 0.0
        for val, fc, hatch, tc in [(present, C[pcol], None, "white"),
                                   (fabs, flanks_absent_fc, "////", C["ink"]),
                                   (part, partial_fc, None, C["muted"])]:
            if val > 0:
                ax.add_patch(Rectangle((x0, y - hh), val, 2 * hh, facecolor=fc, edgecolor="white",
                                       linewidth=0.6, hatch=hatch, zorder=2))
                ax.text(x0 + val / 2.0, y, str(val), ha="center", va="center", fontsize=5.6,
                        color=tc, zorder=4)
                x0 += val
    # group titles (far left) + brackets between titles and row labels
    ax.text(-13.6, 3.1, "MP-MIP\noxidoreductase", ha="left", va="center", fontsize=6.0,
            fontweight="bold", color=C["mp"], linespacing=1.05)
    ax.text(-13.6, 0.9, "TMI\nmethyl./cupin", ha="left", va="center", fontsize=6.0,
            fontweight="bold", color=C["tmi"], linespacing=1.05)
    ax.plot([-7.0, -7.0], [2.7 - hh, 3.5 + hh], color=C["mp"], lw=1.1)
    ax.plot([-7.0, -7.0], [0.5 - hh, 1.3 + hh], color=C["tmi"], lw=1.1)
    # legend
    handles = [Patch(facecolor=C["mp"], edgecolor="white", label="block + both flanks (lineage colour)"),
               Patch(facecolor=flanks_absent_fc, hatch="////", edgecolor=C["muted"],
                     label="both flanks present, block absent"),
               Patch(facecolor=partial_fc, edgecolor="white", label="partial / divergent context")]
    ax.legend(handles=handles, loc="lower right", fontsize=5.6, handlelength=1.4,
              borderaxespad=0.4, labelspacing=0.35)


def build_fig5(D):
    C, MM, COLW = fs.C, fs.MM, fs.COLW
    A, pres = D["A"], D["pres"]
    fs.setup()
    fig = plt.figure(figsize=(COLW * MM, 140 * MM), facecolor="white")
    ax_a = fig.add_axes([0.070, 0.400, 0.400, 0.500])
    ax_b = fig.add_axes([0.570, 0.400, 0.400, 0.500])
    ax_c = fig.add_axes([0.090, 0.090, 0.870, 0.200])

    _panel_a(ax_a, A, pres, D["dis"], D["agree"], C)
    _panel_b(ax_b, D["B"], C)
    _panel_c(ax_c, D["ctx_count"], C)

    fs.panel_letter(ax_a, "a", dx=-0.13, dy=1.02)
    fs.panel_letter(ax_b, "b", dx=-0.235, dy=1.02)
    fs.panel_letter(ax_c, "c", dx=-0.055, dy=1.05)
    return fig


D5 = load_fig5_data()
fig5 = build_fig5(D5)
fs.save(fig5, "Figure_5")