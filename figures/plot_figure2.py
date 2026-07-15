import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import sys, os
import shutil

# Set up figstyle
sys.path.insert(0, os.getcwd())

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
import figstyle as fs

PAIR_ORDER = ["pairA", "pairB", "pairC", "pairD"]
PAIR_SAMPLES = {"pairA": ("Mi1", "Ma20"), "pairB": ("Mi2", "Mi22"),
                "pairC": ("Mi23", "Mi8"), "pairD": ("Mi25", "Mi4")}
PRIMARY_THRESHOLD, AUXILIARY_THRESHOLD = 161.8, 952.42

def _truth(v):
    return v is True or str(v).strip().lower() == "true"

def _classify(r):
    if _truth(r["either_route_primary_positive"]):
        return 2
    if _truth(r["either_route_auxiliary_positive"]):
        return 1
    return 0

def load_fig4_data():
    ani = pd.read_csv("source_data/Figure2_source_pair_ani.tsv", sep="\t")
    s1 = pd.read_csv("source_data/Figure2_stage1_detection_summary.tsv", sep="\t").sort_values("minor_percent")
    c10 = pd.read_csv("source_data/Figure2_complete_10_percent_conditions.tsv", sep="\t")
    pr = {p: i for i, p in enumerate(PAIR_ORDER)}
    c10["pr"] = c10["pair_id"].map(pr)
    c10["dr"] = c10.apply(lambda r: 0 if r["major"] == PAIR_SAMPLES[r["pair_id"]][0] else 1, axis=1)
    c10 = c10.sort_values(["pr", "dr"]).reset_index(drop=True)
    ratio_cols = [("strict_burden_20_80_per_mbp", PRIMARY_THRESHOLD),
                  ("meta_burden_20_80_per_mbp", PRIMARY_THRESHOLD),
                  ("strict_burden_10_90_per_mbp", AUXILIARY_THRESHOLD),
                  ("meta_burden_10_90_per_mbp", AUXILIARY_THRESHOLD)]
    ratios = np.column_stack([c10[col] / thr for col, thr in ratio_cols])
    c_rowlabels = [f"{r.pair_id[-1]}  {r.major}\u2192{r.minor}" for r in c10.itertuples()]
    c40 = pd.read_csv("source_data/Figure2_complete_40_conditions.tsv", sep="\t")
    seld = c40[c40["pair_id"].isin(["pairA", "pairD"])].copy()
    seld["pr"] = seld["pair_id"].map({"pairA": 0, "pairD": 1})
    seld["dr"] = seld.apply(lambda r: 0 if r["major"] == PAIR_SAMPLES[r["pair_id"]][0] else 1, axis=1)
    rk = (seld[["pair_id", "major", "minor", "pr", "dr", "minor_percent"]].drop_duplicates()
          .sort_values(["pr", "dr", "minor_percent"]).reset_index(drop=True))
    depths = [1_000_000, 2_000_000, 4_000_000]
    dmat = np.full((len(rk), 3), np.nan)
    for i, k in rk.iterrows():
        ss = seld[(seld.pair_id == k.pair_id) & (seld.major == k.major) & (seld.minor == k.minor)
                  & (seld.minor_percent == k.minor_percent)].set_index("total_pairs")
        for j, dp in enumerate(depths):
            dmat[i, j] = _classify(ss.loc[dp])
    d_rowlabels = [f"{k.pair_id[-1]}  {k.major} + {int(k.minor_percent)}% {k.minor}" for k in rk.itertuples()]
    return dict(ani=ani, s1=s1, ratios=ratios, c_rowlabels=c_rowlabels, dmat=dmat, d_rowlabels=d_rowlabels)

def _bounded_grid(ax, ncols, nrows, lw=1.0, color="white"):
    for j in range(1, ncols):
        ax.plot([j - 0.5, j - 0.5], [-0.5, nrows - 0.5], color=color, lw=lw, zorder=4, solid_capstyle="butt")
    for i in range(1, nrows):
        ax.plot([-0.5, ncols - 0.5], [i - 0.5, i - 0.5], color=color, lw=lw, zorder=4, solid_capstyle="butt")

def build_fig4(D):
    C, MM, COLW = fs.C, fs.MM, fs.COLW
    ani, s1 = D["ani"], D["s1"]
    ratios, c_rowlabels = D["ratios"], D["c_rowlabels"]
    dmat, d_rowlabels = D["dmat"], D["d_rowlabels"]
    fs.setup()
    fig = plt.figure(figsize=(COLW * MM, 145 * MM), facecolor="white")
    ax_a = fig.add_axes([0.105, 0.615, 0.185, 0.300])
    ax_b = fig.add_axes([0.440, 0.615, 0.205, 0.300])
    ax_c = fig.add_axes([0.790, 0.605, 0.170, 0.290])
    ax_d = fig.add_axes([0.150, 0.115, 0.700, 0.360])

    # (a) challenge-pair similarity
    for i, row in ani.iterrows():
        col = C[PAIR_ORDER[i]]
        ax_a.plot([row.ani_min, row.ani_max], [i, i], color=col, lw=2.4, solid_capstyle="round", zorder=2)
        ax_a.scatter(row.ani_mean, i, s=28, color=col, edgecolor="white", linewidth=0.5, zorder=3)
        anchor = max(row.ani_max, row.ani_mean)
        ax_a.text(anchor + 0.030, i, f"{row.ani_mean:.2f}", va="center", ha="left", fontsize=6.0, color=C["ink"])
    ax_a.set_yticks(range(4))
    ax_a.set_yticklabels([f"{PAIR_ORDER[i][-1]}  {ani.samples[i]}" for i in range(4)], fontsize=6.2)
    ax_a.invert_yaxis()
    ax_a.set_xlim(98.83, 99.63); ax_a.set_xticks([98.9, 99.1, 99.3, 99.5]); ax_a.set_ylim(3.65, -0.65)
    fs.set_frame(ax_a, "b"); ax_a.tick_params(axis="y", length=0)
    ax_a.set_xlabel("Reciprocal FastANI (%)", fontsize=7)
    ax_a.set_title("Challenge-pair similarity", loc="left", pad=5, fontsize=8)
    ax_a.text(0.0, -0.235, "A, B cross TMI \u00d7 MP-MIP lineages;\nC, D closely related within MP-MIP",
              transform=ax_a.transAxes, ha="left", va="top", fontsize=5.8, color=C["muted"], linespacing=1.25)
    fs.panel_letter(ax_a, "a", dx=-0.34, dy=1.02)

    # (b) dilution detection
    x = s1["minor_percent"].to_numpy()
    ax_b.plot(x, s1["primary_positive"], color=C["primary"], marker="o", ms=4.0, lw=1.5, zorder=3,
              label="Primary rule (20\u201380%)")
    ax_b.plot(x, s1["combined_positive"], color=C["auxiliary"], marker="s", ms=3.8, lw=1.5, zorder=3,
              label="Combined rule")
    ax_b.set_xlim(-2, 52); ax_b.set_ylim(-1.2, 25.5)
    ax_b.set_xticks([0, 5, 10, 20, 30, 50]); ax_b.set_yticks([0, 6, 12, 18, 24])
    fs.set_frame(ax_b, "lb")
    ax_b.set_xlabel("Minor-source reads (%)", fontsize=7)
    ax_b.set_ylabel("Positive dependent windows / 24", fontsize=6.8)
    ax_b.set_title("Clean-reference dilution stage", loc="left", pad=5, fontsize=8)
    ax_b.annotate("combined\n\u2191 at 10%", xy=(10, 24), xytext=(1.5, 18.5), fontsize=5.8,
                  color=C["auxiliary"], ha="left", va="center",
                  arrowprops=dict(arrowstyle="-", color=C["auxiliary"], lw=0.6))
    ax_b.annotate("primary\n\u2191 at 20%", xy=(20, 24), xytext=(22.5, 12.0), fontsize=5.8,
                  color=C["primary"], ha="left", va="center",
                  arrowprops=dict(arrowstyle="-", color=C["primary"], lw=0.6))
    ax_b.legend(loc="lower right", fontsize=5.9, handlelength=1.6, borderaxespad=0.5)
    ax_b.text(0.0, -0.235, "144 windows are dependent (mean overlap 58.4%, max 97.0%):\n"
              "technical windows, not independent replicates",
              transform=ax_b.transAxes, ha="left", va="top", fontsize=5.6, color=C["muted"], linespacing=1.25)
    fs.panel_letter(ax_b, "b", dx=-0.22, dy=1.02)

    # (c) 10% ratio matrix
    logv = np.log2(np.clip(ratios, 0.125, 8.0))
    norm = mpl.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    ax_c.imshow(logv, aspect="auto", cmap=fs.DIVERGING, norm=norm, interpolation="none")
    for i in range(8):
        for j in range(4):
            dark = abs(logv[i, j]) >= 1.5
            ax_c.text(j, i, f"{ratios[i, j]:.1f}", ha="center", va="center", fontsize=6.0,
                      color="white" if dark else C["ink"], zorder=5)
    _bounded_grid(ax_c, 4, 8, lw=0.9)
    ax_c.plot([1.5, 1.5], [-0.5, 7.5], color="white", lw=2.0, zorder=4)
    for hb in (1.5, 3.5, 5.5):
        ax_c.plot([-0.5, 3.5], [hb, hb], color="white", lw=1.6, zorder=4)
    ax_c.set_xticks(range(4)); ax_c.set_xticklabels(["Strict", "Meta", "Strict", "Meta"], fontsize=6.0)
    ax_c.set_yticks(range(8)); ax_c.set_yticklabels(c_rowlabels, fontsize=5.8)
    ax_c.tick_params(length=0, pad=2)
    for s in ax_c.spines.values():
        s.set_visible(False)
    ax_c.set_ylim(9.4, -1.7); ax_c.set_xlim(-0.55, 3.55)
    ax_c.text(0.5, -1.02, "Primary rule", ha="center", va="center", fontsize=6.4, color=C["ink"])
    ax_c.text(2.5, -1.02, "Auxiliary rule", ha="center", va="center", fontsize=6.4, color=C["ink"])
    ax_c.text(1.5, 8.9, "cell = burden \u00f7 threshold;  > 1 = positive", ha="center", va="center",
              fontsize=5.5, color=C["muted"])
    ax_c.set_title("10% reconstruction", loc="left", pad=5, fontsize=8)
    fs.panel_letter(ax_c, "c", dx=-0.46, dy=1.045)

    # (d) nested-depth grid
    cmap_d = mpl.colors.ListedColormap([C["negative"], C["auxiliary"], C["primary"]])
    ax_d.imshow(dmat, aspect="auto", interpolation="none", cmap=cmap_d, vmin=-0.5, vmax=2.5)
    ax_d.set_xticks(range(3)); ax_d.set_xticklabels(["1M pairs", "2M pairs", "4M pairs"], fontsize=6.5)
    ax_d.set_yticks(range(12)); ax_d.set_yticklabels(d_rowlabels, fontsize=6.0)
    ax_d.tick_params(length=0)
    sym = {0: "\u2013", 1: "A", 2: "P"}
    for i in range(12):
        for j in range(3):
            st = int(dmat[i, j])
            ax_d.text(j, i, sym[st], ha="center", va="center", fontsize=6.6, fontweight="bold",
                      color="white" if st == 2 else C["ink"])
    _bounded_grid(ax_d, 3, 12, lw=1.0)
    for hb in (2.5, 5.5, 8.5):
        ax_d.plot([-0.5, 2.5], [hb, hb], color="white", lw=2.2, zorder=4)
    ax_d.plot([-0.5, 2.5], [5.5, 5.5], color=C["ink"], lw=0.7, zorder=5)
    for s in ax_d.spines.values():
        s.set_visible(False)
    ax_d.set_ylim(11.5, -0.5); ax_d.set_xlim(-0.5, 2.5)
    ax_d.set_title("Nested-depth sensitivity of the complete two-route rule", loc="left", pad=5, fontsize=8)
    leg = [Patch(facecolor=C["primary"], label="P  primary positive"),
           Patch(facecolor=C["auxiliary"], label="A  auxiliary-only positive"),
           Patch(facecolor=C["negative"], label="\u2013  neither rule positive")]
    ax_d.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.075),
                ncol=3, fontsize=6.2, handlelength=1.1, columnspacing=1.8, handletextpad=0.5)
    ax_d.text(0.5, -0.175, "Depth tiers are nested prefixes of one seeded window (not independent);\n"
              "they define the tested operating range, not a universal limit.",
              transform=ax_d.transAxes, ha="center", va="top", fontsize=5.6, color=C["muted"], linespacing=1.25)
    fs.panel_letter(ax_d, "d", dx=-0.135, dy=1.02)
    return fig

D4 = load_fig4_data()
fig4 = build_fig4(D4)
fs.save(fig4, "Figure_2")
