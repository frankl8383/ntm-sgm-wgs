import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import shutil, tarfile, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


pd.set_option("display.width", 200); pd.set_option("display.max_columns", 40)

import sys; import figstyle as fs; fs.setup()

GROUP_COLOR = {
    "TMI": fs.C["tmi"], "MP-MIP": fs.C["mp"],
    "chimaera-adjacent": fs.C.get("chimaera","#8A78A8"),
    "yongonense-boundary": fs.C.get("yongonense", fs.C["avium"]),
    "M. avium": fs.C["avium"],
    "M. timonense": fs.C["colombiense"],
    "M. colombiense": fs.C["colombiense"],
    "other": fs.C["other_pub"],
}

def decorate_tree(ax, tree, x, y, order, panel_df, *, label_fs=5.6,
                  timonense_directlabel=False):
    idx = panel_df.set_index("tree_id")
    maxx = max(x.values())
    local_pos = []
    n_local = n_rec = 0
    for tip in order:
        r = idx.loc[tip.name]
        grp = r["analysis_group"]; src = r["source"]
        if src == "anchor":
            ax.scatter([x[tip]],[y[tip]], marker="D", s=17, facecolor=fs.C["ink"],
                       edgecolor="none", zorder=4, clip_on=False)
            ax.annotate("rooting anchor", (x[tip], y[tip]), xytext=(x[tip]+maxx*0.014, y[tip]),
                        fontsize=label_fs-0.3, color=fs.C["muted"], va="center", ha="left",
                        style="italic", clip_on=False)
            continue
        col = GROUP_COLOR.get(grp, fs.C["other_pub"])
        if src == "local":
            ax.scatter([x[tip]],[y[tip]], marker="*", s=42, facecolor=col,
                       edgecolor="white", linewidth=0.4, zorder=5, clip_on=False)
            lab = r["sample_label"]
            if str(r["recovered"]).strip()=="True":
                lab = lab + r"$^{\mathregular{R}}$"; n_rec += 1
            local_pos.append((x[tip], y[tip], lab, col)); n_local += 1
        else:
            ax.scatter([x[tip]],[y[tip]], marker="o", s=6, facecolor=col,
                       edgecolor="none", zorder=3, clip_on=False)
            if timonense_directlabel and grp=="M. timonense":
                ax.annotate("M. timonense", (x[tip], y[tip]),
                            xytext=(x[tip]+maxx*0.012, y[tip]), fontsize=label_fs-0.3,
                            color=col, va="center", ha="left", style="italic", clip_on=False)
    return dict(maxx=maxx, local_pos=local_pos, n_local=n_local, n_rec=n_rec)

def place_local_labels(ax, local_pos, maxx, *, dx=0.012, fs_=5.6, min_gap=1.35):
    lp = sorted(local_pos, key=lambda t: t[1])
    placed_y = []
    for lx, ly, lab, col in lp:
        ty = ly
        for py in placed_y:
            if abs(ty-py) < min_gap: ty = py + min_gap
        placed_y.append(ty)
        ax.annotate(lab, (lx, ly), xytext=(lx+maxx*dx, ty), fontsize=fs_,
                    fontweight="bold", color=col, va="center", ha="left", clip_on=False,
                    annotation_clip=False)

df1 = pd.read_csv("phylogeny/panel_intracellulare_complex_94.tsv", sep="\t", dtype=str).fillna("")
anchor1 = df1[df1.source=="anchor"].tree_id.iat[0]
keep1 = [t for t in df1.tree_id if df1.set_index("tree_id").loc[t,"analysis_group"]!="marseillense"]
tree1 = fs.load_prune_ladderize("phylogeny/intracellulare_complex_94.treefile",
                                keep_ids=keep1, anchor_id=anchor1)

fig1 = plt.figure(figsize=(fs.COLW*fs.MM, 220*fs.MM))
ax1 = fig1.add_axes([0.045, 0.055, 0.90, 0.86])
x1, y1, order1 = fs.draw_tree(ax1, tree1, lw=0.5)
st1 = decorate_tree(ax1, tree1, x1, y1, order1, df1, label_fs=5.8)
place_local_labels(ax1, st1["local_pos"], st1["maxx"], dx=0.013, fs_=5.8, min_gap=1.45)

ax1.set_xlabel("Substitutions per complete-core SNP site", fontsize=7)
ax1.set_xlim(-st1["maxx"]*0.012, st1["maxx"]*1.12)
ax1.set_ylim(-1.5, len(order1)+2.5)
ax1.tick_params(axis="x", length=2.6)

fig1.text(0.045, 0.965, "Current public context of the local M. intracellulare complex",
          fontsize=9, fontweight="bold", ha="left", va="bottom")
fig1.text(0.045, 0.952, "90 of 94 genomes displayed; four distant M. marseillense records omitted for legibility (complete Newick supplied separately)",
          fontsize=6.2, color=fs.C["muted"], ha="left", va="top")
fs.panel_letter(ax1, "a", dx=-0.006, dy=1.0)

handles = [
    Line2D([],[],marker="o",ls="",mfc=fs.C["tmi"],mec="none",ms=5,label="TMI (n=40)"),
    Line2D([],[],marker="o",ls="",mfc=fs.C["mp"],mec="none",ms=5,label="MP-MIP (n=39)"),
    Line2D([],[],marker="o",ls="",mfc=fs.C.get("chimaera","#8A78A8"),mec="none",ms=5,label="chimaera-adjacent (n=6)"),
    Line2D([],[],marker="o",ls="",mfc=fs.C.get("yongonense","#6FA1B0"),mec="none",ms=5,label="yongonense-boundary (n=4)"),
    Line2D([],[],marker="*",ls="",mfc=fs.C["muted"],mec="white",mew=0.4,ms=8,label="local genome (star + label)"),
    Line2D([],[],marker="D",ls="",mfc=fs.C["ink"],mec="none",ms=5,label="rooting anchor"),
]
leg = ax1.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.985),
                 frameon=False, handletextpad=0.5, labelspacing=0.55, fontsize=6.4,
                 borderaxespad=0.0)
ax1.text(0.02, 0.72, r"$^{\mathregular{R}}$ recovered local genome (rescued from mixed/failed assembly)",
         transform=ax1.transAxes, fontsize=6.2, color=fs.C["ink"], va="top", ha="left")

ax1.set_xticks(np.arange(0, 0.26, 0.05))
ax1.set_xlim(-st1["maxx"]*0.012, st1["maxx"]*1.12)

maxx1 = st1["maxx"]

for t in list(ax1.texts):
    t.remove()

anchor_tip = next(t for t in order1 if df1.set_index('tree_id').loc[t.name,'source']=='anchor')
ax1.annotate("rooting anchor", (x1[anchor_tip], y1[anchor_tip]),
             xytext=(x1[anchor_tip]+maxx1*0.014, y1[anchor_tip]),
             fontsize=5.5, color=fs.C["muted"], va="center", ha="left", style="italic", clip_on=False)
place_local_labels(ax1, st1["local_pos"], maxx1, dx=0.028, fs_=5.8, min_gap=1.45)
ax1.text(0.02, 0.72, r"$^{\mathregular{R}}$ recovered local genome (rescued from mixed/failed assembly)",
         transform=ax1.transAxes, fontsize=6.2, color=fs.C["ink"], va="top", ha="left")
ax1.set_xlim(-maxx1*0.012, maxx1*1.16)
ax1.set_xticks(np.arange(0,0.26,0.05))

p1 = fs.save(fig1, "Supplementary_Figure_S1")