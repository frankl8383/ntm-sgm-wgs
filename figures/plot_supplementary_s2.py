import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import shutil, tarfile, os

import pandas as pd
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 40)

import sys; import figstyle as fs; fs.setup()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

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

df2a = pd.read_csv("phylogeny/panel_avium_timonense_64.tsv", sep="\t", dtype=str).fillna("")
df2b = pd.read_csv("phylogeny/panel_colombiense_29.tsv", sep="\t", dtype=str).fillna("")
anc2a = df2a[df2a.source=="anchor"].tree_id.iat[0]
anc2b = df2b[df2b.source=="anchor"].tree_id.iat[0]
OMIT_A = "PUB_GCF_010723675.1"
keep2a = [t for t in df2a.tree_id if t != OMIT_A]

tree2a = fs.load_prune_ladderize("phylogeny/avium_timonense_64.treefile",
                                 keep_ids=keep2a, anchor_id=anc2a)
tree2b = fs.load_prune_ladderize("phylogeny/colombiense_29.treefile",
                                 anchor_id=anc2b)

fig2 = plt.figure(figsize=(fs.COLW*fs.MM, 150*fs.MM))
axA = fig2.add_axes([0.045, 0.145, 0.44, 0.75])
axB = fig2.add_axes([0.555, 0.145, 0.42, 0.75])

xA, yA, ordA = fs.draw_tree(axA, tree2a, lw=0.5)
stA = decorate_tree(axA, tree2a, xA, yA, ordA, df2a, label_fs=6.0, timonense_directlabel=True)
place_local_labels(axA, stA["local_pos"], stA["maxx"], dx=0.03, fs_=6.2, min_gap=1.6)
axA.set_xlabel("Substitutions per complete-core SNP site", fontsize=6.8)
axA.set_xlim(-stA["maxx"]*0.015, stA["maxx"]*1.22)
axA.set_ylim(-1.5, len(ordA)+1.5)
axA.set_xticks(np.arange(0,0.13,0.02)); axA.tick_params(axis="x", length=2.4)
fs.panel_letter(axA, "a", dx=-0.01, dy=1.0)
axA.set_title("M. avium / M. timonense public context", loc="left", pad=10, fontsize=7.6)
axA.text(0, 1.012, "63 of 64 genomes shown; 8,921 variable sites; one distant public M. timonense tip omitted",
         transform=axA.transAxes, fontsize=5.8, color=fs.C["muted"], va="bottom", ha="left")

xB, yB, ordB = fs.draw_tree(axB, tree2b, lw=0.5)
stB = decorate_tree(axB, tree2b, xB, yB, ordB, df2b, label_fs=6.0)
place_local_labels(axB, stB["local_pos"], stB["maxx"], dx=0.03, fs_=6.2, min_gap=1.6)
axB.set_xlabel("Substitutions per complete-core SNP site", fontsize=6.8)
axB.set_xlim(-stB["maxx"]*0.015, stB["maxx"]*1.20)
axB.set_ylim(-1.5, len(ordB)+1.5)
axB.set_xticks(np.arange(0,0.36,0.05)); axB.tick_params(axis="x", length=2.4)
fs.panel_letter(axB, "b", dx=-0.01, dy=1.0)
axB.set_title("M. colombiense public context", loc="left", pad=10, fontsize=7.6)
axB.text(0, 1.012, "29 genomes; 851 complete-core variable sites",
         transform=axB.transAxes, fontsize=5.8, color=fs.C["muted"], va="bottom", ha="left")

handles2 = [
    Line2D([],[],marker="o",ls="",mfc=fs.C["avium"],mec="none",ms=5,label="M. avium label"),
    Line2D([],[],marker="o",ls="",mfc=fs.C["colombiense"],mec="none",ms=5,label="M. timonense (a) / M. colombiense (b) label"),
    Line2D([],[],marker="o",ls="",mfc=fs.C["other_pub"],mec="none",ms=5,label="other public label"),
    Line2D([],[],marker="*",ls="",mfc=fs.C["muted"],mec="white",mew=0.4,ms=8,label="local genome ($^{\\mathregular{R}}$ recovered)"),
    Line2D([],[],marker="D",ls="",mfc=fs.C["ink"],mec="none",ms=5,label="rooting anchor"),
]
fig2.legend(handles=handles2, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.005),
            frameon=False, handletextpad=0.5, columnspacing=1.6, fontsize=6.2)

axA.set_title(""); axB.set_title("")
axA.set_title("", loc="left"); axB.set_title("", loc="left")
axA.set_title("", loc="center"); axB.set_title("", loc="center")
for ax in (axA, axB):
    for t in list(ax.texts):
        if t.get_text().startswith(("63 of 64","29 genomes")):
            t.remove()

fig2.text(0.055, 0.955, "M. avium / M. timonense public context",
          fontsize=7.6, fontweight="bold", ha="left", va="bottom", style="italic")
fig2.text(0.055, 0.935, "63 of 64 genomes shown; 8,921 variable sites; one distant public M. timonense tip omitted",
          fontsize=5.7, color=fs.C["muted"], ha="left", va="top")
fig2.text(0.565, 0.955, "M. colombiense public context",
          fontsize=7.6, fontweight="bold", ha="left", va="bottom", style="italic")
fig2.text(0.565, 0.935, "29 genomes; 851 complete-core variable sites",
          fontsize=5.7, color=fs.C["muted"], ha="left", va="top")

for ax,let in [(axA,"a"),(axB,"b")]:
    for t in list(ax.texts):
        if t.get_text() in ("a","b") and t.get_fontweight()=="bold" and t.get_fontsize()>=10:
            t.remove()
fig2.text(0.012, 0.955, "a", fontsize=11, fontweight="bold", ha="left", va="bottom")
fig2.text(0.522, 0.955, "b", fontsize=11, fontweight="bold", ha="left", va="bottom")

p2 = fs.save(fig2, "Supplementary_Figure_S2")