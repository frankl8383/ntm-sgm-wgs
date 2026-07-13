import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import shutil, tarfile, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.path import Path
from matplotlib.lines import Line2D

pd.set_option("display.width", 200); pd.set_option("display.max_columns", 50)


import sys
import figstyle as fs; fs.setup()

# Load data
fa2_df = pd.read_csv("source_data/Figure2a_atlas_selection_flow.tsv", sep="\t")
comp = pd.read_csv("source_data/Figure2a_atlas_reporting_label_composition.tsv", sep="\t")
b2_df = pd.read_csv("source_data/Figure2b_local_group_composition.tsv", sep="\t")
c2_df = pd.read_csv("source_data/Figure2c_local_public_type_ani.tsv", sep="\t")
d2_df = pd.read_csv("source_data/Figure2d_local_ani_margins.tsv", sep="\t")

C = fs.C

def _lum(hexc):
    hexc = hexc.lstrip("#"); r,g,b = (int(hexc[i:i+2],16)/255 for i in (0,2,4))
    f = lambda c: c/12.92 if c<=0.03928 else ((c+0.055)/1.055)**2.4
    return 0.2126*f(r)+0.7152*f(g)+0.0722*f(b)
def tc_for(bg):
    return "#FFFFFF" if _lum(bg) < 0.42 else fs.C["ink"]

def box(ax, x, y, w, h, text, n, fc, tc=None, fs_txt=7.0, fs_n=8.5):
    tc = tc or tc_for(fc)
    ax.add_patch(FancyBboxPatch((x,y), w, h, boxstyle="round,pad=0.004,rounding_size=0.010",
                 facecolor=fc, edgecolor="none", linewidth=0.7, zorder=2))
    ax.text(x+w/2, y+h*0.60, text, ha="center", va="center", fontsize=fs_txt,
            color=tc, zorder=3, linespacing=1.15)
    ax.text(x+w/2, y+h*0.235, f"n = {n}", ha="center", va="center",
            fontsize=fs_n, fontweight="bold", color=tc, zorder=3)

def elbow(ax, pts, color=None, lw=0.9, head=True, ms=8):
    color = color or fs.C["muted"]
    verts = list(pts); codes = [Path.MOVETO]+[Path.LINETO]*(len(verts)-1)
    p = Path(verts, codes)
    style = "-|>" if head else "-"
    ar = FancyArrowPatch(path=p, arrowstyle=style, mutation_scale=ms, color=color,
                         lw=lw, shrinkA=0, shrinkB=0, joinstyle="miter", capstyle="butt", zorder=1)
    ax.add_patch(ar); return ar

def crop_panels(fig, regions, pad=6, prefix="crop"):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    W,H = fig.canvas.get_width_height()
    buf = np.asarray(fig.canvas.buffer_rgba())
    from PIL import Image
    paths={}
    for name, ax in regions.items():
        bb = ax.get_window_extent(r)
        x0=max(0,int(bb.x0)-pad); x1=min(W,int(bb.x1)+pad)
        y0=max(0,int(bb.y0)-pad); y1=min(H,int(bb.y1)+pad)
        row0=H-y1; row1=H-y0
        Image.fromarray(buf[row0:row1, x0:x1]).save(f"{prefix}_{name}.png")
        paths[name]=f"{prefix}_{name}.png"
    return paths

# ---------- Figure 2 data prep ----------
fa2 = fa2_df.set_index("stage")["n"].to_dict()
comp["short"] = comp.reporting_species.str.replace("Mycobacterium ","M. ", regex=False)

named = ["Mycobacterium avium","Mycobacterium intracellulare","Mycobacterium colombiense",
         "Mycobacterium paraintracellulare","Mycobacterium mantenii","Mycobacterium marseillense"]
comp_named = comp[comp.reporting_species.isin(named)].set_index("reporting_species").reindex(named)
tail_n = int(comp[~comp.reporting_species.isin(named)].n.sum())
seg_labels = ["M. avium","M. intracellulare","M. colombiense","M. paraintracellulare",
              "M. mantenii","M. marseillense","Other MAC labels"]
seg_n = list(comp_named.n.astype(int)) + [tail_n]
seg_col = [C["avium"], C["tmi"], C["colombiense"], C["mp"], C["muted"], C["branch"], C["grid"]]

GC = {"M_intracellulare_complex":C["tmi"], "M_avium_timonense_boundary":C["avium"], "M_colombiense":C["colombiense"]}
GL = {"M_intracellulare_complex":"M. intracellulare complex",
      "M_avium_timonense_boundary":"M. avium public context", "M_colombiense":"M. colombiense"}

b2 = b2_df.set_index("broad_analysis_panel").reindex(
     ["M_intracellulare_complex","M_avium_timonense_boundary","M_colombiense"])

c2 = c2_df[
     ["sample_id","broad_analysis_panel","best_current_public_ani","best_direct_type_ani"]
     ].sort_values("best_current_public_ani").reset_index(drop=True)
d2 = d2_df[
     ["sample_id","broad_analysis_panel","public_distinct_species_ani_margin"]
     ].sort_values("public_distinct_species_ani_margin").reset_index(drop=True)

plt.close("all")
fig = plt.figure(figsize=(fs.COLW*fs.MM, 140*fs.MM))
outer = fig.add_gridspec(2,1, height_ratios=[0.46,1.0], left=0.055, right=0.985,
                         top=0.95, bottom=0.115, hspace=0.34)
ax_a = fig.add_subplot(outer[0])
bot = outer[1].subgridspec(1,3, width_ratios=[0.80,1.06,1.14], wspace=0.42)
ax_b = fig.add_subplot(bot[0,0]); ax_c = fig.add_subplot(bot[0,1]); ax_d = fig.add_subplot(bot[0,2])

# ================= (a) atlas construction + composition =================
ax_a.set_xlim(0,1); ax_a.set_ylim(0,1); ax_a.axis("off")
fs.panel_letter(ax_a,"a",dx=0.004,dy=1.04)
ax_a.set_title("Current public MAC atlas: construction and composition", loc="left", x=0.036, pad=2)
fy, fh = 0.31, 0.56
box(ax_a,0.010,fy,0.150,fh,"Current\nrecords",1283,C["neutral"],fs_txt=6.7,fs_n=8.2)
box(ax_a,0.210,fy,0.150,fh,"Assembly-QC\npass",957,C["neutral"],fs_txt=6.7,fs_n=8.2)
box(ax_a,0.410,fy,0.150,fh,"QC-first pair +\nBioSample dedup.",481,C["tmi_light"],fs_txt=6.7,fs_n=8.2)
yc_a = fy+fh/2
elbow(ax_a,[(0.160,yc_a),(0.210,yc_a)])
elbow(ax_a,[(0.360,yc_a),(0.410,yc_a)])
bx0, bw, bary, barh = 0.635, 0.355, 0.60, 0.20
tot = sum(seg_n); cur = bx0
for lab,n,col in zip(seg_labels,seg_n,seg_col):
    seg = bw*n/tot
    ax_a.add_patch(Rectangle((cur,bary),seg,barh,facecolor=col,edgecolor="white",linewidth=0.6,zorder=2))
    if seg > 0.045:
        ax_a.text(cur+seg/2,bary+barh/2,str(n),ha="center",va="center",fontsize=6.3,
                  color="white" if col in (C["avium"],C["tmi"]) else C["ink"],zorder=3)
    cur += seg
ax_a.text(bx0,bary+barh+0.06,"Retained atlas by reporting label (n = 481)",fontsize=6.7,ha="left",va="bottom")
lx0, ly0, dx_col, dy_row = bx0, bary-0.16, 0.185, 0.135
for i,(lab,n,col) in enumerate(zip(seg_labels,seg_n,seg_col)):
    r,cix = i//2, i%2
    xx = lx0 + cix*dx_col; yy = ly0 - r*dy_row
    ax_a.add_patch(Rectangle((xx,yy-0.035),0.028,0.07,facecolor=col,edgecolor="none",zorder=3))
    ax_a.text(xx+0.036,yy,f"{lab} ({n})",va="center",ha="left",fontsize=5.7,color=C["ink"])

# ================= (b) local group composition =================
fs.set_frame(ax_b,"b"); ax_b.spines["left"].set_visible(False)
fs.panel_letter(ax_b,"b",dx=-0.02,dy=1.05)
ax_b.set_title("Updated local cohort", loc="left", x=0.0, pad=6)
y = np.arange(len(b2))[::-1]
ax_b.barh(y, b2.direct, color=C["direct"], height=0.5, label="direct", zorder=2)
ax_b.barh(y, b2.recovered, left=b2.direct, color=C["recovered"], height=0.5, label="recovered", zorder=2)
for i,g in enumerate(b2.index):
    t = int(b2.direct[g]+b2.recovered[g])
    ax_b.text(t+0.4, y[i], str(t), va="center", fontsize=7.5, fontweight="bold", color=C["ink"])
    ax_b.text(0.0, y[i]+0.46, GL[g], va="bottom", ha="left", fontsize=6.4, color=C["ink"])
ax_b.set_yticks([]); ax_b.set_ylim(-0.6, len(b2)-0.1)
ax_b.set_xlim(0,19); ax_b.set_xticks([0,5,10,15]); ax_b.set_xlabel("Local genomes (n)")
ax_b.legend(loc="lower right", handlelength=1.1, handletextpad=0.5, borderaxespad=0.4)

# ================= (c) ANI to nearest public vs type anchor =================
fs.set_frame(ax_c,"b"); ax_c.spines["left"].set_visible(False)
fs.panel_letter(ax_c,"c",dx=-0.05,dy=1.05)
ax_c.set_title("Nearest public genome and type anchor", loc="left", x=0.0, pad=6)
yc = np.arange(len(c2))[::-1]
for i,row in c2.iterrows():
    col = GC[row.broad_analysis_panel]
    ax_c.plot([row.best_direct_type_ani,row.best_current_public_ani],[yc[i],yc[i]],
              color=C["leader"], lw=0.8, zorder=1)
    ax_c.scatter(row.best_direct_type_ani, yc[i], s=20, marker="D", facecolor="white",
                 edgecolor=col, linewidth=0.9, zorder=3)
    ax_c.scatter(row.best_current_public_ani, yc[i], s=21, marker="o", facecolor=col,
                 edgecolor="white", linewidth=0.4, zorder=3)
ax_c.set_yticks(yc); ax_c.set_yticklabels(c2.sample_id)
for tick,g in zip(ax_c.get_yticklabels(), c2.broad_analysis_panel):
    tick.set_color(GC[g]); tick.set_fontweight("bold")
ax_c.set_xlim(97.5,100.05); ax_c.set_xticks([97.5,98,98.5,99,99.5,100])
ax_c.set_ylim(-0.8, len(c2)-0.2)
ax_c.set_xlabel("ANI to reference (%)"); ax_c.tick_params(axis="y",length=0)
leg_c = [Line2D([],[],marker="o",ls="",mfc="#666",mec="white",ms=5,label="nearest current\npublic genome"),
         Line2D([],[],marker="D",ls="",mfc="white",mec="#666",ms=5,label="best direct\ntype anchor")]
ax_c.legend(handles=leg_c, loc="center left", bbox_to_anchor=(0.0,0.52),
            handletextpad=0.3, labelspacing=0.6, borderaxespad=0.3, fontsize=5.6)

# ================= (d) separation margins =================
fs.set_frame(ax_d,"b"); ax_d.spines["left"].set_visible(False)
fs.panel_letter(ax_d,"d",dx=-0.05,dy=1.05)
ax_d.set_title("Separation from nearest distinct label", loc="left", x=0.0, pad=6)
yd = np.arange(len(d2))[::-1]
NARROW = 0.15
for i,row in d2.iterrows():
    col = GC[row.broad_analysis_panel]
    m = row.public_distinct_species_ani_margin
    ax_d.plot([0.003, m],[yd[i],yd[i]], color=C["leader"], lw=0.8, zorder=1)
    ax_d.scatter(m, yd[i], s=22, color=col, edgecolor="white", linewidth=0.4, zorder=3)
ax_d.axvline(NARROW, color=C["muted"], linestyle=(0,(4,2)), lw=0.8, zorder=2)
ax_d.set_xscale("log"); ax_d.set_xlim(0.003,9)
ax_d.xaxis.set_major_locator(mticker.FixedLocator([0.01,0.1,1]))
ax_d.xaxis.set_major_formatter(mticker.FixedFormatter(["0.01","0.1","1"]))
ax_d.set_yticks(yd); ax_d.set_yticklabels(d2.sample_id)
for tick,g in zip(ax_d.get_yticklabels(), d2.broad_analysis_panel):
    tick.set_color(GC[g]); tick.set_fontweight("bold")
ax_d.set_ylim(-0.8, len(d2)+1.6)
ax_d.set_xlabel("ANI margin to nearest distinct\npublic label (percentage points, log)")
ax_d.tick_params(axis="y",length=0)
ax_d.text(0.12, len(d2)+0.55, "narrow-margin\nprompt region", fontsize=5.6, color=C["muted"],
          ha="right", va="center", linespacing=1.05)
leg_g = [Line2D([],[],marker="o",ls="",mfc=GC[g],mec="white",ms=5,label=GL[g]) for g in GC]
ax_d.legend(handles=leg_g, loc="upper right", handletextpad=0.3, labelspacing=0.3, borderaxespad=0.4, fontsize=5.4)

probs = fs.qa_report(fig)
print("QA Figure 2:", probs)
if probs:
    raise RuntimeError(f"Figure 2 layout QA failed: {probs}")
# Per-panel crops for perceptual QA review (opt-in): set QA_CROPS=1 to emit f2_*.png
if os.environ.get("QA_CROPS"):
    paths = crop_panels(fig, {"a":ax_a,"b":ax_b,"c":ax_c,"d":ax_d}, prefix="f2")

pdf2, png2 = fs.save(fig, "Figure_2")
print("saved:", pdf2, png2)
