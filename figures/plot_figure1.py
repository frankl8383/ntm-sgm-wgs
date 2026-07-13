import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.path import Path
from matplotlib.lines import Line2D

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)

# Load figstyle
sys.path.insert(0, ".")
import shutil
import figstyle as fs
fs.setup()

# Load data
dfs = {}
dfs["Figure1a_cohort_flow.tsv"] = pd.read_csv("source_data/Figure1a_cohort_flow.tsv", sep="\t")
dfs["Figure1b_synthetic_mixture_benchmark.tsv"] = pd.read_csv("source_data/Figure1b_synthetic_mixture_benchmark.tsv", sep="\t")
dfs["Figure1c_residual_mixture_burden.tsv"] = pd.read_csv("source_data/Figure1c_residual_mixture_burden.tsv", sep="\t")

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

THR = 161.826223
b1 = dfs["Figure1b_synthetic_mixture_benchmark.tsv"].copy()
c1 = dfs["Figure1c_residual_mixture_burden.tsv"].copy()
c1 = c1.sort_values("meta_mixed_sites_per_mbp").reset_index(drop=True)

plt.close("all")

fig = plt.figure(figsize=(fs.COLW*fs.MM, 122*fs.MM))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.06], width_ratios=[0.92, 1.08],
                      left=0.075, right=0.985, top=0.945, bottom=0.135, hspace=0.42, wspace=0.30)
ax_a = fig.add_subplot(gs[0, :]); ax_b = fig.add_subplot(gs[1, 0]); ax_c = fig.add_subplot(gs[1, 1])

# ============================ (a) cohort flow ============================
ax_a.set_xlim(0,1); ax_a.set_ylim(0,1); ax_a.axis("off")
fs.panel_letter(ax_a, "a", dx=-0.005, dy=1.03)
ax_a.set_title("Benchmark-gated recovery of interpretable genomes", loc="left", x=0.028, pad=2)
H=0.30
def bx(x,yc,w,text,n,fc): box(ax_a,x,yc-H/2,w,H,text,n,fc,fs_txt=6.7,fs_n=8.2)
bx(0.005,0.50,0.155,"Presumed\nSGM–NTM",38,C["neutral"])
bx(0.230,0.83,0.170,"Directly\nretained",13,C["direct"])
bx(0.230,0.50,0.170,"Recovery\ncandidates",11,C["neutral"])
bx(0.230,0.17,0.170,"Not advanced\nto recovery",14,C["neutral"])
bx(0.500,0.67,0.180,"Recovered\n(both routes)",8,C["recovered"])
bx(0.500,0.25,0.180,"Residual mixture\n(excluded)",3,C["excluded"])
bx(0.795,0.67,0.200,"Final interpretable\nMAC/SGM genomes",21,C["neutral"])
elbow(ax_a,[(0.160,0.50),(0.195,0.50),(0.195,0.83),(0.230,0.83)])
elbow(ax_a,[(0.160,0.50),(0.230,0.50)])
elbow(ax_a,[(0.160,0.50),(0.195,0.50),(0.195,0.17),(0.230,0.17)])
elbow(ax_a,[(0.400,0.50),(0.455,0.50),(0.455,0.67),(0.500,0.67)])
elbow(ax_a,[(0.400,0.50),(0.455,0.50),(0.455,0.25),(0.500,0.25)])
elbow(ax_a,[(0.400,0.83),(0.755,0.83),(0.755,0.735),(0.795,0.735)])
elbow(ax_a,[(0.680,0.67),(0.795,0.67)])
ax_a.text(0.588,0.95,"13 retained + 8 recovered",ha="center",va="center",fontsize=6.2,color=C["muted"])
ax_a.text(0.590,0.045,"excluded before lineage & accessory analyses",ha="center",va="center",fontsize=6.0,color=C["muted"])

# ============================ (b) benchmark scatter ============================
fs.set_frame(ax_b,"lb"); fs.panel_letter(ax_b,"b",dx=-0.155,dy=1.03)
ax_b.set_title("Synthetic-mixture recruitment benchmark",loc="left",x=-0.11,pad=6)
cross=b1[b1.benchmark.str.startswith("cross")]; near=b1[b1.benchmark.eq("nearMAC50")]
ax_b.scatter(cross.sensitivity_percent,cross.precision_percent,s=30,marker="o",facecolor=C["primary"],edgecolor="white",linewidth=0.5,zorder=3)
ax_b.scatter(near.sensitivity_percent,near.precision_percent,s=42,marker="D",facecolor=C["excluded"],edgecolor="white",linewidth=0.5,zorder=3)
cx,cy=cross.sensitivity_percent.mean(),cross.precision_percent.mean()
ax_b.annotate("cross-genus controls\n25 / 65 / 95% target",(cx,cy),xytext=(-6,-16),textcoords="offset points",fontsize=6.1,ha="right",va="top",arrowprops={"arrowstyle":"-","color":C["leader"],"lw":0.5})
nr=near.iloc[0]
ax_b.annotate("50:50 within-MAC\nfailure control",(nr.sensitivity_percent,nr.precision_percent),xytext=(8,2),textcoords="offset points",fontsize=6.1,ha="left",va="center")
ax_b.set_xlim(84.5,92.5); ax_b.set_ylim(40,103)
ax_b.set_xticks([85,86,87,88,89,90,91,92]); ax_b.set_yticks([40,60,80,100])
ax_b.set_xlabel("Target-pair sensitivity (%)"); ax_b.set_ylabel("Recruited-pair precision (%)")

# ============================ (c) residual mixture dumbbells ============================
fs.set_frame(ax_c,"lb"); fs.panel_letter(ax_c,"c",dx=-0.13,dy=1.03)
ax_c.set_title("Residual within-MAC mixture vs fixed threshold",loc="left",x=-0.09,pad=6)
for i,row in c1.iterrows():
    col=C["recovered"] if row.recovered else C["excluded"]
    s_,m_=row.strict_mixed_sites_per_mbp,row.meta_mixed_sites_per_mbp
    ax_c.plot([s_,m_],[i,i],color=C["leader"],lw=0.9,zorder=1,solid_capstyle="round")
    ax_c.scatter(s_,i,s=26,marker="o",facecolor=col,edgecolor="white",linewidth=0.45,zorder=3)
    ax_c.scatter(m_,i,s=27,marker="s",facecolor=col,edgecolor="white",linewidth=0.45,zorder=3)
ax_c.axvline(THR,color=C["auxiliary"],linestyle=(0,(4,2)),lw=1.0,zorder=2)
ax_c.set_xscale("log"); ax_c.set_xlim(12,4600); ax_c.set_ylim(-0.7,len(c1)-0.3)
ax_c.xaxis.set_major_locator(mticker.FixedLocator([100,1000]))
ax_c.xaxis.set_major_formatter(mticker.FixedFormatter(["100","1000"]))
ax_c.set_yticks(np.arange(len(c1)))
ax_c.set_yticklabels([f"{s}†" if s=="Mi32" else s for s in c1.sample_id])
for tick,rec in zip(ax_c.get_yticklabels(),c1.recovered):
    tick.set_color(C["recovered"] if rec else C["excluded"])
ax_c.set_xlabel("Intermediate-frequency sites per callable Mb (log scale)")
ax_c.tick_params(axis="y",length=0)
ax_c.text(THR*1.15,len(c1)-1.0,"fixed 161.8 sites/Mb\nthreshold",fontsize=6.0,color="#8A651E",ha="left",va="top")
leg=[Line2D([],[],marker="o",ls="",mfc="#777",mec="white",ms=5,label="strict route"),
     Line2D([],[],marker="s",ls="",mfc="#777",mec="white",ms=5,label="metaSPAdes bin"),
     Line2D([],[],marker="o",ls="",mfc=C["recovered"],mec="white",ms=5,label="recovered"),
     Line2D([],[],marker="o",ls="",mfc=C["excluded"],mec="white",ms=5,label="residual mixture"),
     Line2D([],[],marker="",ls="",label="† type-boundary caution")]
ax_c.legend(handles=leg,loc="upper left",ncol=1,handletextpad=0.3,labelspacing=0.32,borderaxespad=0.3,bbox_to_anchor=(0.0,1.0))

ax_c.get_legend().remove()
ax_c.legend(handles=leg, loc="lower right", ncol=1, handletextpad=0.3, labelspacing=0.32,
            borderaxespad=0.6, bbox_to_anchor=(1.0, 0.0))

pdf1, png1 = fs.save(fig, "Figure_1")
print("saved:", pdf1, png1)