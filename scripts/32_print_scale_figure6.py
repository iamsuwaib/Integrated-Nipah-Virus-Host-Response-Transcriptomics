#!/usr/bin/env python3
"""
32_print_scale_figure6.py

Rebuilds Figure 6 at print size (6.5 in wide, 600 dpi) directly from the refined secretome/surfaceome tables
written by 08_secretome_surfaceome.R (advanced_analyses/tables):
  secretome_surfaceome_candidate_table_refined.csv
  secretome_surfaceome_category_summary_refined.csv
Figure 6A: candidates with at least one significantly upregulated context, ranked by the study-level priority score
(GSE32902 and GSE33133 counted as one HUVEC study); Figure 6B: candidate genes per physiological class.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT = os.environ.get("NIPAH_PROJECT_DIR", "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics")
TAB = os.path.join(PROJECT, "advanced_analyses", "tables")
OUT = os.path.join(PROJECT, "Revised_manuscript_R2", "02_figures")
DPI = 600
plt.rcParams.update({"font.family": "DejaVu Sans"})
comp_col = {"secreted_or_surface-associated": "#009E73", "secretome": "#D55E00", "surfaceome": "#0072B2"}
comp_lab = {"secreted_or_surface-associated": "Secreted or surface-associated", "secretome": "Secretome", "surfaceome": "Surfaceome"}

a = pd.read_csv(os.path.join(TAB, "secretome_surfaceome_candidate_table_refined.csv"))
a = a[a.evidence_level != "detected_not_significant"].sort_values("priority_score", ascending=True).reset_index(drop=True)
mk = {"high": "o", "moderate": "^", "supportive": "s", "limited": "P"}
n = len(a)
W, H = 6.5, 0.45 + n * 0.205 + 0.62
fig = plt.figure(figsize=(W, H), dpi=DPI)
ax = fig.add_axes([0.10, 0.62 / H, 0.60, n * 0.205 / H])
yp = np.arange(n)
for y, r in zip(yp, a.itertuples()):
    c = comp_col[r.compartment]
    ax.hlines(y, 0, r.priority_score, color=c, lw=1.0, alpha=0.7, zorder=1)
    ax.scatter(r.priority_score, y, marker=mk[r.evidence_level], s=10 + r.n_significant_up * 9, c=c,
               edgecolor="none" if r.evidence_level == "limited" else "#222222", lw=0.3, zorder=3)
ax.set_yticks(yp); ax.set_yticklabels(a.gene, fontsize=7, fontweight="bold")
ax.set_xlim(0, a.priority_score.max() * 1.1); ax.set_ylim(-0.7, n - 0.3)
ax.grid(axis="x", color="#E6E6E6", lw=0.5); ax.set_axisbelow(True)
ax.tick_params(labelsize=7, length=2, width=0.5, pad=1.5)
ax.set_xlabel("Priority score", fontsize=7.5, fontweight="bold")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.text(0.01, 0.985, "A", fontsize=15, fontweight="bold", va="top")
fig.text(0.07, 0.975, "Extracellular and cell-surface mediators in NiV host-response signatures", fontsize=7.8, fontweight="bold", va="top")
lx = 0.74
y0 = 0.90
def block(title, handles, y):
    fig.text(lx, y, title, fontsize=7, fontweight="bold", va="top")
    leg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(lx - 0.01, y - 0.02 * 5 / H * 4.4), frameon=False, fontsize=6.6,
                     handletextpad=0.3, labelspacing=0.45, borderpad=0)
    fig.add_artist(leg)
h1 = [Line2D([], [], marker=m, ls="", mfc="#777777", mec="#222222", mew=0.3, ms=5.5, label=k) for k, m in mk.items()]
h2 = [Line2D([], [], marker="o", ls="-", color=c, mfc=c, mec="none", ms=5, lw=1.0, label=comp_lab[k]) for k, c in comp_col.items()]
h3 = [Line2D([], [], marker="o", ls="", mfc="#777777", mec="#222222", mew=0.3, ms=np.sqrt(10 + k * 9) * 0.95, label=str(k)) for k in (1, 3, 5, 6)]
top_in = H - 0.5
def at(y_in): return y_in / H
fig.text(lx, at(top_in), "Evidence level", fontsize=7, fontweight="bold", va="top")
fig.add_artist(fig.legend(handles=h1, loc="upper left", bbox_to_anchor=(lx - 0.01, at(top_in - 0.2)), frameon=False, fontsize=6.6, handletextpad=0.3, labelspacing=0.4, borderpad=0))
fig.text(lx, at(top_in - 1.15), "Compartment", fontsize=7, fontweight="bold", va="top")
fig.add_artist(fig.legend(handles=h2, loc="upper left", bbox_to_anchor=(lx - 0.01, at(top_in - 1.35)), frameon=False, fontsize=6.6, handletextpad=0.3, labelspacing=0.4, borderpad=0))
fig.text(lx, at(top_in - 2.1), "Significant upregulated\ncontexts", fontsize=7, fontweight="bold", va="top", linespacing=1.1)
fig.add_artist(fig.legend(handles=h3, loc="upper left", bbox_to_anchor=(lx - 0.01, at(top_in - 2.5)), frameon=False, fontsize=6.6, handletextpad=0.3, labelspacing=0.5, borderpad=0))
fig.savefig(os.path.join(OUT, "Figure_6A_candidate_ranking_print.png"), dpi=DPI, facecolor="white")
plt.close(fig)

b = pd.read_csv(os.path.join(TAB, "secretome_surfaceome_category_summary_refined.csv")).sort_values("n_genes", ascending=True).reset_index(drop=True)
H2 = 0.45 + len(b) * 0.3 + 0.7
fig = plt.figure(figsize=(6.5, H2), dpi=DPI)
ax = fig.add_axes([0.40, 0.55 / H2, 0.56, len(b) * 0.3 / H2])
ax.barh(np.arange(len(b)), b.n_genes, color=b.compartment.map(comp_col), edgecolor="#222222", lw=0.4, height=0.7, zorder=3)
for i, r in enumerate(b.itertuples()):
    ax.text(r.n_genes + 0.15, i, str(int(r.n_genes)), fontsize=6.8, va="center")
ax.set_yticks(np.arange(len(b))); ax.set_yticklabels([c[0].upper()+c[1:] for c in b.physiology_class], fontsize=7)
ax.set_xlim(0, b.n_genes.max() * 1.15); ax.set_ylim(-0.6, len(b) - 0.4)
ax.grid(axis="x", color="#E6E6E6", lw=0.5, zorder=0); ax.set_axisbelow(True)
ax.tick_params(labelsize=7, length=2, width=0.5, pad=1.5)
ax.set_xlabel("Number of candidate genes", fontsize=7.5, fontweight="bold")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.text(0.01, 0.985, "B", fontsize=15, fontweight="bold", va="top")
fig.text(0.07, 0.975, "Extracellular and cell-surface mediator classes", fontsize=7.8, fontweight="bold", va="top")
hh = [Patch(fc=c, ec="#222222", lw=0.4, label=comp_lab[k]) for k, c in comp_col.items()] if False else None
from matplotlib.patches import Patch
hh = [Patch(fc=c, ec="#222222", lw=0.4, label=comp_lab[k]) for k, c in comp_col.items()]
ax.legend(handles=hh, loc="lower right", title="Compartment", frameon=False, fontsize=6.6, title_fontsize=7, handletextpad=0.4, labelspacing=0.4, borderpad=0.2)
fig.savefig(os.path.join(OUT, "Figure_6B_mediator_classes_print.png"), dpi=DPI, facecolor="white")
print(a[["gene", "priority_score", "evidence_level"]].iloc[::-1].to_string())
