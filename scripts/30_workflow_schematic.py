#!/usr/bin/env python3
"""30_workflow_schematic.py - Figure 1 (study design / workflow) at print size."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

PROJECT = os.environ.get("NIPAH_PROJECT_DIR", "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics")
OUT = os.path.join(PROJECT, "Revised_manuscript_R2", "02_figures")
os.makedirs(OUT, exist_ok=True)
W, H = 6.5, 3.9
fig = plt.figure(figsize=(W, H), dpi=600)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
NAVY = "#1F4E79"
def box(x, y, w, h, fc, ec, title, lines, tsize=8, lsize=7, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.07", fc=fc, ec=ec, lw=0.8,
                                ls="--" if dashed else "-"))
    ax.text(x + 0.1, y + h - 0.1, title, fontsize=tsize, fontweight="bold", color=NAVY, va="top")
    for i, l in enumerate(lines):
        ax.text(x + 0.1, y + h - 0.1 - 0.2 - i * 0.155, l, fontsize=lsize, va="top", color="#222222")
def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=1.1, color="#4472C4"))

# left: HUVEC experiment group
ax.add_patch(FancyBboxPatch((0.08, 1.72), 1.72, 2.1, boxstyle="round,pad=0.0,rounding_size=0.08", fc="none", ec="#7F9CC4", lw=0.8, ls="--"))
ax.text(0.16, 3.76, "One underlying HUVEC experiment", fontsize=6.6, color="#335A8A", va="top", fontweight="bold")
ax.text(0.16, 3.61, "(two related GEO series)", fontsize=6.6, color="#335A8A", va="top")
box(0.17, 2.72, 1.54, 0.78, "#DCEBFA", "#7F9CC4", "GSE32902", ["Human HUVEC", "NiV vs mock"])
box(0.17, 1.82, 1.54, 0.78, "#DCEBFA", "#7F9CC4", "GSE33133", ["Human HUVEC", "NiV, NiV-dC, mock"])
box(0.08, 0.35, 1.72, 1.1, "#E3F1D8", "#8DB56F", "GSE310471", ["Independent in vivo study", "AGM lung and tonsil", "baseline, 3, 4, 5 DPI"])
# middle column
box(2.35, 2.85, 1.6, 0.85, "#F4F6F8", "#8A93A0", "Differential expression", ["limma (microarray)", "DESeq2 (RNA-seq)"], tsize=7.6)
box(2.35, 1.95, 1.6, 0.6, "#E6F5EE", "#7DB8A0", "Ranked GSEA", ["GO BP pathway enrichment"])
box(2.35, 0.35, 1.6, 1.3, "#FCEBD2", "#D9A760", "Cross-dataset integration", ["signature matrix", "candidate prioritization", "regulatory-signature scoring", "secretome / surfaceome"], lsize=6.8, tsize=7.3)
arrow(3.15, 2.85, 3.15, 2.58); arrow(3.15, 1.95, 3.15, 1.68)
# merge arrows from left
arrow(1.82, 3.1, 2.33, 3.3); arrow(1.82, 2.2, 2.33, 3.1); arrow(1.82, 0.9, 2.33, 3.0)
# right column
box(4.5, 2.6, 1.9, 1.05, "#E8DDF3", "#A58BC4", "Co-expression analysis", ["WGCNA (lung, tonsil)", "module-trait correlation", "stability and power checks"])
box(4.5, 1.3, 1.9, 0.95, "#F8DCDD", "#D68A8E", "Biological model", ["conserved IFN-stimulated core", "tissue-associated programs"])
box(4.5, 0.15, 1.9, 0.85, "#EEEEEE", "#A0A0A0", "Robustness checks", ["study-level scoring, permutation", "null, orthology verification"], lsize=6.6)
arrow(3.97, 1.0, 4.48, 1.75); arrow(3.97, 3.2, 4.48, 3.1); arrow(5.45, 2.58, 5.45, 2.27)
fig.savefig(os.path.join(OUT, "Figure_1_workflow_schematic_print.png"), dpi=600, facecolor="white")
