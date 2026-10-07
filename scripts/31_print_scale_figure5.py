#!/usr/bin/env python3
"""
31_print_scale_figure5.py

Rebuilds Figure 5 at print size (6.5 in wide, 600 dpi) from Supplementary Table S11 (targeted GSEA)
and Supplementary Table S12 (module over-representation analysis):
  Figure_5A_targeted_GSEA_print.png     : NES (colour) and -log10 FDR (size) for 20 curated GO BP terms
  Figure_5B_module_ORA_print.png        : GO BP over-representation in lung tan/ME12 and tonsil blue/ME2
Terms with FDR > 0.25 (GSEA) are not drawn, as in the original figure.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

PROJECT = os.environ.get("NIPAH_PROJECT_DIR", "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics")
XLSX = os.path.join(PROJECT, "Revised_manuscript_R2", "03_supplementary", "Supplementary_Tables_CSBJ.xlsx")
OUT = os.path.join(PROJECT, "Revised_manuscript_R2", "02_figures")
DPI = 600
plt.rcParams.update({"font.family": "DejaVu Sans"})

g = pd.read_excel(XLSX, sheet_name="S11_Targeted_GSEA", header=1)
g = g[g["p.adjust"] <= 0.25].copy()
cols = [("GSE32902", "NiV_vs_Mock", "NiV/mock\n(32902)"),
        ("GSE33133", "NiV_vs_Mock", "NiV/mock\n(33133)*"),
        ("GSE33133", "NiVdC_vs_Mock", "NiV-dC/\nmock"),
        ("GSE33133", "NiVdC_vs_NiV", "NiV-dC/\nNiV"),
        ("GSE310471", "Lung_3DPI_vs_baseline", "Lung\n3 DPI"), ("GSE310471", "Lung_4DPI_vs_baseline", "Lung\n4 DPI"),
        ("GSE310471", "Lung_5DPI_vs_baseline", "Lung\n5 DPI"),
        ("GSE310471", "Tonsil_3DPI_vs_baseline", "Tonsil\n3 DPI"), ("GSE310471", "Tonsil_4DPI_vs_baseline", "Tonsil\n4 DPI"),
        ("GSE310471", "Tonsil_5DPI_vs_baseline", "Tonsil\n5 DPI")]
themes = [("Antiviral/IFN", "Antiviral\nIFN", ["Defense response to virus", "Response to virus", "Response to type I interferon",
                                               "Type I IFN signaling", "Interferon signaling", "Antiviral innate immune response"]),
          ("Complement", "Complement", ["Complement activation", "Complement, classical pathway", "Complement, alternative pathway"]),
          ("Coagulation/hemostasis", "Coagulation/\nhemostasis", ["Regulation of blood coagulation", "Blood coagulation", "Hemostasis", "Fibrinolysis"]),
          ("Vascular/leukocyte", "Vascular/\nleukocyte", ["Vascular process", "Leukocyte migration", "Endothelial cell migration"]),
          ("Metabolic remodeling", "Metabolic\nremodeling", ["Oxidative phosphorylation", "Mitochondrial ATP synthesis", "Aerobic respiration", "Ribosome biogenesis"])]
rows = [(th_short, t) for _, th_short, ts in themes for t in ts]
nrow = len(rows)
W, H = 6.5, 5.45
fig = plt.figure(figsize=(W, H), dpi=DPI)
left, cw, rh = 2.5, 0.295, 0.215
x_gap = {4: 0.10, 7: 0.10}
xs = []; x = left
for j in range(10):
    if j in x_gap: x += x_gap[j]
    xs.append(x + cw / 2); x += cw
top = H - 0.55
bottom = top - nrow * rh
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-3.5, 3.5)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
# theme blocks
y = top
for th, th_short, ts in themes:
    h = len(ts) * rh
    ax.add_patch(plt.Rectangle((0.04, y - h + 0.02), 0.6, h - 0.04, fc="#F1F2F4", ec="#BBBBBB", lw=0.5))
    ax.text(0.34, y - h / 2, th_short, fontsize=6, fontweight="bold", ha="center", va="center", linespacing=1.15)
    for blk in [(0, 3), (4, 6), (7, 9)]:
        x0 = xs[blk[0]] - cw / 2; x1 = xs[blk[1]] + cw / 2
        ax.add_patch(plt.Rectangle((x0, y - h + 0.02), x1 - x0, h - 0.04, fc="white", ec="#444444", lw=0.5, zorder=1))
    y -= h
for blk, name in [((0, 3), "HUVEC"), ((4, 6), "Lung"), ((7, 9), "Tonsil")]:
    x0 = xs[blk[0]] - cw / 2; x1 = xs[blk[1]] + cw / 2
    ax.add_patch(plt.Rectangle((x0, top + 0.03), x1 - x0, 0.24, fc="#F1F2F4", ec="#444444", lw=0.5))
    ax.text((x0 + x1) / 2, top + 0.15, name, fontsize=7.5, fontweight="bold", ha="center", va="center")
for i, (_, t) in enumerate(rows):
    yy = top - (i + 0.5) * rh
    ax.text(left - 0.08, yy, t, fontsize=6.5, ha="right", va="center")
    for j, (ds, ct, lab) in enumerate(cols):
        r = g[(g.dataset == ds) & (g.contrast == ct) & (g.display_term == t)]
        if len(r):
            r = r.iloc[0]
            s = -np.log10(max(r["p.adjust"], 1e-12))
            ax.scatter(xs[j], yy, s=8 + s * 11, c=[cmap(norm(r.NES))], edgecolor="none", alpha=0.95, zorder=3)
for j, (_, _, lab) in enumerate(cols):
    ax.text(xs[j], bottom - 0.06, lab, fontsize=6.2, ha="right", va="top", rotation=45, rotation_mode="anchor")
# legends
lx = 5.85
cax = fig.add_axes([lx / W, (top - 1.05) / H, 0.1 / W, 1.0 / H])
cax.imshow(np.linspace(3.5, -3.5, 256)[:, None], cmap=cmap, aspect="auto", extent=[0, 1, -3.5, 3.5])
cax.set_xticks([]); cax.yaxis.tick_right(); cax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
cax.tick_params(labelsize=6, length=1.5, pad=1.5)
for s in cax.spines.values(): s.set_visible(False)
ax.text(lx - 0.02, top + 0.0, "NES", fontsize=6.8, fontweight="bold", va="bottom")
ax.text(lx - 0.02, top - 1.35, "-log10\n(FDR)", fontsize=6.3, fontweight="bold", va="top")
for k, v in enumerate([2, 5]):
    ax.scatter(lx + 0.05, top - 1.62 - k * 0.2, s=8 + v * 11, c="#555555", edgecolor="none")
    ax.text(lx + 0.18, top - 1.62 - k * 0.2, str(v), fontsize=6, va="center")
ax.text(0.04, H - 0.08, "A", fontsize=15, fontweight="bold", va="top")
ax.text(left - 0.08, H - 0.14, "* GSE33133 NiV vs mock re-lists the GSE32902 sample records.", fontsize=6, color="#444444", ha="left", va="top")
fig.savefig(os.path.join(OUT, "Figure_5A_targeted_GSEA_print.png"), dpi=DPI, facecolor="white")
plt.close(fig)

# ----------------------------------------------------------- Figure 5B
o = pd.read_excel(XLSX, sheet_name="S12_Module_ORA", header=1)
panels = [("Lung tan/ME12: antiviral IFN module", "Lung tan/ME12: antiviral IFN module", "#3A78B5"),
          ("Tonsil blue/ME2: complement/coagulation module", "Tonsil blue/ME2: complement/coagulation-associated module", "#BC2B3D")]
sub = [o[(o.module_panel == p) & (o["p.adjust"] < 0.05)].assign(nl=lambda d: -np.log10(d["p.adjust"])).sort_values("nl", ascending=False) for p, _, _ in panels]
nb = [len(s) for s in sub]
bh = 0.17
W, H = 6.5, 0.1 + 0.28 + nb[0] * bh + 0.42 + 0.28 + nb[1] * bh + 0.75
fig = plt.figure(figsize=(W, H), dpi=DPI)
ymax = max(s.nl.max() for s in sub) * 1.12
ycur = H - 0.08
for k, ((p, title, col), s) in enumerate(zip(panels, sub)):
    ph = nb[k] * bh
    axp = fig.add_axes([2.55 / W, (ycur - 0.28 - ph) / H, 3.7 / W, ph / H])
    yp = np.arange(len(s))[::-1]
    axp.barh(yp, s.nl, height=0.72, color=col, zorder=3)
    for yv, r in zip(yp, s.itertuples()):
        axp.text(r.nl + ymax * 0.01, yv, str(int(r.Count)), fontsize=6.3, va="center", ha="left")
    axp.set_yticks(yp); axp.set_yticklabels(s.display_term, fontsize=6.6)
    axp.set_xlim(0, ymax); axp.set_ylim(-0.6, len(s) - 0.4)
    axp.grid(axis="x", color="#E3E3E3", lw=0.5, zorder=0); axp.set_axisbelow(True)
    axp.tick_params(length=2, width=0.5, labelsize=6.5, pad=1.5)
    for sp in ["top", "right"]: axp.spines[sp].set_visible(False)
    if k == 0: axp.set_xticklabels([])
    else: axp.set_xlabel("-log10(FDR)", fontsize=7.2, fontweight="bold")
    fig.text(0.45 / W, (ycur - 0.06) / H, title, fontsize=7.6, fontweight="bold", va="top")
    ycur -= 0.28 + ph + 0.42
fig.text(0.008, 0.992, "B", fontsize=15, fontweight="bold", va="top")
fig.text(0.45 / W, 0.08 / H, "Numbers at bar ends: module genes annotated to the term; terms with FDR < 0.05 are shown.", fontsize=6, color="#444444")
fig.savefig(os.path.join(OUT, "Figure_5B_module_ORA_print.png"), dpi=DPI, facecolor="white")
print(nb)
