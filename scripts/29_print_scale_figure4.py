#!/usr/bin/env python3
"""
29_print_scale_figure4.py

Rebuilds Figure 4 as two print-scale panels (6.5 in wide, 600 dpi, 6.5-8.5 pt text):
  Figure_4A_wgcna_module_membership_print.png : lung tan (ME12) and tonsil blue (ME2) module kME panels
  Figure_4B_regulator_axis_heatmap_print.png  : regulatory-signature scores across contrasts
Inputs (advanced_analyses/tables): focused_wgcna_candidate_modules.csv, tf_upstream_regulator_signature_scores.csv
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram

PROJECT = os.environ.get("NIPAH_PROJECT_DIR", "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics")
TAB = os.path.join(PROJECT, "advanced_analyses", "tables")
OUT = os.path.join(PROJECT, "Revised_manuscript_R2", "02_figures")
os.makedirs(OUT, exist_ok=True)
DPI = 600
plt.rcParams.update({"font.family": "DejaVu Sans"})

# ------------------------------------------------------------------ Figure 4A
m = pd.read_csv(os.path.join(TAB, "focused_wgcna_candidate_modules.csv"))
panels = [("Lung antiviral IFN module", "Lung antiviral IFN module (tan, ME12)",
           "High-membership antiviral genes in the infection/progression-associated module"),
          ("Tonsil complement/coagulation module", "Tonsil complement/coagulation module (blue, ME2)",
           "High-membership disease genes in the infection/progression-associated module")]
colors = {"Conserved antiviral/IFN core": "#2B7BB9",
          "Mainly in vivo immune-associated module": "#7B3294",
          "In vivo complement/coagulation disease module": "#D95F02"}
legend_names = {"Conserved antiviral/IFN core": "Conserved antiviral/IFN core",
                "Mainly in vivo immune-associated module": "In vivo IFN/progression support",
                "In vivo complement/coagulation disease module": "Complement/coagulation disease module"}
W = 6.5
rh = 0.150
n1 = (m.discovery_panel == panels[0][0]).sum(); n2 = (m.discovery_panel == panels[1][0]).sum()
H = 0.55 + n1 * rh + 0.55 + 0.45 + n2 * rh + 0.55 + 0.55
fig = plt.figure(figsize=(W, H), dpi=DPI)
y_cursor = H - 0.08
for pi, (key, title, sub) in enumerate(panels):
    t = m[m.discovery_panel == key].sort_values("kME", ascending=True).reset_index(drop=True)
    n = len(t)
    ph = n * rh
    y_title = y_cursor
    fig.text(0.05, y_title / H, title, fontsize=8.5, fontweight="bold", va="top")
    fig.text(0.05, (y_title - 0.19) / H, sub, fontsize=7, color="#444444", va="top")
    ax = fig.add_axes([0.14, (y_title - 0.45 - ph) / H, 0.83, ph / H])
    yp = np.arange(n)
    ax.hlines(yp, 0.78, t.kME, color="#B8B8B8", lw=0.8, zorder=1)
    ax.scatter(t.kME, yp, s=t.priority_score * 2.0 + 6, c=t.manuscript_module.map(colors),
               edgecolor="#222222", lw=0.4, zorder=3)
    ax.axvline(0.80, color="#555555", ls="--", lw=0.7)
    ax.set_xlim(0.78, 1.0); ax.set_ylim(-0.7, n - 0.3)
    ax.set_yticks(yp); ax.set_yticklabels(t.gene, fontsize=7, fontweight="bold")
    ax.set_xticks([0.80, 0.85, 0.90, 0.95, 1.00]); ax.set_xticklabels(["0.80", "0.85", "0.90", "0.95", "1.00"], fontsize=7)
    ax.tick_params(length=2, width=0.5, pad=1.5)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_linewidth(0.5)
    ax.set_xlabel("|kME| module membership", fontsize=7.5, fontweight="bold", labelpad=2)
    r1, r2 = t.infected_module_cor.iloc[0], t.dpi_module_cor.iloc[0]
    ax.text(0.795, n - 0.55, f"module-infection r = {r1:.2f}\nmodule-DPI r = {r2:.2f}", fontsize=6.5, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#333333", lw=0.5), linespacing=1.4)
    y_cursor = y_title - 0.45 - ph - 0.55 - 0.1
fig.text(0.008, (H - 0.06) / H, "A", fontsize=15, fontweight="bold", va="top")
# legend
used = []
for k in ["Conserved antiviral/IFN core", "Mainly in vivo immune-associated module", "In vivo complement/coagulation disease module"]:
    used.append(Line2D([], [], marker="o", ls="", mfc=colors[k], mec="#222222", mew=0.4, ms=6, label=legend_names[k]))
fig.legend(handles=used, loc="lower left", bbox_to_anchor=(0.01, 0.0), ncol=3, fontsize=6.8, frameon=False,
           handletextpad=0.2, columnspacing=1.2)
fig.text(0.01, 0.30 / H, "Point size is proportional to the cross-dataset priority score; dashed line, |kME| = 0.80.",
         fontsize=6.5, color="#444444", va="center")
fig.savefig(os.path.join(OUT, "Figure_4A_wgcna_module_membership_print.png"), dpi=DPI, facecolor="white")
plt.close(fig)

# ------------------------------------------------------------------ Figure 4B
d = pd.read_csv(os.path.join(TAB, "tf_upstream_regulator_signature_scores.csv"))
mat = d.pivot(index="regulator", columns="contrast", values="score")
cols = ["GSE310471_Lung_3DPI_vs_baseline", "GSE310471_Lung_4DPI_vs_baseline", "GSE310471_Lung_5DPI_vs_baseline",
        "GSE310471_Tonsil_3DPI_vs_baseline", "GSE310471_Tonsil_4DPI_vs_baseline", "GSE310471_Tonsil_5DPI_vs_baseline",
        "GSE32902_HUVEC_NiV_vs_Mock", "GSE33133_HUVEC_NiV_vs_Mock", "GSE33133_HUVEC_NiVdC_vs_Mock",
        "GSE33133_HUVEC_NiVdC_vs_NiV"]
clab = ["Lung 3 DPI", "Lung 4 DPI", "Lung 5 DPI", "Tonsil 3 DPI", "Tonsil 4 DPI", "Tonsil 5 DPI",
        "HUVEC NiV vs Mock (GSE32902)", "HUVEC NiV vs Mock (GSE33133)*", "HUVEC NiV-dC vs Mock", "HUVEC NiV-dC vs NiV"]
mat = mat[cols]
dn = dendrogram(linkage(mat.values, "complete"), no_plot=True)
rows = [mat.index[i] for i in dn["leaves"]]
mat = mat.loc[rows]
n = len(rows)
W, H = 6.5, 3.95
fig = plt.figure(figsize=(W, H), dpi=DPI)
x_d, w_d = 0.05, 0.38
x_h, cw = 0.48, 0.325
rh2 = 0.30
y_h = 1.40
axd = fig.add_axes([x_d / W, y_h / H, w_d / W, n * rh2 / H])
for ic, dc in zip(dn["icoord"], dn["dcoord"]):
    axd.plot(dc, [(y - 5) / 10 + 0.5 for y in ic], color="k", lw=0.6)
axd.set_ylim(n, 0); axd.invert_xaxis(); axd.axis("off")
axh = fig.add_axes([x_h / W, y_h / H, cw * 10 / W, n * rh2 / H])
vmin, vmax = -2.2, 14.4
norm = TwoSlopeNorm(vmin=vmin, vcenter=6.1, vmax=vmax)
im = axh.imshow(mat.values, aspect="auto", cmap="RdBu_r", norm=norm)
for i in range(n):
    for j in range(10):
        v = mat.values[i, j]
        axh.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5, color="white" if (v > 11 or v < 0.3) else "black")
axh.set_xticks(range(10))
axh.set_xticklabels(clab, rotation=45, ha="right", rotation_mode="anchor", fontsize=6.8)
axh.yaxis.tick_right(); axh.set_yticks(range(n)); axh.set_yticklabels(rows, fontsize=7)
axh.tick_params(length=0, pad=2)
for s in axh.spines.values():
    s.set_visible(False)
fig.text(0.05, (H - 0.1) / H, "Regulatory-signature scores across Nipah contrasts", fontsize=8.5, fontweight="bold", va="top")
cax = fig.add_axes([x_h / W, (H - 0.52) / H, 1.6 / W, 0.09 / H])
cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.outline.set_visible(False)
cb.set_ticks([-2, 0, 2, 4, 6, 8, 10, 12, 14]); cb.ax.tick_params(labelsize=6.5, length=1.5, pad=1.5)
fig.text((x_h + 1.7) / W, (H - 0.475) / H, "Regulatory-signature score", fontsize=6.8, va="center")
fig.text(0.012, 0.09 / H, "* GSE33133 NiV vs Mock is shown for reference and excluded from axis means because it duplicates the GSE32902 sample records.",
         fontsize=6, color="#333333", va="center")
fig.text(0.008, (H - 0.06) / H, "B", fontsize=15, fontweight="bold", va="top")
fig.savefig(os.path.join(OUT, "Figure_4B_regulator_axis_heatmap_print.png"), dpi=DPI, facecolor="white")
print(rows)
