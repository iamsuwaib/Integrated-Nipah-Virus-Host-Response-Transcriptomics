#!/usr/bin/env python3
"""
28_print_scale_figures.py

Rebuilds Figure 2A (integrated signature heatmap) and Figure 2B (conserved-signature
lollipop) directly from integrated_results/tables/integrated_signature_long.csv at true
print size (6.5 in wide, 600 dpi, 6.5-8 pt text), so that all gene and contrast labels
remain legible when the figures are placed at page width in the manuscript.
Scoring conventions are identical to 04_integrated_cross_dataset_signature_heatmap.R:
log2FC capped at +/-5, grey = gene not detected, GSE33133 NiV vs Mock shown with an
asterisk and excluded from recurrence counting.

Usage:  python3 28_print_scale_figures.py   (NIPAH_PROJECT_DIR optional)
Outputs (Revised_manuscript_R2/02_figures):
  Figure_2A_integrated_signature_heatmap_print.png
  Figure_2B_conserved_signature_lollipop_print.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib as mpl

PROJECT = os.environ.get("NIPAH_PROJECT_DIR", "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics")
LONG = os.path.join(PROJECT, "integrated_results", "tables", "integrated_signature_long.csv")
OUT = os.path.join(PROJECT, "Revised_manuscript_R2", "02_figures")
os.makedirs(OUT, exist_ok=True)
DPI = 600
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})

d = pd.read_csv(LONG)
contrast_order = [
    "GSE32902_HUVEC_NiV_vs_Mock", "GSE33133_HUVEC_NiV_vs_Mock",
    "GSE33133_HUVEC_NiVdC_vs_Mock", "GSE33133_HUVEC_NiVdC_vs_NiV",
    "GSE310471_Lung_3DPI_vs_baseline", "GSE310471_Lung_4DPI_vs_baseline",
    "GSE310471_Lung_5DPI_vs_baseline", "GSE310471_Tonsil_3DPI_vs_baseline",
    "GSE310471_Tonsil_4DPI_vs_baseline", "GSE310471_Tonsil_5DPI_vs_baseline"]
labels = ["NiV vs Mock (32902)", "NiV vs Mock (33133)*", "NiV-dC vs Mock", "NiV-dC vs NiV",
          "Lung 3 DPI", "Lung 4 DPI", "Lung 5 DPI", "Tonsil 3 DPI", "Tonsil 4 DPI", "Tonsil 5 DPI"]
cat_order = ["ISG_antiviral", "cytokine_chemokine", "endothelial_activation", "complement", "coagulation"]
cat_names = {"ISG_antiviral": "ISG / antiviral", "cytokine_chemokine": "Cytokine / chemokine",
             "endothelial_activation": "Endothelial activation", "complement": "Complement",
             "coagulation": "Coagulation"}
cat_col = {"ISG_antiviral": "#4E79A7", "cytokine_chemokine": "#E15759", "endothelial_activation": "#F28E2B",
           "complement": "#59A14F", "coagulation": "#B07AA1"}
meta = d.drop_duplicates("contrast").set_index("contrast").loc[contrast_order]
model_col = {"Human endothelial cell": "#A6775F", "African green monkey in vivo": "#E8C23A"}
tissue_col = {"HUVEC": "#76BDB5", "Lung": "#F28E2B", "Tonsil": "#B07AA1"}
dataset_col = {"GSE32902": "#4E79A7", "GSE33133": "#59A14F", "GSE310471": "#E15759"}
model_short = {"Human endothelial cell": "Human endothelial cell", "African green monkey in vivo": "African green monkey in vivo"}

# ------------------------------------------------------------------ Figure 2A
genes = (d.drop_duplicates("gene")[["gene", "category"]]
         .assign(c=lambda x: x.category.map({k: i for i, k in enumerate(cat_order)}))
         .sort_values(["c", "gene"]))
mat = d.pivot_table(index="gene", columns="contrast", values="log2FC", aggfunc="first")
mat = mat.reindex(index=genes.gene, columns=contrast_order)
capped = mat.clip(-5, 5)

W, H = 6.5, 8.15
fig = plt.figure(figsize=(W, H), dpi=DPI)
def ax_in(x, y, w, h):  # inches from left/bottom
    return fig.add_axes([x / W, y / H, w / W, h / H])

x0, cw = 0.98, 0.325
cell_w = cw * 10
rh = 0.1135
n = len(genes)
cell_h = rh * n
y_top = H - 0.18
bar_h = 0.105
# annotation bars
for i, (lab, cols, key) in enumerate([("Model", model_col, "model"), ("Tissue", tissue_col, "tissue"),
                                       ("Dataset", dataset_col, "dataset")]):
    ax = ax_in(x0, y_top - (i + 1) * bar_h, cell_w, bar_h)
    ax.imshow([[list(matplotlib.colors.to_rgb(cols[v])) for v in meta[key]]], aspect="auto")
    ax.set_axis_off()
    fig.text((x0 + cell_w + 0.06) / W, (y_top - (i + 0.5) * bar_h) / H, lab, fontsize=6.5, fontweight="bold",
             va="center", ha="left")
y_cells_top = y_top - 3 * bar_h - 0.04
axh = ax_in(x0, y_cells_top - cell_h, cell_w, cell_h)
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#E5E5E5")
axh.imshow(np.ma.masked_invalid(capped.values), cmap=cmap, norm=Normalize(-5, 5), aspect="auto",
           interpolation="nearest")
axh.set_xticks(np.arange(10))
axh.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=7)
axh.xaxis.set_ticks_position("bottom")
axh.set_yticks(np.arange(n))
axh.yaxis.tick_right()
axh.set_yticklabels(genes.gene, fontsize=6.8)
axh.tick_params(axis="both", length=1.5, width=0.4, pad=1.5)
for s in axh.spines.values():
    s.set_visible(False)
# vertical separators between model/dataset groups (thin white)
for xv in [0.5, 3.5]:
    axh.axvline(xv, color="white", lw=1.2)
# category strip + separators
axc = ax_in(x0 - 0.12, y_cells_top - cell_h, 0.09, cell_h)
axc.imshow([[list(matplotlib.colors.to_rgb(cat_col[c]))] for c in genes.category], aspect="auto")
axc.set_axis_off()
cum = 0
for c in cat_order:
    k = (genes.category == c).sum()
    cum += k
    if cum < n:
        axh.axhline(cum - 0.5, color="white", lw=1.0)

# legends (right column)
lx = 5.0
yy = H - 0.2
def legend_block(title, items, y):
    fig.text(lx / W, y / H, title, fontsize=7, fontweight="bold", va="top")
    y -= 0.16
    for name, col in items:
        axp = ax_in(lx, y - 0.1, 0.12, 0.1)
        axp.add_patch(plt.Rectangle((0, 0), 1, 1, color=col)); axp.set_xlim(0, 1); axp.set_ylim(0, 1); axp.set_axis_off()
        fig.text((lx + 0.17) / W, (y - 0.05) / H, name, fontsize=6.5, va="center")
        y -= 0.145
    return y - 0.06
# colorbar
cb_ax = ax_in(lx, yy - 0.95, 0.13, 0.85)
grad = np.linspace(5, -5, 256)[:, None]
cb_ax.imshow(grad, cmap="RdBu_r", norm=Normalize(-5, 5), aspect="auto", extent=[0, 1, -5, 5])
cb_ax.set_xticks([]); cb_ax.yaxis.tick_right()
cb_ax.set_yticks([-4, -2, 0, 2, 4]); cb_ax.tick_params(labelsize=6.5, length=1.5, width=0.4, pad=1.5)
for s in cb_ax.spines.values():
    s.set_linewidth(0.4)
fig.text((lx + 0.0) / W, (yy + 0.0) / H, "log2 fold change", fontsize=7, fontweight="bold", va="top")
# move colorbar below title
cb_ax.set_position([lx / W, (yy - 1.12) / H, 0.13 / W, 0.9 / H])
y = yy - 1.28
y = legend_block("Model", [("Human endothelial cell", model_col["Human endothelial cell"]), ("African green monkey", model_col["African green monkey in vivo"])], y)
y = legend_block("Tissue", list(tissue_col.items()), y)
y = legend_block("Dataset", list(dataset_col.items()), y)
y = legend_block("Gene category", [(cat_names[c], cat_col[c]) for c in cat_order], y)
axp = ax_in(lx, y - 0.1, 0.12, 0.1)
axp.add_patch(plt.Rectangle((0, 0), 1, 1, color="#E5E5E5")); axp.set_axis_off(); axp.set_xlim(0, 1); axp.set_ylim(0, 1)
fig.text((lx + 0.17) / W, (y - 0.05) / H, "Not detected", fontsize=6.5, va="center")
fig.text(0.012, 0.992, "A", fontsize=15, fontweight="bold", va="top")
fig.savefig(os.path.join(OUT, "Figure_2A_integrated_signature_heatmap_print.png"), dpi=DPI, facecolor="white")
plt.close(fig)

# ------------------------------------------------------------------ Figure 2B
s = d[d.scoring_included]
g = s.groupby(["gene", "category"]).apply(lambda t: pd.Series({
    "n_up": int((t.direction == "up").sum()),
    "n_huvec": int(((t.tissue == "HUVEC") & (t.direction == "up")).sum()),
    "n_vivo": int((t.tissue.isin(["Lung", "Tonsil"]) & (t.direction == "up")).sum()),
    "mean_fc": t.log2FC.mean(), "min_padj": t.padj.min()}), include_groups=False).reset_index()
g = g.sort_values(["n_up", "n_huvec", "n_vivo", "min_padj"], ascending=[False, False, False, True])
g = g[g.n_up > 0].head(30).sort_values(["n_up", "mean_fc"], ascending=[True, True]).reset_index(drop=True)
def call(r):
    if r.n_huvec >= 1 and r.n_vivo >= 2: return "conserved"
    if r.n_huvec >= 1: return "huvec"
    if r.n_vivo >= 2: return "vivo"
    return "none"
g["call"] = g.apply(call, axis=1)
ccol = {"conserved": "#762A83", "vivo": "#E66101", "huvec": "#1F78B4", "none": "#999999"}
cname = {"conserved": "Conserved in HUVEC and in vivo", "vivo": "Mainly in vivo",
         "huvec": "Mainly HUVEC", "none": "Not strongly conserved"}
W2, H2 = 6.5, 5.3
fig = plt.figure(figsize=(W2, H2), dpi=DPI)
ax = fig.add_axes([0.13, 0.17, 0.84, 0.80])
ypos = np.arange(len(g))
ax.hlines(ypos, 0, g.n_up, color="#B8B8B8", lw=1.0, zorder=1)
ax.scatter(g.n_up, ypos, s=(g.n_vivo * 8 + 14), c=g.call.map(ccol), edgecolor="#333333", lw=0.4, zorder=3)
for yv, r in zip(ypos, g.itertuples()):
    ax.text(r.n_up + 0.22, yv, f"HUVEC: {int(r.n_huvec)} | in vivo: {int(r.n_vivo)}", fontsize=6.5, va="center", color="#222222")
ax.set_yticks(ypos); ax.set_yticklabels(g.gene, fontsize=7)
ax.set_xlim(0, g.n_up.max() + 3.2)
ax.set_xticks(range(0, int(g.n_up.max()) + 2))
ax.tick_params(axis="x", labelsize=7.5); ax.tick_params(length=2, width=0.5)
ax.set_xlabel("Number of significant upregulated contrasts", fontsize=8, fontweight="bold")
ax.set_ylim(-0.7, len(g) - 0.3)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
for sp in ["left", "bottom"]:
    ax.spines[sp].set_linewidth(0.5)
present = [k for k in ["conserved", "vivo", "huvec", "none"] if (g.call == k).any()]
h1 = [Line2D([], [], marker="o", ls="", mfc=ccol[k], mec="#333333", mew=0.4, ms=6, label=cname[k]) for k in present]
l1 = fig.legend(handles=h1, loc="lower left", bbox_to_anchor=(0.03, 0.005), ncol=1, fontsize=7, frameon=False,
                handletextpad=0.3, labelspacing=0.35)
h2 = [Line2D([], [], marker="o", ls="", mfc="white", mec="#333333", mew=0.5,
             ms=np.sqrt(k * 8 + 14), label=str(k)) for k in range(1, 7)]
fig.legend(handles=h2, loc="lower right", bbox_to_anchor=(0.97, 0.005), ncol=6, fontsize=7, frameon=False,
           title="In vivo significant upregulated contrasts", title_fontsize=7, handletextpad=0.2, columnspacing=1.0)
fig.text(0.012, 0.992, "B", fontsize=15, fontweight="bold", va="top")
fig.savefig(os.path.join(OUT, "Figure_2B_conserved_signature_lollipop_print.png"), dpi=DPI, facecolor="white")
print(g[["gene", "n_up", "n_huvec", "n_vivo", "call"]].iloc[::-1].to_string())
