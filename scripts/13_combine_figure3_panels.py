#!/usr/bin/env python3
"""
13_combine_figure3_panels.py

Combines the two regenerated HUVEC-signature heatmap panels (Lung, Tonsil)
into the final two-panel Figure 3 used in the manuscript, with bold "A"/"B"
panel labels, matching the layout of the original Figure_3_tissue_signature_panels.png.

Run this AFTER re-running GSE310471_DESeq2_analysis_png_only.R (which
regenerates Heatmap_HUVEC_signature_Lung.png and
Heatmap_HUVEC_signature_Tonsil.png with the revised, larger fonts).

Usage (from the Nipah_transcriptomics project root):
    python3 13_combine_figure3_panels.py
"""

import os
from PIL import Image, ImageDraw, ImageFont

PROJECT_DIR = "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
LUNG_PNG = os.path.join(PROJECT_DIR, "GSE310471", "results", "figures", "Heatmap_HUVEC_signature_Lung.png")
TONSIL_PNG = os.path.join(PROJECT_DIR, "GSE310471", "results", "figures", "Heatmap_HUVEC_signature_Tonsil.png")

OUT_DIR = os.path.join(PROJECT_DIR, "Revised_manuscript", "02_figures")
OUT_LUNG = os.path.join(OUT_DIR, "Figure_3_lung_HUVEC_signature.png")
OUT_TONSIL = os.path.join(OUT_DIR, "Figure_3_tonsil_HUVEC_signature.png")
OUT_COMBINED = os.path.join(OUT_DIR, "Figure_3_tissue_signature_panels.png")

LABEL_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/local/lib/python3.10/dist-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSans-Bold.ttf",
]

GAP_PX = 60          # horizontal gap between the two panels
MARGIN_PX = 20        # outer margin around the combined figure
LABEL_MARGIN = 30     # inset of the A/B label from each panel's top-left corner


def get_font(size):
    for path in LABEL_FONT_CANDIDATES:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def main():
    lung = Image.open(LUNG_PNG).convert("RGB")
    tonsil = Image.open(TONSIL_PNG).convert("RGB")

    # Save (possibly updated) individual panels into the manuscript figures folder
    lung.save(OUT_LUNG)
    tonsil.save(OUT_TONSIL)

    h = max(lung.height, tonsil.height)
    w = lung.width + GAP_PX + tonsil.width

    canvas = Image.new("RGB", (w + 2 * MARGIN_PX, h + 2 * MARGIN_PX), "white")
    canvas.paste(lung, (MARGIN_PX, MARGIN_PX))
    canvas.paste(tonsil, (MARGIN_PX + lung.width + GAP_PX, MARGIN_PX))

    draw = ImageDraw.Draw(canvas)
    label_size = max(48, int(min(lung.width, tonsil.width) * 0.035))
    font = get_font(label_size)

    draw.text((MARGIN_PX + LABEL_MARGIN, MARGIN_PX + LABEL_MARGIN), "A", fill="black", font=font)
    draw.text((MARGIN_PX + lung.width + GAP_PX + LABEL_MARGIN, MARGIN_PX + LABEL_MARGIN), "B", fill="black", font=font)

    canvas.save(OUT_COMBINED)
    print("Wrote:")
    print(" ", OUT_LUNG)
    print(" ", OUT_TONSIL)
    print(" ", OUT_COMBINED, canvas.size)


if __name__ == "__main__":
    main()
