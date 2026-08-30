#!/usr/bin/env python3
"""
15_combine_figure4_panels.py

Combines the regenerated focused-WGCNA discovery panel (lung + tonsil, panel A)
and the upstream regulator-axis score heatmap (panel B) into the two-panel
Figure 4 used in the manuscript, with bold "A"/"B" panel labels.

Run this AFTER re-running 04_focused_wgcna_module_discovery_figure.R and
02_tf_upstream_regulator_signature_scoring.R (which regenerate
focused_wgcna_candidate_modules_panel.png and
tf_upstream_regulator_score_heatmap.png with the GSE33133 duplicate excluded
from recurrence/priority-score counting).

Usage (from the Nipah_transcriptomics project root):
    python3 15_combine_figure4_panels.py
"""

import os
from PIL import Image, ImageDraw, ImageFont

PROJECT_DIR = "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
WGCNA_PANEL_PNG = os.path.join(PROJECT_DIR, "advanced_analyses", "figures", "focused_wgcna_candidate_modules_panel.png")
REGULATOR_HEATMAP_PNG = os.path.join(PROJECT_DIR, "advanced_analyses", "figures", "tf_upstream_regulator_score_heatmap.png")

OUT_DIR = os.path.join(PROJECT_DIR, "Revised_manuscript", "02_figures")
OUT_WGCNA = os.path.join(OUT_DIR, "Figure_4_focused_WGCNA_modules.png")
OUT_REGULATOR = os.path.join(OUT_DIR, "Figure_4_regulator_heatmap.png")
OUT_COMBINED = os.path.join(OUT_DIR, "Figure_4_network_regulator_panels.png")

LABEL_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/local/lib/python3.10/dist-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSans-Bold.ttf",
]

TARGET_WIDTH = 3600
LABEL_MARGIN = 30


def get_font(size):
    for path in LABEL_FONT_CANDIDATES:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def resize_to_width(im, width):
    w, h = im.size
    new_h = round(h * (width / w))
    return im.resize((width, new_h), Image.LANCZOS)


def main():
    wgcna = Image.open(WGCNA_PANEL_PNG).convert("RGB")
    regulator = Image.open(REGULATOR_HEATMAP_PNG).convert("RGB")

    wgcna_r = resize_to_width(wgcna, TARGET_WIDTH)
    regulator_r = resize_to_width(regulator, TARGET_WIDTH)

    wgcna_r.save(OUT_WGCNA)
    regulator_r.save(OUT_REGULATOR)

    total_h = wgcna_r.height + regulator_r.height
    canvas = Image.new("RGB", (TARGET_WIDTH, total_h), "white")
    canvas.paste(wgcna_r, (0, 0))
    canvas.paste(regulator_r, (0, wgcna_r.height))

    draw = ImageDraw.Draw(canvas)
    label_size = max(48, int(TARGET_WIDTH * 0.022))
    font = get_font(label_size)

    draw.text((LABEL_MARGIN, LABEL_MARGIN), "A", fill="black", font=font)
    draw.text((LABEL_MARGIN, wgcna_r.height + LABEL_MARGIN), "B", fill="black", font=font)

    canvas.save(OUT_COMBINED)
    print("Wrote:")
    print(" ", OUT_WGCNA, wgcna_r.size)
    print(" ", OUT_REGULATOR, regulator_r.size)
    print(" ", OUT_COMBINED, canvas.size)


if __name__ == "__main__":
    main()
