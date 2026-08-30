#!/usr/bin/env python3
"""
14_combine_figure2_panels.py

Combines the regenerated integrated heatmap (panel A) and conserved-signature
lollipop plot (panel B) into the two-panel Figure 2 used in the manuscript,
with bold "A"/"B" panel labels.

Run this AFTER re-running integrated_cross_dataset_signature_heatmap.R (which
regenerates integrated_HUVEC_in_vivo_signature_heatmap.png and
integrated_conserved_signature_lollipop.png with the GSE33133 duplicate
excluded from recurrence/priority-score counting; see Response to Reviewers R2-1).

NOTE on the panel-A title: pheatmap's own `main =` title in
integrated_HUVEC_in_vivo_signature_heatmap.png is clipped at the top edge of
the PNG canvas (the title row overflows the fixed-size png() device before
pheatmap draws it, a known pheatmap layout quirk when a `main` title is
combined with multi-row column annotations). Rather than depend on getting
pheatmap's internal spacing exactly right, this script crops off that broken
title band entirely (the heatmap content itself, including the annotation
colour bars and all gene rows, is never clipped) and draws a clean
replacement title of its own above the cropped heatmap. This is robust to
any future re-run of the R script and does not require pixel-perfect canvas
sizing on the R side.

Usage (from the Nipah_transcriptomics project root):
    python3 14_combine_figure2_panels.py
"""

import os
from PIL import Image, ImageDraw, ImageFont

PROJECT_DIR = "D:/Postdoc_Data/Vorolgia/Nipah_transcriptomics"
HEATMAP_PNG = os.path.join(PROJECT_DIR, "integrated_results", "figures", "integrated_HUVEC_in_vivo_signature_heatmap.png")
LOLLIPOP_PNG = os.path.join(PROJECT_DIR, "integrated_results", "figures", "integrated_conserved_signature_lollipop.png")

OUT_DIR = os.path.join(PROJECT_DIR, "Revised_manuscript", "02_figures")
OUT_HEATMAP = os.path.join(OUT_DIR, "Figure_2_integrated_signature_heatmap.png")
OUT_LOLLIPOP = os.path.join(OUT_DIR, "Figure_2b_conserved_signature_lollipop.png")
OUT_COMBINED = os.path.join(OUT_DIR, "Figure_2_integrated_signature_panels.png")

LABEL_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/local/lib/python3.10/dist-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSans-Bold.ttf",
]

TARGET_WIDTH = 3400   # common width both panels are resized to
LABEL_MARGIN = 30

# Panel-A title text, redrawn cleanly rather than relying on pheatmap's clipped one.
PANEL_A_TITLE = "Integrated Nipah host-response signature: HUVEC and in vivo tissues"
TITLE_FONT_SIZE = 40
TITLE_BAR_HEIGHT = 150  # white band above the (cropped) heatmap holding the "A" label + title

# Row (in the ORIGINAL, un-resized heatmap PNG's pixel coordinates) where pheatmap's
# column-annotation colour bars begin, i.e. everything above this row is the broken/
# clipped title band and gets discarded. Re-check this if the R script's fontsize,
# cellheight, or annotation rows change materially.
HEATMAP_TITLE_CROP_PX = 56  # re-measured after R2-11 fontsize bump (was 59 at fontsize=9/8/8)


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
    heatmap_raw = Image.open(HEATMAP_PNG).convert("RGB")
    lollipop = Image.open(LOLLIPOP_PNG).convert("RGB")

    # Drop the clipped pheatmap title band; keep everything from the colour bars down.
    w, h = heatmap_raw.size
    heatmap_cropped = heatmap_raw.crop((0, HEATMAP_TITLE_CROP_PX, w, h))

    heatmap_r = resize_to_width(heatmap_cropped, TARGET_WIDTH)
    lollipop_r = resize_to_width(lollipop, TARGET_WIDTH)

    # Build panel A = fresh title bar + cropped/resized heatmap
    panel_a = Image.new("RGB", (TARGET_WIDTH, TITLE_BAR_HEIGHT + heatmap_r.height), "white")
    panel_a.paste(heatmap_r, (0, TITLE_BAR_HEIGHT))

    draw_a = ImageDraw.Draw(panel_a)
    title_font = get_font(TITLE_FONT_SIZE)
    title_bbox = draw_a.textbbox((0, 0), PANEL_A_TITLE, font=title_font)
    title_w = title_bbox[2] - title_bbox[0]
    title_h = title_bbox[3] - title_bbox[1]
    title_x = (TARGET_WIDTH - title_w) // 2
    title_y = (TITLE_BAR_HEIGHT - title_h) // 2 - title_bbox[1]
    draw_a.text((title_x, title_y), PANEL_A_TITLE, fill="black", font=title_font)

    # Save (possibly updated) individual panels into the manuscript figures folder
    panel_a.save(OUT_HEATMAP)
    lollipop_r.save(OUT_LOLLIPOP)

    total_h = panel_a.height + lollipop_r.height
    canvas = Image.new("RGB", (TARGET_WIDTH, total_h), "white")
    canvas.paste(panel_a, (0, 0))
    canvas.paste(lollipop_r, (0, panel_a.height))

    draw = ImageDraw.Draw(canvas)
    label_size = max(48, int(TARGET_WIDTH * 0.022))
    font = get_font(label_size)

    draw.text((LABEL_MARGIN, LABEL_MARGIN), "A", fill="black", font=font)
    draw.text((LABEL_MARGIN, panel_a.height + LABEL_MARGIN), "B", fill="black", font=font)

    canvas.save(OUT_COMBINED)
    print("Wrote:")
    print(" ", OUT_HEATMAP, panel_a.size)
    print(" ", OUT_LOLLIPOP, lollipop_r.size)
    print(" ", OUT_COMBINED, canvas.size)


if __name__ == "__main__":
    main()
