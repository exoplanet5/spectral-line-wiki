#!/usr/bin/env python
"""Make image crops of NBS Monograph 61 (Moore, Minnaert & Houtgast 1966) pages around given
wavelengths so that table entries can be transcribed/verified by eye.

Usage: python mmh_crops.py <pdf> <outdir> page:wl:label [page:wl:label ...]
Uses pdftotext -bbox to locate the OCR token closest to wl on that page, then crops a band of
the rendered page (pdftoppm, 110 dpi) around it.
"""
import os, re, subprocess, sys
from PIL import Image

pdf, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
DPI = 110
for spec in sys.argv[3:]:
    page, wl, label = spec.split(':')
    page = int(page); wl = float(wl)
    html = subprocess.run(['pdftotext', '-f', str(page), '-l', str(page), '-bbox', pdf, '-'],
                          capture_output=True, text=True).stdout
    m = re.search(r'<page width="([\d.]+)" height="([\d.]+)"', html)
    pw, ph = float(m.group(1)), float(m.group(2))
    best = None
    for mm in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', html):
        t = re.match(r'(\d{4}\.\d{1,3})', mm.group(5))
        if not t:
            continue
        d = abs(float(t.group(1)) - wl)
        if best is None or d < best[0]:
            best = (d, float(mm.group(2)), float(mm.group(1)), mm.group(5))
    png = os.path.join(outdir, f'page{page}')
    if not os.path.exists(png + f'-{page}.png') and not os.path.exists(png + f'-{page:03d}.png'):
        subprocess.run(['pdftoppm', '-f', str(page), '-l', str(page), '-r', str(DPI), '-png', pdf, png])
    fn = [f for f in os.listdir(outdir) if f.startswith(f'page{page}-')][0]
    im = Image.open(os.path.join(outdir, fn))
    W, H = im.size
    y = best[1] / ph * H if best else H / 2
    x = best[2] / pw * W if best else 0
    left = 0 if x < W / 2 else W // 2 - 30
    crop = im.crop((0, max(0, int(y) - 110), W, min(H, int(y) + 110)))
    crop.save(os.path.join(outdir, f'crop_{label}_p{page}.png'))
    print(label, page, 'nearest OCR token', best[3] if best else None, 'dist', round(best[0], 3) if best else None)
