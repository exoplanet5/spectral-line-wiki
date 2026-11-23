#!/bin/zsh
# usage: ocr_pdf.sh <pdf> <out.txt> <workdir> [dpi]
# Renders every page of a scanned PDF to PNG (PyMuPDF, run isolated) and OCRs it with the
# macOS Vision framework (compiled from ocr_vision.swift).
set -e
PDF=$1; OUT=$2; WORK=$3; DPI=${4:-200}
HERE=${0:a:h}
mkdir -p $WORK
OCR=$WORK/../ocr_vision
[[ -x $OCR ]] || swiftc -O -o $OCR $HERE/ocr_vision.swift
~/.venvs/astro313/bin/python -I -c "
import fitz,sys
d=fitz.open(sys.argv[1])
for i in range(d.page_count):
    d[i].get_pixmap(dpi=int(sys.argv[3]), colorspace=fitz.csGRAY).save(sys.argv[2]+f'/p{i+1:03d}.png')
" $PDF $WORK $DPI
$OCR $OUT $WORK/p*.png
echo "OCR done: $OUT"
