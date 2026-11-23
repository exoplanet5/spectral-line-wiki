#!/bin/zsh
# Download selected arXiv PDFs (untrusted) into OUTDIR and convert them to text.
# Usage: fetch_arxiv_pdfs.sh OUTDIR id1 id2 ...
set -u
outdir=$1; shift
mkdir -p "$outdir"
for id in "$@"; do
  fn=${id//\//_}
  if [[ ! -s "$outdir/$fn.pdf" ]]; then
    curl -sS -L -m 120 -A "spectral-line-wiki-research/0.1" -o "$outdir/$fn.pdf" "https://arxiv.org/pdf/$id" || echo "FAIL $id"
    sleep 3
  fi
  pdftotext -layout "$outdir/$fn.pdf" "$outdir/$fn.txt" 2>/dev/null && echo "ok $id $(wc -c < "$outdir/$fn.txt") chars"
done
