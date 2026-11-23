#!/bin/bash
# Fill remaining O2 windows (0.74-2.5 um) for isotopologues 36 (16O2), 37 (16O18O), 38 (16O17O).
# Usage: fetch_hitran_o2_fill.sh <outdir>. Empty ranges return HTTP 404 (saved as .nodata_404.html).
OUT="$1"; mkdir -p "$OUT"
fetch () {
  url="https://hitran.org/lbl/api?iso_ids_list=$1&numin=$2&numax=$3"
  curl -s -L -m 600 -A "Mozilla/5.0 (research script)" "$url" -o "$OUT/$4.par"; rc=$?
  printf "%s\n%s\ncurl_rc=%s\n" "$url" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" > "$OUT/$4.url"
  if head -c 20 "$OUT/$4.par" | grep -q "<"; then mv "$OUT/$4.par" "$OUT/$4.nodata_404.html"; echo "$4: no data"; else echo "$4: rc=$rc $(wc -l < "$OUT/$4.par") lines"; fi
}
for iso in 36 37 38; do
  n=O2_16; [ $iso = 37 ] && n=O2_iso37; [ $iso = 38 ] && n=O2_iso38
  fetch $iso 3990 6000 ${n}_3990-6000cm
  fetch $iso 6600 7600 ${n}_6600-7600cm
  fetch $iso 8200 11300 ${n}_8200-11300cm
  fetch $iso 11700 12800 ${n}_11700-12800cm
  fetch $iso 13250 13530 ${n}_13250-13530cm
done
