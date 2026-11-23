#!/bin/bash
# Fetch HITRAN line-by-line .par records (160-char HITRAN2004+ format) from HITRANonline
# public API (the endpoint used by HAPI's queryHITRAN "old-fashioned .par search").
# Usage: fetch_hitran_airglow.sh <outdir>
# Global isotopologue IDs (HAPI ISO table): 48 = (16O)H, 36 = (16O)2
OUT="$1"; mkdir -p "$OUT"
fetch () { # iso numin numax name
  url="https://hitran.org/lbl/api?iso_ids_list=$1&numin=$2&numax=$3"
  for try in 1 2 3; do
    curl -s -L -m 600 -A "Mozilla/5.0 (research script)" "$url" -o "$OUT/$4.par"; rc=$?
    [ $rc -eq 0 ] && break
    echo "  retry $try (curl rc=$rc)"; sleep 5
  done
  printf "%s\n%s\ncurl_rc=%s\n" "$url" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" > "$OUT/$4.url"
  echo "$4: rc=$rc $(wc -l < "$OUT/$4.par") lines"
}
rm -f "$OUT"/OH16_3800-11400cm.*
for a in 3800 4800 5800 6800 7800 8800 9800 10800; do
  b=$((a+1000)); [ $b -gt 11400 ] && b=11400
  fetch 48 $a $b OH16_${a}-${b}cm
done
fetch 36 6000 6600 O2_16_6000-6600cm
fetch 36 7600 8200 O2_16_7600-8200cm
fetch 36 11300 11700 O2_16_11300-11700cm
fetch 36 12800 13250 O2_16_12800-13250cm
# O2 minor isotopologues: 37 = (16O)(18O), 38 = (16O)(17O)
for iso in 37 38; do
  fetch $iso 6000 6600 O2_iso${iso}_6000-6600cm
  fetch $iso 7600 8200 O2_iso${iso}_7600-8200cm
  fetch $iso 11300 11700 O2_iso${iso}_11300-11700cm
  fetch $iso 12800 13250 O2_iso${iso}_12800-13250cm
done
