#!/bin/bash
# Same fetch function as fetch_hitran_airglow.sh, only the O2 minor isotopologues (37, 38)
OUT="$1"; mkdir -p "$OUT"
fetch () {
  url="https://hitran.org/lbl/api?iso_ids_list=$1&numin=$2&numax=$3"
  curl -s -L -m 600 -A "Mozilla/5.0 (research script)" "$url" -o "$OUT/$4.par"; rc=$?
  printf "%s\n%s\ncurl_rc=%s\n" "$url" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" > "$OUT/$4.url"
  echo "$4: rc=$rc $(wc -l < "$OUT/$4.par") lines"
}
for iso in 37 38; do
  fetch $iso 6000 6600 O2_iso${iso}_6000-6600cm
  fetch $iso 7600 8200 O2_iso${iso}_7600-8200cm
  fetch $iso 11300 11700 O2_iso${iso}_11300-11700cm
  fetch $iso 12800 13250 O2_iso${iso}_12800-13250cm
done
