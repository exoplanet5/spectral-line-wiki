#!/bin/bash
# 18OH (HITRAN global iso id 49) in the NIR airglow range. Usage: fetch_hitran_18OH.sh <outdir>
OUT="$1"; mkdir -p "$OUT"
url="https://hitran.org/lbl/api?iso_ids_list=49&numin=3800&numax=11400"
curl -s -L -m 600 -A "Mozilla/5.0 (research script)" "$url" -o "$OUT/OH18_3800-11400cm.par"; rc=$?
printf "%s\n%s\ncurl_rc=%s\n" "$url" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rc" > "$OUT/OH18_3800-11400cm.url"
echo "rc=$rc $(wc -l < "$OUT/OH18_3800-11400cm.par") lines"
