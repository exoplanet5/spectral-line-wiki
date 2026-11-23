#!/usr/bin/env python
"""Extract the composite-transition rows (atomic-weight rows) of Tables 1-11 of
Murphy & Berengut (2014, MNRAS 438, 388; arXiv:1311.2949v3) from the
pdftotext -layout output, and compute the many-multiplet velocity shift for
a given Delta alpha / alpha:   Delta v_i = -2 c (Delta alpha/alpha) q_i / omega_i
(eq. 1 of the paper inverted).

Strings are kept verbatim (wavenumber, vacuum wavelength with uncertainty in
parentheses, delta-v, q with uncertainty).

Usage (input is untrusted download):
    python -I parse_murphy_berengut_2014.py 1311.2949.txt OUT.json
"""
import json
import re
import sys

C = 299_792_458.0

ROW = re.compile(
    r"^\s*(?P<ion>[A-Z][a-z]?\s+[iv]+)?\s*"
    r"(?P<tran>\d{4}(?:\.\d)?)\s+"
    r"(?P<A>\d{2,3}\.\d{2,4})\s+"
    r"(?P<w0>\d{4,5}\.\d+\(\d+\))\S*\s+"
    r"(?P<X>[0-4])\s+"
    r"(?P<lam>\d{4}\.\d+\(\d+\))\s+"
    r"(?P<dv>\d+(?:\.\d+)?)\s+"
    r"(?P<rest>.*)$")
QTAIL = re.compile(r"\s(?P<q>[−-]?\d+)[a-z, ]*\((?P<dq>\d+)\)\s*$")


def main():
    txt = open(sys.argv[1], encoding="utf-8").read().splitlines()
    rows, ion = [], None
    for line in txt:
        m = ROW.match(line)
        if not m:
            continue
        if m.group("ion"):
            ion = " ".join(m.group("ion").split())
        rest = m.group("rest")
        qm = QTAIL.search(rest)
        q = dq = None
        if qm:
            q = int(qm.group("q").replace("−", "-"))
            dq = int(qm.group("dq"))
        w0s = m.group("w0")
        w0 = float(w0s.split("(")[0])
        rec = {
            "ion": ion, "transition": m.group("tran"), "A_composite": m.group("A"),
            "omega0_cm-1": w0s, "lambda0_vac_A": m.group("lam"), "delta_v_lab_m_s": m.group("dv"),
            "q_cm-1": q, "q_unc_cm-1": dq,
        }
        if q is not None:
            # velocity shift for Delta alpha/alpha = 1e-6 (1 ppm)
            rec["dv_for_dalpha_1ppm_m_s"] = round(-2 * C * 1e-6 * q / w0, 3)
            lam = float(m.group("lam").split("(")[0])
            rec["dlambda_for_1ppm_A"] = round(lam * rec["dv_for_dalpha_1ppm_m_s"] / C, 7)
        rec["lab_unc_A_equiv"] = round(float(m.group("lam").split("(")[0]) * float(m.group("dv")) / C, 7)
        rows.append(rec)
    for r in rows:
        print(f"{r['ion']:7s} {r['transition']:7s} w0={r['omega0_cm-1']:20s} lam={r['lambda0_vac_A']:18s} "
              f"dv_lab={r['delta_v_lab_m_s']:5s} q={r['q_cm-1']}({r['q_unc_cm-1']}) "
              f"dv(1ppm)={r.get('dv_for_dalpha_1ppm_m_s')} m/s  dlam(1ppm)={r.get('dlambda_for_1ppm_A')} A")
    json.dump(rows, open(sys.argv[2], "w"), indent=1, ensure_ascii=False)
    print(len(rows), "rows ->", sys.argv[2])


if __name__ == "__main__":
    main()
