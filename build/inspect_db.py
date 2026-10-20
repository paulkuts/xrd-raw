# -*- coding: utf-8 -*-
"""Быстрая инспекция собранного cod_phases.jsonl."""
import json
from pathlib import Path

p = Path("cod_phases.jsonl")
rs = [json.loads(l) for l in open(p, encoding="utf-8")]
print("всего COD-фаз:", len(rs))


def top(r):
    return max(r["peaks"], key=lambda x: x[1])


ca = [r for r in rs if "Ca" in r["elements"] and "H" in r["elements"]]
print("\nCa + H фаз:", len(ca))
for r in ca[:20]:
    print("  {:<22} {:<14} {:<24} d1={:.3f}".format(r["formula"], r["id"], r["name"][:24], top(r)[0]))

big = [r for r in rs if top(r)[0] > 7.0]
print("\nфаз с сильнейшей линией больше 7 A:", len(big))
for r in big[:15]:
    print("  {:<22} {:<14} {:<24} d1={:.3f}".format(r["formula"], r["id"], r["name"][:24], top(r)[0]))

want = ("CaH4SO6", "Ca(HO)2", "CaH2O2", "Al(HO)3", "Mg(HO)2", "Ca3Al2(HO)12")
print("\nцелевые фазы:")
for r in rs:
    f = r["formula"].replace("(", "").replace(")", "")
    for w in want:
        if f == w.replace("(", "").replace(")", ""):
            print("  {:<22} {:<14} {:<28} d1={:.3f}".format(r["formula"], r["id"], r["name"][:28], top(r)[0]))
