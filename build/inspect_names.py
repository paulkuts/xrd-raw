# -*- coding: utf-8 -*-
"""Поиск нужных фаз в собранном cod_phases.jsonl по названию/формуле."""
import json
from pathlib import Path

KEYS = [
    "gypsum", "bassanite", "anhydrite", "ettringite", "thaumasite",
    "portlandite", "hatrurite", "larnite", "alite", "belite", "brownmillerite",
    "mayenite", "tricalcium", "dicalcium", "hydrocalumite", "katoite",
    "calcite", "aragonite", "vaterite", "dolomite", "magnesite", "brucite",
    "goethite", "hematite", "magnetite", "quartz", "cristobalite", "tridymite",
    "gibbsite", "boehmite", "diaspore", "anatase", "rutile", "halite",
    "syngenite", "thenardite", "mirabilite", "arcanite", "aphthitalite",
    "wollastonite", "xonotlite", "tobermorite", "afwillite", "hydrogarnet",
    "ettringite", "monosulfate", "friedel", "periclase", "lime", "portland",
]


def top(r):
    return max(r["peaks"], key=lambda x: x[1])


def main():
    rs = [json.loads(l) for l in open("cod_phases.jsonl", encoding="utf-8")]
    print("всего фаз в дампе:", len(rs))
    for key in KEYS:
        hits = [r for r in rs if key in (r.get("name") or "").lower()]
        if not hits:
            continue
        print("\n[{}] найдено {}".format(key, len(hits)))
        for r in hits[:6]:
            print("   {:<20} {:<14} {:<28} d1={:6.3f}".format(
                r["formula"][:20], r["id"], (r.get("name") or "")[:28], top(r)[0]))


if __name__ == "__main__":
    main()
