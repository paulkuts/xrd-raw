#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Сборка компактной базы эталонных фаз для портативного XRD-конвертера.

Вход:  cod_phases.jsonl    (COD, CC0) — экспериментальные структуры, приоритет 1
       jarvis_phases.jsonl (JARVIS-DFT, NIST, CC BY 4.0) — расчётные, приоритет 2
Выход: phases_db.json.gz — {"meta": {...}, "phases": [...]}
       Пики в базе: [2θ(Cu Kα), d (Å), I отн. (%)]

После сборки печатается проверка по контрольным линиям известных фаз: если
эталон сдвинут относительно эксперимента, видно сразу.
"""
import gzip
import json
import math
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT = BASE / "phases_db.json.gz"
WAVELENGTH = 1.5418      # Cu Kα средневзвешенная, как в peakpick.py

# контрольные фазы: (формула без скобок, опорные линии d Å, подпись)
CONTROLS = [
    ("SiO2", [3.343, 4.257, 1.817, 2.457], "кварц alpha"),
    ("SiO2", [4.040, 2.485, 2.840], "кристобалит"),
    ("CaCO3", [3.035, 2.285, 2.095, 1.913, 1.875], "кальцит"),
    ("CaCO3", [3.396, 3.273, 2.700, 2.341], "арагонит"),
    ("CaHO2", [4.900, 2.628, 1.927, 1.796], "портландит"),
    ("MgHO2", [4.770, 2.365, 1.794], "брусит"),
    ("CaH4SO6", [7.630, 4.280, 3.060, 2.870], "гипс"),
    ("CaSO4", [3.492, 2.849, 2.328], "ангидрит"),
    ("Al2O3", [2.552, 2.085, 1.601], "корунд"),
    ("TiO2", [3.248, 1.687, 2.488], "рутил"),
    ("Fe2O3", [2.700, 2.519, 1.694], "гематит"),
    ("Fe3O4", [2.530, 1.485, 1.616], "магнетит"),
    ("MgO", [2.106, 1.489], "периклаз"),
    ("CaO", [2.405, 1.700], "известь"),
    ("Cu", [2.088, 1.808, 1.278], "медь"),
    ("Al", [2.338, 2.024, 1.431], "алюминий"),
    ("Si", [3.135, 1.920, 1.637], "кремний"),
    ("NaCl", [2.821, 1.994, 1.628], "галит"),
    ("Ca3SiO5", [2.780, 2.744, 2.610, 1.760], "алит C3S"),
    ("Ca2SiO4", [2.791, 2.745, 2.606, 2.875], "белит C2S"),
    ("Ca3Al2O6", [2.696, 1.905, 2.686], "C3A"),
    ("Ca2Fe2O5", [2.644, 1.918, 2.623], "браунмиллерит C4AF"),
    ("Ca6Al2S3O49|ettringite", [9.730, 5.610, 3.880, 2.564, 2.209], "эттрингит"),
]


def tt_of_d(d):
    x = WAVELENGTH / (2.0 * d)
    if x >= 1.0 or x <= 0:
        return None
    return 2.0 * math.degrees(math.asin(x))


def tidy_name(name, formula):
    name = re.sub(r"\s+", " ", (name or "")).strip()
    if name and re.search(r"[A-Za-zА-Яа-я]{3}", name) and name not in (";", "-"):
        return name[:60]
    return formula


def load(name, prio):
    p = BASE / name
    if not p.exists():
        print("нет файла {}".format(name))
        return []
    out = []
    with open(p, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            peaks = []
            for d, inten in r.get("peaks", []):
                tt = tt_of_d(d)
                if tt is None or tt < 5 or tt > 100:
                    continue
                peaks.append([round(tt, 3), round(d, 4), round(inten, 1)])
            if len(peaks) < 2:
                continue
            peaks.sort(key=lambda x: x[0])
            r["peaks"] = peaks[:24]
            r["prio"] = prio
            r["name"] = tidy_name(r.get("name"), r["formula"])
            r.pop("term", None)
            if not r.get("id") and r.get("src") == "JARVIS-DFT":
                r["id"] = r.get("id") or ""
            out.append(r)
    print("{}: {}".format(name, len(out)))
    return out


def dedupe(phases):
    """Дедуп: та же формула + близкий набор сильных линий; COD важнее JARVIS."""
    kept = []
    by_formula = {}
    for p in phases:
        by_formula.setdefault(p["formula"], []).append(p)
    for formula, group in by_formula.items():
        group.sort(key=lambda p: p["prio"])
        taken = []
        for p in group:
            sig = [round(x[1], 2) for x in sorted(p["peaks"], key=lambda y: -y[2])[:6]]
            dup = False
            for q in taken:
                qsig = [round(x[1], 2) for x in sorted(q["peaks"], key=lambda y: -y[2])[:6]]
                same = sum(1 for a in sig if any(abs(a - b) <= 0.05 for b in qsig))
                if same >= max(4, len(sig) - 1):
                    dup = True
                    break
            if not dup:
                taken.append(p)
                kept.append(p)
    return kept


def norm_formula(f):
    """Формула без скобок, чтобы Ca(HO)2 и CaH2O2 считались одной фазой."""
    return (f or "").replace("(", "").replace(")", "")


def check_controls(phases):
    print("\nконтроль эталонов (доля найденных опорных линий):")
    ok = part = bad = miss = 0
    for formula, refs, label in CONTROLS:
        cands = [p for p in phases
                 if norm_formula(p["formula"]) == formula.split("|")[0]
                 or (len(formula.split("|")) > 1
                     and formula.split("|")[1].lower() in (p.get("name") or "").lower())]
        if not cands:
            print("  --  {:22s} {:18s} нет в базе".format(label, formula))
            miss += 1
            continue
        best_ph, best_hit = None, -1
        for p in cands:
            dd = [x[1] for x in p["peaks"]]
            hit = sum(1 for r in refs if any(abs(r - d) <= 0.05 for d in dd))
            if hit > best_hit:
                best_ph, best_hit = p, hit
        share = best_hit / len(refs)
        if share == 1.0:
            ok += 1
            mark = "OK "
        elif share >= 0.6:
            part += 1
            mark = "~  "
        else:
            bad += 1
            mark = "!! "
        print("  {} {:22s} {:18s} линий {}/{} ({}), {}".format(
            mark, label, formula[:18], best_hit, len(refs), best_ph.get("id", ""),
            best_ph.get("src", "")))
    print("  итог: полностью {}, частично {}, плохо {}, нет в базе {}".format(
        ok, part, bad, miss))


def main():
    cod = load("cod_phases.jsonl", 1)
    cod += load("cod_extra.jsonl", 1)
    jar = load("jarvis_phases.jsonl", 2)
    phases = dedupe(cod + jar)
    print("после дедупликации: {}".format(len(phases)))

    n_cod = sum(1 for p in phases if p["prio"] == 1)
    els = set()
    for p in phases:
        els.update(p["elements"])
    meta = {
        "version": "2.0",
        "wavelength": "CuKa 1.5418 A",
        "count": len(phases),
        "from_COD": n_cod,
        "from_JARVIS": len(phases) - n_cod,
        "sources": ("COD (crystallography.net, CC0) - experimental CIF; "
                    "JARVIS-DFT (NIST, CC BY 4.0) - DFT structures with ICSD reference"),
        "elements": sorted(els),
        "note": ("Peaks calculated from crystal structures (pymatgen XRDCalculator). "
                 "Not the PDF ICDD database; legally valid identification against the "
                 "PDF card file requires an ICDD license."),
    }
    db = {"meta": meta, "phases": phases}
    with gzip.open(str(OUT), "wt", encoding="utf-8", compresslevel=9) as f:
        json.dump(db, f, ensure_ascii=False, separators=(",", ":"))
    print("записано: {} ({:.1f} МБ)".format(OUT, OUT.stat().st_size / 1024 / 1024))
    print("фаз: {} (COD {}, JARVIS {})".format(len(phases), n_cod, len(phases) - n_cod))
    print("элементов в базе: {}".format(len(els)))
    check_controls(phases)


if __name__ == "__main__":
    main()
