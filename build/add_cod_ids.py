# -*- coding: utf-8 -*-
"""Добор отдельных фаз из COD по номерам записей (CIF) -> cod_extra.jsonl.

Нужен для фаз, которых нет в готовом дампе CIF (например, гипс, цементные
клинкерные фазы). Качает по одному CIF, считает паттерн (Cu Kα).

  python add_cod_ids.py --ids "1010981 1011074" --out cod_extra.jsonl
  python add_cod_ids.py --file ids.txt
"""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

URL = "https://www.crystallography.net/cod/{}.cif"
UA = {"User-Agent": "XRD-phase-build/1.0 (research use; COD is CC0)"}
RE_FORM = re.compile(r"_chemical_formula_sum\s+'?([^'\n]+)")
RE_MINERAL = re.compile(r"_chemical_name_mineral\s+'?([^'\n]+)")
RE_SYSTEM = re.compile(r"_chemical_name_systematic\s+'?([^'\n]+)")
RE_COMMON = re.compile(r"_chemical_name_common\s+'?([^'\n]+)")
RE_SG = re.compile(r"_symmetry_space_group_name_H-M\s+'?([^'\n]+)")


def carbon_count(formula):
    m = re.search(r"C(\d*)(?![a-z])", formula or "")
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


def clean_name(text):
    for rx in (RE_MINERAL, RE_SYSTEM, RE_COMMON):
        m = rx.search(text)
        if m:
            name = re.sub(r"\s+", " ", m.group(1).strip().strip('"').strip("'")).strip()
            if len(name) >= 3 and re.search(r"[A-Za-z]{3}", name) and name not in (";", "-"):
                return name[:80]
    return ""


def fetch(cid, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL.format(cid), headers=UA)
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read().decode("utf-8", "replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            print("  {} попытка {}: {}".format(cid, attempt + 1, e), flush=True)
            time.sleep(3 + attempt * 3)
    return None


def main(ids, out_path, delay=1.5):
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    xrd = XRDCalculator(wavelength="CuKa")
    out = open(out_path, "a", encoding="utf-8", buffering=1)
    ok = skipped = failed = 0
    for cid in ids:
        cid = str(cid).strip()
        if not cid:
            continue
        text = fetch(cid)
        time.sleep(delay)
        if not text or "_cell_length_a" not in text:
            print("{}: нет данных".format(cid), flush=True)
            failed += 1
            continue
        m = RE_FORM.search(text)
        if m and "C" in m.group(1) and "H" in m.group(1) and carbon_count(m.group(1)) > 2:
            print("{}: органика, пропуск".format(cid), flush=True)
            skipped += 1
            continue
        try:
            s = Structure.from_str(text, fmt="cif")
            pat = xrd.get_pattern(s, two_theta_range=(5.0, 90.0), scaled=True)
        except Exception as e:  # noqa: BLE001
            print("{}: разбор не удался ({})".format(cid, e), flush=True)
            failed += 1
            continue
        peaks = []
        for inten, d in zip(pat.y, pat.d_hkls):
            if inten < 1.0:
                continue
            peaks.append([round(float(d), 4), round(float(inten), 1)])
            if len(peaks) >= 24:
                break
        if len(peaks) < 2:
            failed += 1
            continue
        sm = RE_SG.search(text)
        rec = {
            "src": "COD", "id": "COD {}".format(cid), "name": clean_name(text),
            "formula": s.composition.reduced_formula,
            "elements": sorted(e.symbol for e in s.composition.elements),
            "sg": re.sub(r"\s+", " ", sm.group(1)).strip()[:24] if sm else "",
            "term": "", "peaks": peaks,
        }
        out.write(json.dumps(rec, ensure_ascii=False) + "\n")
        top = max(peaks, key=lambda x: x[1])
        print("{}: {} {} | пик {} d={:.3f}".format(
            cid, rec["formula"], rec["name"][:28], len(peaks), top[0]), flush=True)
        ok += 1
    out.close()
    print("\nдобавлено {}, пропущено {}, ошибок {}".format(ok, skipped, failed))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="")
    ap.add_argument("--file", default="")
    ap.add_argument("--out", default="cod_extra.jsonl")
    ap.add_argument("--delay", type=float, default=1.5)
    a = ap.parse_args()
    ids = a.ids.split() if a.ids else []
    if a.file and Path(a.file).exists():
        ids += Path(a.file).read_text(encoding="utf-8").split()
    if not ids:
        sys.exit("нет ID")
    main(ids, a.out, a.delay)
