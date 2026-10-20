#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Сборка эталонных фаз из дампа CIF Crystallography Open Database — многопроцессно.

Источник: sreekanth22/COD_CIF (HuggingFace; сами данные COD — CC0).
Архив читается потоком в главном процессе, разбор CIF и расчёт паттерна
(Cu Kα) идут в пуле процессов. На диск пишется только результат.

  python build_cod_dump_mp.py --workers 10
"""
import argparse
import json
import re
import sys
import tarfile
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

warnings.filterwarnings("ignore")

BASE = Path(__file__).resolve().parent
ARCHIVE = BASE / "cod_cifs.tar.gz"
OUT = BASE / "cod_phases.jsonl"

RE_FORM = re.compile(r"_chemical_formula_sum\s+'?([^'\n]+)")
RE_MINERAL = re.compile(r"_chemical_name_mineral\s+'?([^'\n]+)")
RE_SYSTEM = re.compile(r"_chemical_name_systematic\s+'?([^'\n]+)")
RE_COMMON = re.compile(r"_chemical_name_common\s+'?([^'\n]+)")
RE_SG = re.compile(r"_symmetry_space_group_name_H-M\s+'?([^'\n]+)")

_XRD = None


def _init_worker():
    global _XRD
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    _XRD = XRDCalculator(wavelength="CuKa")


def carbon_count(formula):
    m = re.search(r"C(\d*)(?![a-z])", formula or "")
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


def looks_organic(formula):
    if not formula:
        return False
    s = set(re.findall(r"[A-Z][a-z]?", formula.replace(" ", "")))
    if "C" not in s or "H" not in s:
        return False
    return carbon_count(formula) > 2 or len(s) > 6


def clean_name(text):
    for rx in (RE_MINERAL, RE_SYSTEM, RE_COMMON):
        m = rx.search(text)
        if m:
            name = re.sub(r"\s+", " ", m.group(1).strip().strip('"').strip("'")).strip()
            if len(name) >= 3 and re.search(r"[A-Za-z]{3}", name) and name not in (";", "-"):
                return name[:80]
    return ""


def process(item):
    """Один CIF -> запись фазы или код пропуска."""
    cid, text = item
    m = RE_FORM.search(text)
    if m and looks_organic(m.group(1).strip()):
        return None, "organic"
    from pymatgen.core import Structure
    try:
        s = Structure.from_str(text, fmt="cif")
    except Exception:  # noqa: BLE001
        return None, "bad"
    if len(s) > 300:
        return None, "big"
    comp = s.composition
    els = sorted(e.symbol for e in comp.elements)
    formula = comp.reduced_formula
    if "C" in els and "H" in els and carbon_count(formula) > 2:
        return None, "organic"
    try:
        pat = _XRD.get_pattern(s, two_theta_range=(5.0, 90.0), scaled=True)
    except Exception:  # noqa: BLE001
        return None, "bad"
    peaks = []
    for inten, d in zip(pat.y, pat.d_hkls):
        if inten < 1.0:
            continue
        peaks.append([round(float(d), 4), round(float(inten), 1)])
        if len(peaks) >= 24:
            break
    if len(peaks) < 2:
        return None, "few"
    m = RE_SG.search(text)
    sg = re.sub(r"\s+", " ", m.group(1)).strip()[:24] if m else ""
    return {
        "src": "COD", "id": "COD {}".format(cid), "name": clean_name(text),
        "formula": formula, "elements": els, "sg": sg, "term": "", "peaks": peaks,
    }, "ok"


def iter_cifs(limit_files=0):
    with tarfile.open(str(ARCHIVE), "r:gz") as tf:
        n = 0
        for member in tf:
            if not member.isfile() or not member.name.lower().endswith(".cif"):
                continue
            n += 1
            if limit_files and n > limit_files:
                break
            try:
                raw = tf.extractfile(member).read()
            except Exception:  # noqa: BLE001
                continue
            yield Path(member.name).stem, raw.decode("utf-8", "replace")


def main(workers=10, max_phases=0, limit_files=0):
    if not ARCHIVE.exists():
        sys.exit("нет архива {}".format(ARCHIVE))
    seen = set()
    stats = {}
    n_in = n_ok = 0
    t0 = time.time()
    with open(OUT, "w", encoding="utf-8") as out, Pool(workers, initializer=_init_worker) as pool:
        for rec, status in pool.imap(process, iter_cifs(limit_files), chunksize=40):
            n_in += 1
            stats[status] = stats.get(status, 0) + 1
            if rec is None:
                continue
            key = (rec["formula"], tuple(round(p[0], 2) for p in rec["peaks"][:5]))
            if key in seen:
                stats["dup"] = stats.get("dup", 0) + 1
                continue
            seen.add(key)
            out.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n_ok += 1
            if n_ok % 2000 == 0:
                el = time.time() - t0
                print("файлов {}: фаз {}, {:.0f} фаз/мин, пропуски {}".format(
                    n_in, n_ok, n_ok / max(el, 1) * 60, stats), flush=True)
            if max_phases and n_ok >= max_phases:
                print("лимит фаз {}".format(max_phases), flush=True)
                pool.terminate()
                break
    print("готово: файлов {}, фаз {}, время {:.1f} мин, {}".format(
        n_in, n_ok, (time.time() - t0) / 60, stats), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--max-phases", type=int, default=0)
    ap.add_argument("--limit-files", type=int, default=0)
    a = ap.parse_args()
    main(a.workers, a.max_phases, a.limit_files)
