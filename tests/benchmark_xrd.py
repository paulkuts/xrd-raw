#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Бенчмарк конвейера xrd-raw на синтетических фикстурах (исполнитель D).

Считает измеримые метрики, которые не зависят от «красоты» отчёта:

  Фазы
    true_phases      — сколько фаз реально есть в образце (из truth-JSON);
    found_true       — сколько из них найдено в результате explain();
    recall           — found_true / true_phases;
    candidates       — сколько фаз вернул конвейер (выше min_fom);
    false_candidates — из них не совпавших ни с одной истинной фазой;
    precision        — (candidates - false_candidates) / candidates.

  Пики
    peaks            — сколько пиков нашёл peakpick.pick();
    true_lines       — сколько истинных линий попало в диапазон скана (I >= 3 %);
    lines_detected   — сколько из них имеют пик в пределах 0.20°;
    spurious_peaks   — пики дальше 0.20° от любой истинной линии;
    unexplained      — пики, которые explain() не отнёс ни к одной фазе.

  Время
    t_pick, t_explain — секунды (полная база 32 491 фаза).

estimate_shift НЕ запускается: на полной базе это ~65 с на один вызов
(см. docs/VALIDATION.md, раздел «Ограничения методики»).

Запуск:
    /root/work/.venv/bin/python tests/benchmark_xrd.py --out tmp/bench.json
    /root/work/.venv/bin/python tests/benchmark_xrd.py --only mix_qcg

Вывод: таблица в консоль + машиночитаемый JSON (--out).
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TOOLS = ROOT / "tools"
FIX = HERE / "fixtures"
sys.path.insert(0, str(TOOLS))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp"))

for _stream in (sys.stdout, sys.stderr):     # кириллица в cp866/cp1251 консоли
    try:
        _stream.reconfigure(errors="replace")
    except Exception:                        # noqa: BLE001
        pass

import peakpick  # noqa: E402
import phases  # noqa: E402
import winxrd2xy  # noqa: E402

TOL = 0.25            # допуск сопоставления фаз (как phases.TOL_DEFAULT)
PEAK_MATCH = 0.20     # допуск «пик соответствует истинной линии», градусы
MIN_FOM = 55.0
MIN_LINE_I = 3.0      # линии слабее 3 % в синтетике не видны

FIXTURES = [
    "quartz_pure", "mix_qcg", "overlap_calcite_gypsum", "noise_slope",
    "shift_plus010", "narrow_20_35", "wrongwl_co_as_cu", "noise_only",
    "dft_weak",
]


def formula_elements(formula):
    import re

    return set(re.findall(r"[A-Z][a-z]?", formula or ""))


def same_phase(cand, truth_phase, tol=TOL, min_frac=0.6, min_int=10.0):
    """Та же логика, что в tests/test_regression.py (набор линий + химия)."""
    el_t = formula_elements(truth_phase.get("formula", ""))
    el_c = set(cand.get("elements") or [])
    if el_t and el_c and el_t != el_c:
        return False
    lines = truth_phase.get("lines_in_range")
    if lines is None:
        pts = truth_phase.get("all_lines_d_I") or []
        lines = [{"two_theta": a, "I_percent": b} for a, b in pts]
    pos = [l["two_theta"] for l in lines if l["I_percent"] >= min_int] or \
          [l["two_theta"] for l in lines]
    if not pos:
        return False
    refs = [p[0] for p in cand.get("peaks", [])]
    hits = sum(1 for tt in pos if any(abs(tt - r) <= tol for r in refs))
    return hits / float(len(pos)) >= min_frac


def truth_line_positions(truth_phase, lam=1.54184):
    """Позиции 2θ истинных линий (d -> 2θ, если в truth записаны только d)."""
    lines = truth_phase.get("lines_in_range")
    if lines is not None:
        return [l["two_theta"] for l in lines if l["I_percent"] >= MIN_LINE_I]
    out = []
    for d, rel in truth_phase.get("all_lines_d_I") or []:
        if rel < MIN_LINE_I:
            continue
        x = lam / (2.0 * d)
        if x < 1.0:
            out.append(2.0 * math.degrees(math.asin(x)))
    return out


def run_fixture(name):
    doc = json.loads((FIX / (name + ".json")).read_text(encoding="utf-8"))
    truth = doc["truth"]
    t_parse = time.time()
    data, step, meta = winxrd2xy.parse(str(FIX / (name + ".raw")))
    t_parse = time.time() - t_parse
    x = [a for a, _ in data]
    y = [i for _, i in data]

    t0 = time.time()
    peaks = peakpick.pick(x, y)
    t_pick = time.time() - t0

    db = phases.load_db()
    t0 = time.time()
    res = phases.explain(peaks, db, tol=TOL, min_fom=MIN_FOM)
    t_explain = time.time() - t0

    # --- фазы
    present = truth["present"]
    matched = []
    for tp in present:
        for sc in res["phases"]:
            if same_phase(sc["phase"], tp):
                matched.append((tp.get("key") or tp.get("name"),
                                round(sc["fom"], 1)))
                break
    found_keys = {k for k, _ in matched}
    n_true = len(present)
    n_found = len(found_keys)
    cand = res["phases"]
    false_c = []
    for sc in cand:
        if not any(same_phase(sc["phase"], tp) for tp in present):
            false_c.append({"name": (sc["phase"].get("name") or "")[:40],
                            "formula": sc["phase"].get("formula", ""),
                            "fom": round(sc["fom"], 1),
                            "src": sc["phase"].get("src", "")})

    # --- пики
    lines = []
    for tp in present:
        lines.extend(truth_line_positions(tp, doc["wavelength_A"]))
    lines.sort()
    hit_lines, spurious = 0, []
    for tt in lines:
        if any(abs(p["two_theta"] - tt) <= PEAK_MATCH for p in peaks):
            hit_lines += 1
    for p in peaks:
        if not any(abs(p["two_theta"] - tt) <= PEAK_MATCH for tt in lines):
            spurious.append(round(p["two_theta"], 2))

    per_peak = res.get("per_peak") or {}
    unexplained = [round(peaks[i]["two_theta"], 2) for i, p in enumerate(peaks)
                   if p["intensity"] >= phases.MIN_OBS_INT and i not in per_peak]

    row = {
        "fixture": name,
        "range": [doc["two_theta"]["start"], doc["two_theta"]["stop"]],
        "points": len(data),
        "phases_true": [tp.get("key") or tp.get("name") for tp in present],
        "phases_found": sorted(found_keys),
        "phases_found_rank_fom": {k: v for k, v in matched},
        "true_phases": n_true,
        "found_true": n_found,
        "recall": (n_found / float(n_true)) if n_true else None,
        "candidates": len(cand),
        "false_candidates": len(false_c),
        "false_candidate_list": false_c,
        "precision": ((len(cand) - len(false_c)) / float(len(cand))
                      if cand else None),
        "peaks": len(peaks),
        "n_obs": res.get("n_obs"),
        "true_lines": len(lines),
        "lines_detected": hit_lines,
        "line_recall": (hit_lines / float(len(lines))) if lines else None,
        "spurious_peaks": len(spurious),
        "spurious_positions": spurious[:15],
        "unexplained": len(unexplained),
        "t_parse_s": round(t_parse, 3),
        "t_pick_s": round(t_pick, 3),
        "t_explain_s": round(t_explain, 2),
        "expect_no_confident_candidate": bool(
            truth.get("expect_no_confident_candidate")),
    }
    return row


def aggregate(rows):
    n_true = sum(r["true_phases"] for r in rows)
    n_found = sum(r["found_true"] for r in rows)
    n_cand = sum(r["candidates"] for r in rows)
    n_false = sum(r["false_candidates"] for r in rows)
    lines = sum(r["true_lines"] for r in rows)
    hit = sum(r["lines_detected"] for r in rows)
    return {
        "fixtures": len(rows),
        "true_phases": n_true,
        "found_true": n_found,
        "recall": round(n_found / float(n_true), 3) if n_true else None,
        "candidates": n_cand,
        "false_candidates": n_false,
        "precision": round((n_cand - n_false) / float(n_cand), 3) if n_cand else None,
        "true_lines": lines,
        "lines_detected": hit,
        "line_recall": round(hit / float(lines), 3) if lines else None,
        "spurious_peaks": sum(r["spurious_peaks"] for r in rows),
        "unexplained": sum(r["unexplained"] for r in rows),
        "runtime_s": round(sum(r["t_parse_s"] + r["t_pick_s"] + r["t_explain_s"]
                               for r in rows), 2),
    }


def print_table(rows, tot):
    head = ("{:<22s} {:>5s} {:>6s} {:>7s} {:>7s} {:>8s} {:>6s} {:>6s} {:>8s} {:>7s}")
    print(head.format("фикстура", "точек", "пиков", "линии", "линии", "пики",
                      "фазы", "ложные", "необъяс", "explain"))
    print(head.format("", "", "", "всего", "найд.", "вне линий", "найд.",
                      "кандид.", "пики", "с"))
    print("-" * 92)
    for r in rows:
        print("{:<22s} {:>5d} {:>6d} {:>7d} {:>7d} {:>8d} {:>3d}/{:<2d} "
              "{:>6d} {:>8d} {:>7.2f}".format(
                  r["fixture"], r["points"], r["peaks"], r["true_lines"],
                  r["lines_detected"], r["spurious_peaks"],
                  r["found_true"], max(1, r["true_phases"]),
                  r["false_candidates"], r["unexplained"], r["t_explain_s"]))
    print("-" * 92)
    print("итого: фаз найдено {}/{} (recall {}), кандидатов {}, из них ложных {} "
          "(precision {}), линий {}/{} ({}), шумовых пиков {}, время {:.1f} с".format(
              tot["found_true"], tot["true_phases"], tot["recall"],
              tot["candidates"], tot["false_candidates"], tot["precision"],
              tot["lines_detected"], tot["true_lines"], tot["line_recall"],
              tot["spurious_peaks"], tot["runtime_s"]))


def self_check(rows):
    """Предупреждения о метриках, которые нельзя трактовать как «качество»."""
    notes = []
    for r in rows:
        if r["fixture"] == "wrongwl_co_as_cu" and r["found_true"]:
            notes.append("wrongwl: при неверной длине волны «найдены» истинные "
                         "фазы — вероятно, случайное совпадение; см. VALIDATION.md")
        if r["true_phases"] == 0 and r["candidates"]:
            notes.append("{}: {} кандидатов на образце без фаз (ложное "
                         "определение)".format(r["fixture"], r["candidates"]))
        if r["peaks"] >= 120:
            notes.append("{}: {} пиков — упор в max_peaks (шум)".format(
                r["fixture"], r["peaks"]))
    return notes


def main(argv=None):
    ap = argparse.ArgumentParser(description="Бенчмарк xrd-raw на фикстурах")
    ap.add_argument("--out", default=str(ROOT / "tmp" / "bench.json"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    names = [a.only] if a.only else FIXTURES
    rows = []
    for name in names:
        row = run_fixture(name)
        rows.append(row)
        if not a.quiet:
            print("  {:<22s} explain {:5.2f} c, пиков {}, фаз {}".format(
                name, row["t_explain_s"], row["peaks"], row["candidates"]))

    tot = aggregate(rows)
    db = phases.load_db()
    meta, n_db = phases.db_info(db)
    out = {
        "generated_by": "tests/benchmark_xrd.py",
        "synthetic_only": True,
        "warnings": self_check(rows),
        "code_state_note": (
            "Снимок текущего кода ветки feat/agent-D. Ветки A (phases), "
            "B (parser/peaks), C (export) правятся параллельно: числа "
            "«ДО» относятся к состоянию до их интеграции. Финальное «ПОСЛЕ» "
            "заполняет ведущий после слияния."),
        "database": {"phases": n_db,
                     "sources": meta.get("sources", ""),
                     "wavelength": meta.get("wavelength", "")},
        "settings": {"tol_deg": TOL, "peak_match_deg": PEAK_MATCH,
                     "min_fom": MIN_FOM, "min_line_intensity": MIN_LINE_I},
        "fixtures": rows,
        "totals": tot,
    }
    p = Path(a.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n",
                 encoding="utf-8")
    if not a.quiet:
        print()
        print_table(rows, tot)
        for w in out["warnings"]:
            print("  ! " + w)
        print("\nJSON: {}".format(p))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
