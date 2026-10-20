#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Идентификация кристаллических фаз по пикам дифрактограммы.

Работает офлайн: база эталонных порошковых паттернов лежит рядом
(phases_db.json.gz). Логика — как в поиске по картотеке PDF: сравнение
d-позиций и относительных интенсивностей, отбор по сильнейшим линиям
(Hanawalt-подобный), затем жадное объяснение всего набора пиков.
"""

import gzip
import json
import math
from pathlib import Path

TOL_DEFAULT = 0.25          # допуск по 2θ, градусы
MIN_OBS_INT = 3.0           # пики слабее (% от макс.) не участвуют в поиске фаз
MAX_PHASES = 8              # сколько фаз максимум объясняют образец
DFT_PENALTY = 8.0           # штраф фазам из DFT-базы (нет экспериментального эталона)
CACHE = {}


# ---------------------------------------------------------------- загрузка базы
def load_db(path=None):
    """Загрузить базу фаз (кэшируется в памяти процесса)."""
    if path is None:
        path = Path(__file__).resolve().parent / "phases_db.json.gz"
    else:
        path = Path(path)
    key = str(path)
    if key in CACHE:
        return CACHE[key]
    if not path.exists():
        raise FileNotFoundError("Нет файла базы фаз: {}".format(path))
    with gzip.open(str(path), "rt", encoding="utf-8") as f:
        db = json.load(f)
    phases = db.get("phases", [])
    for p in phases:
        p.setdefault("elements", [])
        p.setdefault("peaks", [])
    db["phases"] = phases
    CACHE[key] = db
    return db


def db_info(db):
    return db.get("meta", {}), len(db.get("phases", []))


# ------------------------------------------------------------------- поиск фаз
def _index_by_angle(phases, step=0.1):
    """Индекс: угловой бакет -> список (индекс фазы, индекс рефлекса)."""
    idx = {}
    for pi, ph in enumerate(phases):
        for ri, pk in enumerate(ph["peaks"]):
            b = int(pk[0] // step)
            idx.setdefault(b, []).append((pi, ri))
    return idx, step


def _candidate_phases(obs_peaks, index, step, tol):
    """Фазы, у которых есть хоть один рефлекс рядом с наблюдаемым пиком."""
    cand = {}
    span = int(math.ceil(tol / step)) + 1
    for oi, op in enumerate(obs_peaks):
        b0 = int(op["two_theta"] // step)
        for b in range(b0 - span, b0 + span + 1):
            for pi, ri in index.get(b, ()):
                cand.setdefault(pi, set()).add(oi)
    return cand


def _score_phase(ph, obs_peaks, tol, allowed_obs=None, unmatched_penalty=True,
                 max_ref_lines=12, shift=0.0):
    """Сопоставить рефлексы фазы с наблюдаемыми пиками.

    В расчёт берём не все рефлексы, а сильнейшие (по умолчанию 12): слабые
    линии в реальном спектре тонут в фоне, и требовать их совпадения нельзя —
    иначе сильные фазы проигрывают экзотике со множеством слабых отражений.
    shift — систематический сдвиг эталона по 2θ (нулевая точка съёмки).
    Возвращает dict с метриками или None, если фаза не проходит отбор
    по сильнейшим линиям.
    """
    all_refs = ph["peaks"]
    if not all_refs:
        return None
    order = sorted(range(len(all_refs)), key=lambda i: -all_refs[i][2])[:max_ref_lines]
    order.sort()
    refs = [all_refs[i] for i in order]
    if allowed_obs is None:
        allowed = range(len(obs_peaks))
    else:
        allowed = allowed_obs

    used_obs = set()
    matched = []                       # (ref_index, obs_index, dtheta)
    for ri, (rtt, rd, rint) in enumerate(refs):
        best = None
        for oi in allowed:
            if oi in used_obs:
                continue
            dth = abs(obs_peaks[oi]["two_theta"] - (rtt + shift))
            if dth <= tol and (best is None or dth < best[0]):
                best = (dth, oi)
        if best is not None:
            used_obs.add(best[1])
            matched.append((ri, best[1], best[0]))

    if not matched:
        return None

    sum_ref = sum(max(r[2], 0.5) for r in refs)
    sum_obs = sum(max(obs_peaks[oi]["intensity"], 0.5) for oi in allowed) or 1.0
    cov_ref = sum(max(refs[ri][2], 0.5) for ri, _, _ in matched) / sum_ref
    cov_obs = sum(max(obs_peaks[oi]["intensity"], 0.5) for _, oi, _ in matched) / sum_obs

    # сильнейшие линии эталона (Hanawalt): три самых интенсивных рефлекса
    key_idx = sorted(range(len(refs)), key=lambda i: -refs[i][2])[:3]
    matched_idx = {ri for ri, _, _ in matched}
    key_hits = sum(1 for i in key_idx if i in matched_idx)
    if key_hits == 0:
        return None
    # жёсткий отбор: экзотика иначе цепляется за пару случайных пиков
    if len(matched) < 2:
        return None
    if len(refs) >= 6 and len(matched) < 3:
        return None
    if key_hits < 2 and cov_ref < 0.45:
        return None
    if key_hits < 3 and cov_ref < 0.40:
        return None

    fom = 100.0 * (0.50 * cov_ref + 0.40 * cov_obs + 0.10 * (key_hits / float(len(key_idx))))

    if unmatched_penalty:
        strong_missed = [refs[i] for i in key_idx if i not in matched_idx and refs[i][2] >= 40.0]
        fom -= 10.0 * len(strong_missed)

    mean_dtheta = sum(d for _, _, d in matched) / len(matched)
    pen = 0.0 if ph.get("prio", 2) == 1 else DFT_PENALTY
    return {
        "fom": max(0.0, fom - pen),
        "fom_raw": fom,
        "cov_ref": cov_ref,
        "cov_obs": cov_obs,
        "key_hits": key_hits,
        "key_total": len(key_idx),
        "n_matched": len(matched),
        "n_ref": len(refs),
        "mean_dtheta": mean_dtheta,
        "matched": matched,
    }


def estimate_shift(peaks, db, tol=TOL_DEFAULT, span=0.40, step=0.05, top_n=3):
    """Подобрать систематический сдвиг нуля по 2θ (плоскость образца, юстировка).

    Перебираем сдвиги эталонов и берём тот, при котором лучшие фазы объясняют
    образец заметно лучше всего. Возвращает (shift, gain).
    """
    obs = [p for p in peaks if p["intensity"] >= MIN_OBS_INT] or peaks
    phases = db["phases"]
    index, idx_step = _index_by_angle(phases)
    cand = list(_candidate_phases(obs, index, idx_step, tol + span).keys())
    if not cand:
        return 0.0, 0.0

    def score(sh):
        best = []
        for pi in cand:
            sc = _score_phase(phases[pi], obs, tol, shift=sh)
            if sc is not None:
                best.append(sc["fom"])
        best.sort(reverse=True)
        return sum(best[:top_n])

    base = score(0.0)
    best_sh, best_val = 0.0, base
    n = int(round(span / step))
    for k in range(-n, n + 1):
        if k == 0:
            continue
        sh = round(k * step, 3)
        v = score(sh)
        if v > best_val + 1e-9:
            best_val, best_sh = v, sh
    # принимаем сдвиг, только если он даёт заметный выигрыш (>8%)
    if base > 0 and best_val < base * 1.08:
        return 0.0, 0.0
    return best_sh, best_val - base


def search_phases(peaks, db, tol=TOL_DEFAULT, min_fom=45.0, top=30,
                  elements_only=None, shift=0.0):
    """Ранжированный список фаз-кандидатов по всему набору пиков."""
    phases = db["phases"]
    index, step = _index_by_angle(phases)
    obs = [p for p in peaks if p["intensity"] >= MIN_OBS_INT]
    if not obs:
        obs = peaks
    cand = _candidate_phases(obs, index, step, tol)

    results = []
    for pi, obs_hits in cand.items():
        ph = phases[pi]
        if elements_only:
            if not set(elements_only) & set(ph["elements"]):
                continue
        sc = _score_phase(ph, obs, tol, shift=shift)
        if sc is None or sc["fom"] < min_fom:
            continue
        sc["phase"] = ph
        results.append(sc)

    results.sort(key=lambda r: -r["fom"])
    # прореживание почти одинаковых фаз (та же формула и близкий FOM)
    out, seen = [], set()
    for r in results:
        key = (r["phase"]["formula"], round(r["fom"] / 10.0))
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
        if len(out) >= top:
            break
    return out


def explain(peaks, db, tol=TOL_DEFAULT, min_fom=55.0, max_phases=MAX_PHASES,
            min_share=0.25, shift=0.0):
    """Жадно объяснить набор пиков минимальным набором фаз."""
    obs = [p for p in peaks if p["intensity"] >= MIN_OBS_INT]
    if not obs:
        obs = peaks
    if shift == "auto":
        shift, _gain = estimate_shift(peaks, db, tol=tol)
    shift = float(shift or 0.0)
    phases = db["phases"]
    index, step = _index_by_angle(phases)
    cand = _candidate_phases(obs, index, step, tol + abs(shift or 0.0))
    cand_list = list(cand.keys())

    unassigned = set(range(len(obs)))
    chosen = []
    for _ in range(max_phases):
        if not unassigned:
            break
        best = None
        for pi in cand_list:
            ph = phases[pi]
            sc = _score_phase(ph, obs, tol, allowed_obs=sorted(unassigned), shift=shift)
            if sc is None or sc["fom"] < min_fom:
                continue
            # доля объяснённых среди ещё не объяснённых пиков
            explained = {oi for _, oi, _ in sc["matched"]}
            share = sum(obs[i]["intensity"] for i in explained) / max(
                1e-9, sum(obs[i]["intensity"] for i in unassigned))
            if share < min_share:
                continue
            sc["share"] = share
            sc["phase"] = ph
            if best is None or sc["fom"] * share > best["fom"] * best["share"]:
                best = sc
        if best is None:
            break
        explained = {oi for _, oi, _ in best["matched"]}
        unassigned -= explained
        best["explained_peaks"] = sorted(explained)
        chosen.append(best)

    # итоговая привязка пик -> фаза
    per_peak = {}
    for ph_no, sc in enumerate(chosen, 1):
        for _, oi, dth in sc["matched"]:
            per_peak.setdefault(oi, []).append((ph_no, dth))
    elements = {}
    for sc in chosen:
        for el in sc["phase"]["elements"]:
            elements.setdefault(el, set()).add(sc["phase"]["formula"])
    return {
        "phases": chosen,
        "per_peak": per_peak,
        "unassigned": sorted(unassigned),
        "elements": {k: sorted(v) for k, v in elements.items()},
        "n_obs": len(obs),
        "shift": shift,
    }


def format_report(res, peaks, filename="", max_cand=3):
    """Текстовый отчёт для консоли и .txt."""
    lines = []
    if filename:
        lines.append("Файл: {}".format(filename))
    lines.append("Пиков найдено: {} (учтено при поиске фаз: {})".format(
        len(peaks), res["n_obs"]))
    sh = res.get("shift") or 0.0
    if abs(sh) >= 0.03:
        lines.append("Систематический сдвиг по 2θ (учтён): {:+.2f}°".format(sh))
    lines.append("")
    lines.append("ФАЗЫ (по убыванию достоверности):")
    if not res["phases"]:
        lines.append("  ничего не идентифицировано — ослабьте порог или расширьте базу")
    for i, sc in enumerate(res["phases"], 1):
        ph = sc["phase"]
        name = ph.get("name") or ph["formula"]
        lines.append(
            "  {}. {:<34s} {:<16s} FOM {:5.1f}  линий {}/{}  пиков объяснено {}".format(
                i, name[:34], ph["formula"][:16], sc["fom"],
                sc["n_matched"], sc["n_ref"], len(sc.get("explained_peaks", []))))
        lines.append("     элементы: {}   | источник: {} {}".format(
            ", ".join(ph["elements"]), ph.get("src", ""), ph.get("id", "")))
    lines.append("")
    if res["elements"]:
        lines.append("ЭЛЕМЕНТЫ (по идентифицированным фазам):")
        for el in sorted(res["elements"]):
            lines.append("  {:3s} — {}".format(el, ", ".join(res["elements"][el][:4])))
    else:
        lines.append("ЭЛЕМЕНТЫ: не определены")
    lines.append("")
    lines.append("ПИКИ -> ФАЗА:")
    for i, p in enumerate(peaks):
        ass = res["per_peak"].get(i)
        if ass:
            tag = ", ".join("#{} {}".format(n, res["phases"][n - 1]["phase"]["formula"])
                            for n, _ in sorted(ass)[:max_cand])
        else:
            tag = "—"
        fw = "  FWHM {:.3f}°".format(p["fwhm"]) if p.get("fwhm") else ""
        lines.append("  {:6.2f}°  d={:6.3f} Å  I={:5.1f}%{}   {}".format(
            p["two_theta"], p["d"], p["intensity"], fw, tag))
    if res["unassigned"]:
        lines.append("  не объяснены пики: {}".format(
            ", ".join("{:.2f}°".format(peaks[i]["two_theta"]) for i in res["unassigned"][:12])))
    return "\n".join(lines)
