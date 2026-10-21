#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Идентификация фаз и элементов по пикам дифрактограммы .raw.

Связывает всё вместе: .raw -> пики -> база эталонов -> фазы, элементы,
график с подписями, отчёт .txt, листы Excel.

Единые источники:
  * излучение — tools/radiation.py (никаких собственных пересчётов и «1.5406 Å»);
  * индексы пиков и статусы — tools/result_view.py (структура docs/CONTRACT.md, п. 3–4);
  * отрисовка — tools/plot_xrd.py.
"""

import hashlib
import os
import sys
from datetime import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "tmp"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parent))

import peakpick  # noqa: E402
import phases  # noqa: E402
import plot_xrd  # noqa: E402
import radiation  # noqa: E402
import result_view  # noqa: E402
import xlsx_out  # noqa: E402
from winxrd2xy import parse  # noqa: E402

REPORT_WIDTH = 78

# Параметры съёмки, которые полезно показать в отчёте (ключ в meta -> подпись)
SURVEY_KEYS = (
    ("software", "Прибор / ПО"),
    ("anode", "Анод (в файле)"),
    ("detector", "Детектор"),
    ("atmosphere", "Атмосфера"),
    ("dwell", "Экспозиция, с"),
    ("step", "Шаг, °"),
    ("points", "Точек"),
    ("date", "Дата съёмки"),
    ("monochromator", "Монохроматор"),
    ("wavelength", "Длина волны из файла, Å"),
)

# То, чего в .raw обычно нет — эти параметры программа принимает по умолчанию
ASSUMED_KEYS = (
    ("monochromator", "монохроматор не указан"),
    ("date", "дата съёмки не указана"),
    ("wavelength", "длина волны берётся из анода, а не из файла"),
)


# ------------------------------------------------------------------ утилиты
def sha256_file(path):
    """SHA-256 исходного файла (потоково, без загрузки в память)."""
    digest = hashlib.sha256()
    with open(str(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detect_run_mode():
    """Как запущено: меню, CLI или из кода — попадает в отчёт."""
    name = Path(sys.argv[0]).stem.lower() if sys.argv and sys.argv[0] else ""
    if name in ("identify",):
        return "CLI (tools/identify.py)"
    if name in ("menu", "start", "__main__"):
        return "меню (tools/menu.py)"
    return "вызов из кода (библиотека)"


def _settings(tol, min_fom, shift_mode, use_profile, db_filter):
    return {
        "tol": float(tol),
        "min_fom": float(min_fom),
        "min_obs_int": float(getattr(phases, "MIN_OBS_INT", 3.0)),
        "shift_mode": str(shift_mode),
        "use_profile": bool(use_profile),
        "db_filter": str(db_filter),
    }


# ------------------------------------------------------------------- анализ
def radiation_short(rad):
    """Короткая строка излучения (для ans/сводок): без переносов и источника."""
    text = (rad.describe() or "").strip().splitlines()
    head = text[0] if text else ""
    return head.split(";")[0].strip()


def _pick(x, y, rad):
    """Пики через peakpick; старые версии без wavelength не ломают вызов."""
    try:
        return peakpick.pick(x, y, wavelength=rad.lam)
    except TypeError:
        return peakpick.pick(x, y)


def _explain(peaks, db, tol, min_fom, shift_mode, use_profile, db_filter):
    """Вызов phases.explain с новыми параметрами, с откатом к старым сигнатурам."""
    attempts = (
        dict(tol=tol, min_fom=min_fom, shift_mode=shift_mode,
             use_profile=use_profile, db_filter=db_filter),
        dict(tol=tol, min_fom=min_fom,
             shift="auto" if shift_mode == "auto" else 0.0),
        dict(tol=tol, min_fom=min_fom),
    )
    last = None
    for kwargs in attempts:
        try:
            return phases.explain(peaks, db, **kwargs)
        except (TypeError, ValueError) as exc:   # старый phases.explain не знает новых kwarg-ов
            last = exc
    raise last


def analyze(path, tol=phases.TOL_DEFAULT, min_fom=55.0, shift_mode="auto",
            use_profile=True, db_filter="all", anode=None, model=None):
    """Пики + фазы + элементы для файла .raw.

    Возвращает dict с ключами peaks/result/meta/data/step/db и добавленными
    radiation (объект Radiation), radiation_info (as_dict), sha256, settings.
    """
    data, step, meta = parse(path)
    x = [a for a, _ in data]
    y = [i for _, i in data]
    rad = radiation.Radiation.from_settings(meta, anode=anode, model=model)
    peaks = _pick(x, y, rad)
    db = phases.load_db()
    res = _explain(peaks, db, tol=tol, min_fom=min_fom, shift_mode=shift_mode,
                   use_profile=use_profile, db_filter=db_filter)
    if not isinstance(res, dict):
        raise TypeError("phases.explain вернул не dict: {!r}".format(type(res)))

    settings = _settings(tol, min_fom, shift_mode, use_profile, db_filter)
    # Какое пространство индексов вернул explain: контракт 1.0 (полный список)
    # или старая версия (отфильтрованный список). Метка нужна потребителям,
    # чтобы не пересчитывать индексы «на глаз» и не терять привязку фаз.
    contract_space = ("peak_status" in res) or ("ambiguous" in res)
    # дополняем структуру недостающими полями (старая версия phases)
    res.setdefault(result_view.SPACE_KEY, "full" if contract_space else "filtered")
    res.setdefault("settings", settings)
    res.setdefault("radiation", rad.as_dict())
    res.setdefault("shift_mode", shift_mode)
    res.setdefault("n_peaks", len(peaks))
    res.setdefault("elements", {})
    res.setdefault("phases", [])
    res.setdefault("per_peak", {})
    res.setdefault("unassigned", [])

    return {"peaks": peaks, "result": res, "meta": meta, "data": data,
            "step": step, "db": db, "radiation": rad,
            "radiation_info": dict(rad.as_dict(),
                                   label=radiation_short(rad),
                                   describe=rad.describe()),
            "sha256": sha256_file(path),
            "settings": settings}


# ------------------------------------------------------------------ график
def label_of(phase, limit=24):
    """Короткая подпись фазы для графика/отчёта."""
    return result_view.phase_name(phase, limit=limit)


def make_graph(ans, path, theme="dark", outdir=None, with_phases=True):
    """График с подписями пиков и фаз (делегируется в plot_xrd.make)."""
    png, svg, _meta = plot_xrd.make(
        path, theme=theme, outdir=outdir, with_phases=with_phases,
        phases_result=ans.get("result"), peaks=ans.get("peaks"),
        radiation=ans.get("radiation"))
    return png, svg


# -------------------------------------------------------------------- отчёт
def _fmt(value, spec="{:.2f}"):
    try:
        return spec.format(float(value))
    except (TypeError, ValueError):
        return "—"


def _survey_block(meta, rad):
    """Подтверждённые и предположенные параметры съёмки."""
    confirmed = []
    for key, title in SURVEY_KEYS:
        value = meta.get(key)
        if value not in (None, ""):
            confirmed.append((title, str(value)))
    titles = dict(SURVEY_KEYS)
    assumed = []
    for key, text in ASSUMED_KEYS:
        if meta.get(key) in (None, ""):
            assumed.append("{} ({}: {})".format(titles.get(key, key), key, text))
    notes = list(rad.notes or [])
    if not rad.confirmed and not any("анод" in n.lower() for n in notes):
        assumed.append("анод — по метаданным не определён, принято {} (ПРЕДПОЛОЖЕНИЕ)".format(rad.anode))
    if not any("модел" in n.lower() or "спектр" in n.lower() for n in notes):
        notes.append("спектральная модель — {} (задана явно)".format(radiation.MODELS[rad.model]))
    for note in notes:
        if note not in assumed:
            assumed.append(note)
    return confirmed, assumed


def _phase_table(scores):
    """Фазы в прежнем формате, расширенном статусом, источником, профилем и флагами."""
    lines = []
    for i, sc in enumerate(scores, 1):
        ph = sc.get("phase") or {}
        name = result_view.phase_name(ph, 34)
        lines.append(
            "  {}. {:<34s} {:<16s} FOM {:5.1f}  линий {}/{}  пиков объяснено {}".format(
                i, name, str(ph.get("formula", ""))[:16],
                float(sc.get("fom") or 0.0), sc.get("n_matched"), sc.get("n_ref"),
                len(sc.get("explained_peaks") or [])))
        extra = ["элементы: {}".format(", ".join(ph.get("elements") or []) or "—"),
                 "источник: {} {}".format(ph.get("src", ""), ph.get("id", "")).strip()]
        lines.append("     " + "   | ".join(extra))
        flags = result_view.phase_flags(sc)
        detail = ["статус: {}".format(result_view.phase_status(sc)),
                  "данные: {}".format(result_view.phase_datasource(sc) or "—")]
        if result_view.phase_prior(sc):
            detail.append("профиль: {}".format(", ".join(result_view.phase_prior(sc))))
        if flags:
            detail.append("флаги: {}".format(", ".join(flags)))
        if sc.get("mean_dtheta") is not None:
            detail.append("средний |Δ2θ| {:.3f}°".format(float(sc["mean_dtheta"])))
        n_ref_total = sc.get("n_ref_total")
        if n_ref_total is not None:
            detail.append("линий в эталоне {}, вне диапазона {}".format(
                n_ref_total, sc.get("n_ref_out_of_range", "—")))
        lines.append("     " + "   | ".join(detail))
        if "FEW_DIAGNOSTIC_LINES" in flags:
            lines.append("     ВНИМАНИЕ: в диапазоне скана мало диагностических линий — "
                         "фаза не может считаться подтверждённой без независимых данных")
    return lines


def _peaks_status_block(view, rad=None):
    """Список пиков со статусами и отклонениями 2θ (индексы полного списка)."""
    scores = view.phase_scores()
    codes = result_view.phase_codes([sc.get("phase") for sc in scores])
    lines = ["ПИКИ: СТАТУСЫ (Δ2θ — отклонение от ближайшего рефлекса фазы)"]
    lines.append("   {:>4s} {:>8s} {:>8s} {:>6s}  {:<15s} {:<34s} {:>8s}".format(
        "№", "2θ, °", "d, Å", "I, %", "статус", "фаза(ы)", "Δ2θ, °"))
    order = view.assigned + view.ambiguous + view.unexplained
    order.sort(key=lambda i: float(view.peaks[i]["two_theta"]))
    for i in order:
        peak = view.peaks[i]
        tt = float(peak["two_theta"])
        status = view.status_of(i)
        entries = view.candidates(i)
        names = []
        for entry in entries[:3]:
            no = entry["phase_no"]
            code = codes[no - 1] if 0 < no <= len(codes) else ""
            score = view.phase(no)
            names.append("{}{} {}".format(
                code, no, result_view.phase_name((score or {}).get("phase"), 16)).strip())
        best_dtheta = entries[0]["dtheta"] if entries else None
        d_value = rad.d(tt) if rad is not None else (peak.get("d") or 0.0)
        lines.append("   {:4d} {:8.2f} {:8.3f} {:6.1f}  {:<15s} {:<34s} {:>8s}".format(
            peak.get("n", i + 1), tt, float(d_value),
            float(peak.get("intensity") or 0.0),
            result_view.STATUS_RU.get(status, status),
            ("; ".join(names) if names else "—")[:34],
            "{:+.3f}".format(best_dtheta) if best_dtheta is not None else "—"))
    if view.below:
        lines.append("   под порогом поиска (below_threshold, НЕ считаются необъяснёнными), "
                     "{} шт.: {}".format(
                         len(view.below),
                         ", ".join("{:.2f}°".format(float(view.peaks[i]["two_theta"]))
                                   for i in view.below[:16])))
    counts = view.counts()
    lines.append("   итого: assigned {}, ambiguous {}, unexplained {}, below_threshold {}".format(
        counts.get(result_view.STATUS_ASSIGNED, 0), counts.get(result_view.STATUS_AMBIGUOUS, 0),
        counts.get(result_view.STATUS_UNEXPLAINED, 0), counts.get(result_view.STATUS_BELOW, 0)))
    return lines


LIMITATIONS = (
    "ОГРАНИЧЕНИЯ МЕТОДА",
    "  • Это предварительная идентификация фаз-кандидатов, а не количественный анализ.",
    "  • Совпадение 1–2 линий недостаточно для достоверной фазы; перекрытия пиков и слабые",
    "    линии могут остаться необъяснёнными — лучше пропуск, чем неверная подпись.",
    "  • Фазы JARVIS-DFT — расчётные структуры без экспериментального эталона (штраф в FOM).",
    "  • Профиль (prior) — мягкий бонус к ранжированию, не доказательство присутствия.",
    "  • Элементы рассчитаны по фазам-кандидатам; XRD сам по себе элементы не определяет.",
    "  • Сдвиг 2θ — диагностика юстировки; автоподбор не является подгонкой под фазу.",
    "  • Результат требует независимой проверки (эталон, литература, другой метод).",
    "",
    "FOM — оценка совпадения (0..100), НЕ вероятность: значения нельзя перемножать,",
    "складывать или трактовать как «процент уверенности».",
)


def build_report(ans, path, run_mode=None, db_meta=None, n_db=None, when=None):
    """Текст отчёта: расширенный формат (старые разделы сохранены)."""
    src = Path(path)
    res = ans["result"]
    peaks = ans["peaks"]
    meta = ans.get("meta") or {}
    rad = ans.get("radiation") or radiation.Radiation.from_settings(meta)
    settings = ans.get("settings") or res.get("settings") or {}
    view = result_view.PeakView(res, peaks)
    run_mode = run_mode or detect_run_mode()
    when = when or datetime.now().strftime("%d.%m.%Y %H:%M")
    if db_meta is None or n_db is None:
        db = ans.get("db")
        if db is not None:
            try:
                db_meta, n_db = phases.db_info(db)
            except Exception:                        # noqa: BLE001
                db_meta, n_db = db_meta or {}, n_db
        else:
            db_meta, n_db = db_meta or {}, n_db
    shift_mode, shift_value, shift_gain = result_view.shift_info(res)

    out = []
    add = out.append
    add("XRD — предварительная идентификация фаз (xrd-raw)")
    add("=" * REPORT_WIDTH)
    add("Файл:            {}".format(src.name))
    add("SHA-256:         {}".format(ans.get("sha256") or sha256_file(src)))
    add("Дата запуска:    {}".format(when))
    add("Режим запуска:   {}".format(run_mode))
    data = ans.get("data") or []
    add("Диапазон 2θ:     {:.2f}–{:.2f}°, шаг {:.3f}°, точек {}".format(
        data[0][0] if data else 0.0, data[-1][0] if data else 0.0,
        float(ans.get("step") or 0.0), len(data)))
    add("Излучение:       {}".format(rad.describe()))
    add("                 база посчитана на {}".format(
        "ту же длину волны" if rad.matches_base() else
        "λ = {} Å — позиции эталонов пересчитаны по d".format(radiation.BASE_WAVELENGTH)))

    confirmed, assumed = _survey_block(meta, rad)
    add("")
    add("ПАРАМЕТРЫ СЪЁМКИ — подтверждены файлом:")
    if confirmed:
        for title, value in confirmed:
            add("  {:<26s} {}".format(title + ":", value))
    else:
        add("  (в файле нет пригодных метаданных)")
    add("ПАРАМЕТРЫ СЪЁМКИ — предположения (НЕ подтверждено файлом):")
    if assumed:
        for text in assumed:
            add("  - {}".format(text))
    else:
        add("  (нет)")

    add("")
    add("НАСТРОЙКИ ДЕТЕКЦИИ И СОПОСТАВЛЕНИЯ:")
    add("  пиков найдено:            {} (из них участвовало в поиске: {})".format(
        view.n, res.get("n_obs", "—")))
    add("  мин. интенсивность пика:  {} % от максимума (пики слабее — below_threshold)".format(
        _fmt(settings.get("min_obs_int"), "{:g}")))
    add("  допуск по 2θ:             ±{}°".format(_fmt(settings.get("tol"), "{:g}")))
    add("  порог FOM:                {}".format(_fmt(settings.get("min_fom"), "{:g}")))
    shift_text = ("{:+.3f}°".format(shift_value) if abs(shift_value) >= 0.005
                  else "не применялся (0.000°)")
    add("  сдвиг нуля 2θ:            {} ({}{})".format(
        shift_mode, shift_text,
        ", выигрыш {:+.2f}".format(float(shift_gain)) if shift_gain else ""))
    add("  профиль (prior):          {}".format(
        "включён — мягкий бонус к ранжированию, не фильтр"
        if settings.get("use_profile", True) else "выключен"))
    add("  отбор базы:               {}".format(settings.get("db_filter") or "all"))
    if db_meta or n_db:
        add("  база эталонов:            {} фаз, версия {}, {}".format(
            n_db or db_meta.get("count", "—"), db_meta.get("version", "—"),
            db_meta.get("wavelength", "—")))
        if db_meta.get("sources"):
            add("  происхождение базы:       {}".format(db_meta["sources"]))
        if db_meta.get("from_COD") or db_meta.get("from_JARVIS"):
            add("                            COD: {}, JARVIS-DFT: {}".format(
                db_meta.get("from_COD", "—"), db_meta.get("from_JARVIS", "—")))

    # --- прежний формат: пики, фазы, элементы (разделы сохранены) ----------
    add("")
    add("-" * REPORT_WIDTH)
    add("Пиков найдено: {} (учтено при поиске фаз: {})".format(view.n, res.get("n_obs", "—")))
    if abs(shift_value) >= 0.03:
        add("Систематический сдвиг по 2θ (учтён): {:+.2f}°".format(shift_value))
    add("")
    add("ФАЗЫ (по убыванию достоверности):")
    scores = view.phase_scores()
    if scores:
        out.extend(_phase_table(scores))
    else:
        add("  ничего не идентифицировано — ослабьте порог или расширьте базу")

    add("")
    elements = res.get("elements") or {}
    if elements:
        add("ЭЛЕМЕНТЫ (по идентифицированным фазам, расчётно — не прямое определение):")
        for el in sorted(elements):
            add("  {:3s} — {}".format(el, ", ".join(elements[el][:4])))
    else:
        add("ЭЛЕМЕНТЫ: не определены")

    add("")
    out.extend(_peaks_status_block(view, rad))

    add("")
    add("-" * REPORT_WIDTH)
    out.extend(LIMITATIONS)
    add("")
    add("Отчёт сформирован автоматически; решение о фазовом составе принимает специалист.")
    return "\n".join(out)


def save_report(ans, path, outdir=None, run_mode=None):
    """Записать <stem>_phases.txt (utf-8-sig) и вернуть (путь, текст)."""
    src = Path(path)
    outdir = Path(outdir) if outdir else src.parent
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    txt = build_report(ans, path, run_mode=run_mode)
    out = outdir / "{}_phases.txt".format(src.stem)
    out.write_text(txt, encoding="utf-8-sig")
    return out, txt


# ---------------------------------------------------------------------- run
def run(path, theme="dark", both_themes=False, tol=phases.TOL_DEFAULT,
        min_fom=55.0, want_xlsx=True, outdir=None, with_phases=True,
        shift_mode="auto", use_profile=True, db_filter="all", anode=None,
        model=None):
    """Полный проход: анализ + отчёт + график(и) + Excel.

    outdir=None — писать рядом с исходником (прежнее поведение).
    with_phases=False — график без подписей фаз (чистый), отчёт и XLSX как обычно.
    """
    ans = analyze(path, tol=tol, min_fom=min_fom, shift_mode=shift_mode,
                  use_profile=use_profile, db_filter=db_filter,
                  anode=anode, model=model)
    outdir = Path(outdir) if outdir else None
    files = []
    themes = ("dark", "light") if both_themes else (theme,)
    for th in themes:
        png, svg = make_graph(ans, path, theme=th, outdir=outdir, with_phases=with_phases)
        files += [png, svg]
    rep, txt = save_report(ans, path, outdir=outdir)
    files.append(rep)
    if want_xlsx:
        target = (outdir / "{}.xlsx".format(Path(path).stem)) if outdir else None
        xlsx = xlsx_out.make_xlsx(path, out=target, phases_result=ans["result"],
                                  peaks=ans["peaks"], radiation=ans["radiation"],
                                  sha256=ans["sha256"], settings=ans["settings"],
                                  db_meta=phases.db_info(ans["db"])[0],
                                  n_phases_db=phases.db_info(ans["db"])[1],
                                  run_mode=detect_run_mode())
        files.append(xlsx)
    return ans, files, txt


def main():
    """CLI: без меню, для проверки и пакетной работы.

    python tools\\identify.py input\\sample.raw [--theme paper] [--both] [--tol 0.3]
    """
    import argparse
    try:
        from console import setup as _setup_console

        _setup_console()
    except Exception:  # noqa: BLE001
        for stream in (sys.stdout, sys.stderr):
            try:
                stream.reconfigure(errors="replace")
            except Exception:  # noqa: BLE001
                pass
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--theme", default="dark", choices=sorted(plot_xrd.THEMES))
    ap.add_argument("--both", action="store_true")
    ap.add_argument("--tol", type=float, default=phases.TOL_DEFAULT)
    ap.add_argument("--min-fom", type=float, default=55.0)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--no-phases", dest="with_phases", action="store_false",
                    help="график без подписей фаз (чистый)")
    ap.add_argument("--shift", default="auto", choices=("auto", "off", "manual"))
    ap.add_argument("--no-profile", dest="use_profile", action="store_false")
    ap.add_argument("--db-filter", default="all")
    ap.add_argument("--anode", default=None)
    ap.add_argument("--model", default=None, choices=sorted(radiation.MODELS))
    a = ap.parse_args()
    for f in a.files:
        p = Path(f)
        if not p.is_file():
            print("нет файла: {}".format(f))
            continue
        print("=" * 60)
        ans, files, txt = run(str(p), theme=a.theme, both_themes=a.both, tol=a.tol,
                              min_fom=a.min_fom, outdir=a.outdir,
                              with_phases=a.with_phases, shift_mode=a.shift,
                              use_profile=a.use_profile, db_filter=a.db_filter,
                              anode=a.anode, model=a.model)
        print(txt)
        print()
        for x in files:
            print("  {}".format(x))


if __name__ == "__main__":
    main()
