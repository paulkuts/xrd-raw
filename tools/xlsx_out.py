#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Сохранение дифрактограммы в Excel (.xlsx): таблица + встроенный график.

Листы:
  «Данные»   — 2θ, d (Å), интенсивность (d считается в едином излучении);
  «График»   — родная диаграмма Excel по этим данным (её можно править в Excel);
  «Сводка»   — параметры съёмки, излучение, сдвиг, настройки, SHA-256, база;
  «Пики»     — все пики со статусами assigned/ambiguous/unexplained/below_threshold
               и всеми фазами-претендентами;
  «Фазы»     — кандидаты: FOM, линии (в диапазоне/всего/вне диапазона),
               источник данных (COD/JARVIS-DFT), профиль, флаги;
  «Элементы» — расчётно по фазам-кандидатам (не прямое определение).

Индексы пиков берутся из результата phases.explain() через result_view:
никаких пересчётов на стороне Excel (docs/CONTRACT.md, п. 3–4).
"""

import hashlib
from datetime import datetime
from pathlib import Path

from openpyxl import Workbook
from openpyxl.chart import LineChart, Reference
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

import result_view
from winxrd2xy import parse

HEADERS = ("2θ, °", "d, Å", "Интенсивность, имп.")
HDR_FILL = "1F4FA8"
STATUS_NOTE = ("Статус: assigned — пик объяснён фазой; ambiguous — на пик претендуют "
               ">=2 фаз (перекрытие); unexplained — не объяснён; "
               "below_threshold — слабее порога поиска, не считается необъяснённым")


def _sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _radiation(radiation, meta):
    try:
        import radiation as _rad
    except Exception:                                   # noqa: BLE001
        return None
    if radiation is not None and hasattr(radiation, "lam"):
        return radiation
    if isinstance(radiation, dict) and radiation.get("anode"):
        try:
            return _rad.Radiation.from_settings(meta or {}, anode=radiation.get("anode"),
                                                model=radiation.get("model"))
        except Exception:                               # noqa: BLE001
            pass
    try:
        return _rad.Radiation.from_settings(meta or {})
    except Exception:                                   # noqa: BLE001
        return None


def _header(ws, titles, widths):
    ws.append(list(titles))
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor=HDR_FILL)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for col, width in zip("ABCDEFGHIJKLMNOPQRSTUV", widths):
        ws.column_dimensions[col].width = width
    ws.freeze_panes = "A2"


def _formula_name(phase):
    phase = phase or {}
    name = (phase.get("name") or "").split(" [ICSD")[0].strip() or phase.get("formula", "")
    return name


def _phase_dict(view, phase_no):
    """Запись фазы из базы по 1-based номеру (без падений на старых структурах)."""
    score = view.phase(phase_no)
    return (score or {}).get("phase") or {}


def make_xlsx(path, out=None, phases_result=None, peaks=None, radiation=None,
              sha256=None, settings=None, db_meta=None, n_phases_db=None,
              run_mode=None):
    """Записать .xlsx. Сигнатура совместима: make_xlsx(path[, out], phases_result=, peaks=)."""
    data, step, meta = parse(path)
    src = Path(path)
    out = Path(out) if out else src.with_suffix(".xlsx")
    rad = _radiation(radiation, meta)
    settings = dict(settings or (phases_result or {}).get("settings") or {})
    view = None
    if phases_result is not None and peaks is not None:
        view = result_view.PeakView(phases_result, peaks)
    if sha256 is None:
        try:
            sha256 = _sha256(src)
        except OSError:
            sha256 = ""

    def d_of(two_theta):
        """Единая длина волны во всех листах — только через radiation."""
        if rad is not None:
            return rad.d(two_theta)
        try:
            import radiation as _rad
            return _rad.d_from_two_theta(two_theta)
        except Exception:                               # noqa: BLE001
            return float("nan")

    wb = Workbook()

    # --- лист с данными -----------------------------------------------------
    ws = wb.active
    ws.title = "Данные"
    ws.append(list(HEADERS))
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor=HDR_FILL)
        cell.alignment = Alignment(horizontal="center")

    for ang, inten in data:
        ws.append([ang, round(d_of(ang), 4), inten])

    ws.column_dimensions["A"].width = 12
    ws.column_dimensions["B"].width = 12
    ws.column_dimensions["C"].width = 20
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, max_col=3):
        row[0].number_format = "0.00"
        row[1].number_format = "0.0000"
        row[2].number_format = "0"
    ws.freeze_panes = "A2"

    # --- родная диаграмма Excel --------------------------------------------
    ws2 = wb.create_sheet("График")
    chart = LineChart()
    chart.title = "Дифрактограмма: {}".format(src.stem)
    chart.style = 2
    chart.x_axis.title = "2θ, градусы"
    chart.y_axis.title = "Интенсивность, имп."
    chart.height = 11
    chart.width = 26
    chart.add_data(Reference(ws, min_col=3, min_row=1, max_row=ws.max_row), titles_from_data=True)
    chart.set_categories(Reference(ws, min_col=1, min_row=2, max_row=ws.max_row))
    chart.series[0].graphicalProperties.line.width = 12000   # 1.2 pt
    chart.series[0].smooth = False
    ws2.add_chart(chart, "B2")

    # --- сводка -------------------------------------------------------------
    ws3 = wb.create_sheet("Сводка")
    imax = max(i for _, i in data)
    top = max(data, key=lambda p: p[1])
    if rad is not None:
        rad_text = rad.describe()
        rad_short = "{} {}, λ = {:.4f} Å".format(
            rad.anode, "Kα1" if getattr(rad, "model", "") == "ka1" else "Kα", rad.lam)
    else:
        rad_text = "Cu Kα, λ = 1.5418 Å (предположение: модуль radiation недоступен)"
        rad_short = "Cu Kα, λ = 1.5418 Å"
    shift_mode, shift_value, shift_gain = result_view.shift_info(phases_result or {})
    rows = [
        ("Файл", src.name),
        ("Папка", str(src.parent)),
        ("SHA-256 исходного файла", sha256 or "—"),
        ("Прибор / ПО", meta.get("software", "WinXRD")),
        ("Излучение", rad_text),
        ("Длина волны, Å", round(rad.lam, 5) if rad is not None else 1.54184),
        ("Анод подтверждён", "да" if (rad is not None and rad.confirmed) else "НЕТ (предположение)"),
        ("Детектор", meta.get("detector", "—")),
        ("Атмосфера", meta.get("atmosphere", "—")),
        ("Точек", len(data)),
        ("2θ начало, °", round(data[0][0], 2)),
        ("2θ конец, °", round(data[-1][0], 2)),
        ("Шаг, °", step),
        ("Максимум, имп.", imax),
        ("Максимум при 2θ, °", round(top[0], 2)),
        ("d максимума, Å", round(d_of(top[0]), 4)),
        ("Сдвиг 2θ (режим)", shift_mode),
        ("Сдвиг 2θ, °", round(shift_value, 3)),
    ]
    if shift_gain:
        rows.append(("Выигрыш от сдвига, %", round(float(shift_gain), 2)))
    rows += [
        ("Допуск по 2θ, °", settings.get("tol", "—")),
        ("Порог FOM", settings.get("min_fom", "—")),
        ("Мин. интенсивность пика, %", settings.get("min_obs_int", "—")),
        ("Профиль (prior)", "вкл" if settings.get("use_profile", True) else "выкл"),
        ("Отбор базы", settings.get("db_filter", "—")),
        ("Режим запуска", run_mode or "—"),
    ]
    if n_phases_db:
        rows.append(("База эталонов, фаз", n_phases_db))
    if db_meta:
        rows.append(("Источник базы", db_meta.get("sources", "")))
        rows.append(("Примечание к базе", db_meta.get("note", "")))
        if db_meta.get("wavelength") or db_meta.get("base_wavelength"):
            rows.append(("Длина волны, в которой посчитана база",
                         db_meta.get("wavelength") or db_meta.get("base_wavelength")))
    if view is not None:
        counts = view.counts()
        rows += [
            ("Пиков всего", view.n),
            ("Пиков assigned", counts.get(result_view.STATUS_ASSIGNED, 0)),
            ("Пиков ambiguous (перекрытие)", counts.get(result_view.STATUS_AMBIGUOUS, 0)),
            ("Пиков unexplained", counts.get(result_view.STATUS_UNEXPLAINED, 0)),
            ("Пиков below_threshold", counts.get(result_view.STATUS_BELOW, 0)),
            ("Фаз-кандидатов", len(view.phase_scores())),
        ]
    rows += [
        ("Создано", datetime.now().strftime("%d.%m.%Y %H:%M")),
        ("FOM", "оценка совпадения (0..100), НЕ вероятность"),
        ("Статусы пиков", STATUS_NOTE),
    ]
    for key, value in rows:
        ws3.append([key, value])
    for row in ws3.iter_rows(min_row=1, max_row=ws3.max_row, max_col=2):
        row[0].font = Font(bold=True)
        row[1].alignment = Alignment(wrap_text=True, vertical="top")
    ws3.column_dimensions["A"].width = 30
    ws3.column_dimensions["B"].width = 72

    # --- листы идентификации фаз -------------------------------------------
    if view is not None:
        scores = view.phase_scores()
        codes = result_view.phase_codes([sc.get("phase") for sc in scores])

        # --- «Пики» ---------------------------------------------------------
        ws4 = wb.create_sheet("Пики")
        head4 = ("№", "2θ, °", "d, Å", "I, %", "FWHM, °", "Статус", "Фаза", "Код",
                 "Формула", "Элементы", "Фазы-претенденты", "Δ2θ лучшей, °")
        _header(ws4, head4, (6, 9, 9, 8, 9, 15, 28, 6, 16, 20, 44, 13))

        def _cand_text(entries):
            parts = []
            for entry in entries:
                no = entry["phase_no"]
                phase = _phase_dict(view, no)
                code = codes[no - 1] if 0 < no <= len(codes) else "?"
                parts.append("{} {} (Δ{:+.3f}°)".format(
                    code, result_view.phase_name(phase, 22), entry["dtheta"]))
            return "; ".join(parts)

        for i, p in enumerate(view.peaks):
            tt = float(p["two_theta"])
            status = view.status_of(i)
            entries = view.candidates(i)
            best = view.best(i)
            row = [p.get("n", i + 1), round(tt, 3), round(d_of(tt), 4),
                   round(float(p.get("intensity") or 0.0), 1),
                   round(p["fwhm"], 3) if p.get("fwhm") else None, status]
            if best is not None:
                phase = _phase_dict(view, best["phase_no"])
                code = codes[best["phase_no"] - 1] if 0 < best["phase_no"] <= len(codes) else ""
                row += [_formula_name(phase), code, phase.get("formula", ""),
                        ", ".join(phase.get("elements") or []),
                        _cand_text(entries), round(best["dtheta"], 3)]
            else:
                row += ["—", "", "", "", "", None]
            ws4.append(row)
        for row in ws4.iter_rows(min_row=2, max_row=ws4.max_row, max_col=len(head4)):
            row[1].number_format = "0.000"
            row[2].number_format = "0.0000"
            row[3].number_format = "0.0"
            row[4].number_format = "0.000"
            row[11].number_format = "+0.000"
        ws4.auto_filter.ref = "A1:{}{}".format(get_column_letter(len(head4)), ws4.max_row)

        # --- «Фазы» ---------------------------------------------------------
        ws5 = wb.create_sheet("Фазы")
        head5 = ("№", "Фаза", "Формула", "Элементы", "FOM", "Статус", "Линий совпало",
                 "Линий в диапазоне", "Линий вне диапазона", "Пиков объяснено",
                 "Средний Δ2θ, °", "Источник данных", "Профиль (prior)", "Флаги", "ID")
        _header(ws5, head5, (5, 30, 16, 18, 8, 12, 13, 14, 15, 14, 13, 15, 20, 30, 14))
        for i, sc in enumerate(scores, 1):
            ph = sc.get("phase") or {}
            n_ref = sc.get("n_ref")
            n_ref_total = sc.get("n_ref_total")
            n_ref_out = sc.get("n_ref_out_of_range")
            if n_ref_out is None and n_ref_total is not None and n_ref is not None:
                n_ref_out = max(0, int(n_ref_total) - int(n_ref))
            ws5.append([i, _formula_name(ph), ph.get("formula", ""),
                        ", ".join(ph.get("elements") or []),
                        round(float(sc.get("fom") or 0.0), 1),
                        result_view.phase_status(sc),
                        sc.get("n_matched"), n_ref, n_ref_out,
                        len(sc.get("explained_peaks") or []),
                        round(float(sc.get("mean_dtheta") or 0.0), 3),
                        result_view.phase_datasource(sc) or ph.get("src", ""),
                        ", ".join(result_view.phase_prior(sc)) or "—",
                        ", ".join(result_view.phase_flags(sc)) or "—",
                        ph.get("id", "")])
        for row in ws5.iter_rows(min_row=2, max_row=ws5.max_row, max_col=len(head5)):
            row[4].number_format = "0.0"
            row[10].number_format = "0.000"
        ws5.auto_filter.ref = "A1:{}{}".format(get_column_letter(len(head5)), ws5.max_row)

        # --- «Элементы» -----------------------------------------------------
        ws6 = wb.create_sheet("Элементы")
        ws6.append(["Элемент (расчётно по фазам-кандидатам)", "Встречается в фазах"])
        for cell in ws6[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor=HDR_FILL)
            cell.alignment = Alignment(horizontal="center")
        for el in sorted((view.result.get("elements") or {})):
            ws6.append([el, ", ".join(view.result["elements"][el])])
        ws6.column_dimensions["A"].width = 38
        ws6.column_dimensions["B"].width = 60

    if out.parent and not out.parent.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out))
    return out
