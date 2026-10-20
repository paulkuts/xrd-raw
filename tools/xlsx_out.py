#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Сохранение дифрактограммы в Excel (.xlsx): таблица + встроенный график.

Лист «Данные»   — 2θ, d (Å), интенсивность.
Лист «График»   — родная диаграмма Excel по этим данным (её можно править в Excel).
Лист «Сводка»   — параметры съёмки и файла.
"""

from datetime import datetime
from pathlib import Path

from openpyxl import Workbook
from openpyxl.chart import LineChart, Reference
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from winxrd2xy import d_spacing, parse

HEADERS = ("2θ, °", "d, Å", "Интенсивность, имп.")


def make_xlsx(path, out=None, phases_result=None, peaks=None):
    data, step, meta = parse(path)
    src = Path(path)
    out = Path(out) if out else src.with_suffix(".xlsx")

    wb = Workbook()

    # --- лист с данными -----------------------------------------------------
    ws = wb.active
    ws.title = "Данные"
    ws.append(list(HEADERS))
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="1F4FA8")
        cell.alignment = Alignment(horizontal="center")

    for ang, inten in data:
        ws.append([ang, round(d_spacing(ang), 4), inten])

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
    rows = [
        ("Файл", src.name),
        ("Папка", str(src.parent)),
        ("Прибор / ПО", meta.get("software", "WinXRD")),
        ("Анод", meta.get("anode", "Cu") + " Kα, λ = 1.5406 Å"),
        ("Детектор", meta.get("detector", "—")),
        ("Атмосфера", meta.get("atmosphere", "—")),
        ("Точек", len(data)),
        ("2θ начало, °", round(data[0][0], 2)),
        ("2θ конец, °", round(data[-1][0], 2)),
        ("Шаг, °", step),
        ("Максимум, имп.", imax),
        ("Максимум при 2θ, °", round(top[0], 2)),
        ("d максимума, Å", round(d_spacing(top[0]), 4)),
        ("Создано", datetime.now().strftime("%d.%m.%Y %H:%M")),
    ]
    for key, value in rows:
        ws3.append([key, value])
    for row in ws3.iter_rows(min_row=1, max_row=ws3.max_row, max_col=2):
        row[0].font = Font(bold=True)
    ws3.column_dimensions["A"].width = 24
    ws3.column_dimensions["B"].width = 46

    # --- листы идентификации фаз -------------------------------------------
    if phases_result is not None and peaks:
        ws4 = wb.create_sheet("Пики")
        head4 = ("№", "2θ, °", "d, Å", "I, %", "FWHM, °", "Фаза", "Формула", "Элементы")
        ws4.append(list(head4))
        for cell in ws4[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4FA8")
            cell.alignment = Alignment(horizontal="center")
        res = phases_result
        for i, p in enumerate(peaks):
            ass = res["per_peak"].get(i)
            if ass:
                no = sorted(ass)[0][0]
                sc = res["phases"][no - 1]
                ph = sc["phase"]
                name = (ph.get("name") or "").split(" [ICSD")[0] or ph["formula"]
                ws4.append([p.get("n", i + 1), round(p["two_theta"], 3), round(p["d"], 4),
                            round(p["intensity"], 1),
                            round(p["fwhm"], 3) if p.get("fwhm") else None,
                            name, ph["formula"], ", ".join(ph["elements"])])
            else:
                ws4.append([p.get("n", i + 1), round(p["two_theta"], 3), round(p["d"], 4),
                            round(p["intensity"], 1),
                            round(p["fwhm"], 3) if p.get("fwhm") else None,
                            "—", "", ""])
        for col, w in zip("ABCDEFGH", (6, 10, 10, 9, 11, 30, 18, 22)):
            ws4.column_dimensions[col].width = w
        ws4.freeze_panes = "A2"

        ws5 = wb.create_sheet("Фазы")
        head5 = ("№", "Фаза", "Формула", "Элементы", "FOM", "Линий совпало",
                 "Линий в эталоне", "Пиков объяснено", "Источник", "ID")
        ws5.append(list(head5))
        for cell in ws5[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4FA8")
            cell.alignment = Alignment(horizontal="center")
        for i, sc in enumerate(res["phases"], 1):
            ph = sc["phase"]
            name = (ph.get("name") or "").split(" [ICSD")[0] or ph["formula"]
            ws5.append([i, name, ph["formula"], ", ".join(ph["elements"]),
                        round(sc["fom"], 1), sc["n_matched"], sc["n_ref"],
                        len(sc.get("explained_peaks", [])),
                        ph.get("src", ""), ph.get("id", "")])
        for col, w in zip("ABCDEFGHIJ", (5, 30, 18, 22, 8, 14, 15, 16, 12, 14)):
            ws5.column_dimensions[col].width = w
        ws5.freeze_panes = "A2"

        ws6 = wb.create_sheet("Элементы")
        ws6.append(["Элемент", "Встречается в фазах"])
        for cell in ws6[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4FA8")
            cell.alignment = Alignment(horizontal="center")
        for el in sorted(res.get("elements", {})):
            ws6.append([el, ", ".join(res["elements"][el])])
        ws6.column_dimensions["A"].width = 14
        ws6.column_dimensions["B"].width = 60

    wb.save(out)
    return out
