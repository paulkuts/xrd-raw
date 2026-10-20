# -*- coding: utf-8 -*-
"""Дымовой тест портатива: парсер, графики, Excel, идентификация."""
import os
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parent
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp"))
sys.stdout.reconfigure(errors="replace")

import identify  # noqa: E402
import plot_xrd  # noqa: E402
import winxrd2xy  # noqa: E402
import xlsx_out  # noqa: E402

src = ROOT / "input" / "sample1.raw"
print("импорты ок")

pts, step, meta = winxrd2xy.parse(src)
print("точек в файле:", len(pts), "шаг", step, "диапазон", pts[0][0], "-", pts[-1][0])

ans, res, txt = identify.run(src, theme="dark")
print("идентификация: файлов", len(res), "| фаз найдено", len(ans["result"].get("phases", [])))
for f in res:
    print("  ", Path(f).name, Path(f).stat().st_size, "байт")

xy = winxrd2xy.write_xy(pts, src, ROOT / "tmp" / "smoke.xy")
csv = winxrd2xy.write_csv(pts, src, ROOT / "tmp" / "smoke.csv")
print("xy/csv:", Path(xy).stat().st_size, Path(csv).stat().st_size)

png, svg, _ = plot_xrd.make(src, theme="light", outdir=ROOT / "tmp")
print("график light:", Path(png).stat().st_size, Path(svg).stat().st_size)

xl = xlsx_out.make_xlsx(src, ROOT / "tmp" / "smoke.xlsx")
print("xlsx:", Path(xl).stat().st_size)

import openpyxl  # noqa: E402

wb = openpyxl.load_workbook(str(xl))
print("листы Excel:", wb.sheetnames)
print("ВСЁ ОК")
