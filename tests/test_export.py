#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Тесты экспорта: TXT / XLSX / PNG / SVG на синтетическом образце.

Ничего не читается из input/ и не нужна сеть: .raw генерируется во временной
папке (tempfile.mkdtemp) в формате WinXRD, база фаз — маленькая синтетическая
(полная база на 32k фаз в тестах не грузится).

Проверяется:
  * единый набор пиков и фаз в TXT, XLSX и подписях графика (по 2θ и индексам);
  * единая длина волны во всех листах XLSX (1.5418, а не 1.5406);
  * PNG >= 300 dpi, SVG остаётся векторным;
  * имена файлов прежние, при указанном outdir рядом с исходником ничего не пишется;
  * устойчивость к старой структуре phases.explain (индексы отфильтрованного
    списка, per_peak кортежами) — регрессия на «уезжающие» подписи фаз.

Запуск: /root/work/.venv/bin/python -m unittest discover -s tests -v
"""

import hashlib
import math
import os
import random
import re
import shutil
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

import identify  # noqa: E402  (ставит MPLCONFIGDIR до импорта matplotlib)
import plot_xrd  # noqa: E402
import radiation  # noqa: E402
import result_view  # noqa: E402
import xlsx_out  # noqa: E402
from winxrd2xy import d_spacing, parse  # noqa: E402

START, STEP, N = 10.0, 0.02, 3001
# (2θ, амплитуда, σ) — кварц, кальцит, портландит
SYNTH_PEAKS = ((18.06, 55.0, 0.10), (20.85, 95.0, 0.10), (26.65, 40.0, 0.10),
               (29.40, 100.0, 0.10), (34.10, 26.0, 0.11), (36.55, 22.0, 0.11),
               (39.41, 58.0, 0.11), (43.15, 46.0, 0.11), (47.12, 20.0, 0.12),
               (47.50, 42.0, 0.12), (50.15, 24.0, 0.12))

# Эталонные линии (d, Å): кварц / кальцит / портландит
REF_LINES = {
    "Quartz": ("SiO2", [4.257, 3.343, 2.457, 1.818]),
    "Calcite": ("CaCO3", [3.035, 2.285, 2.095, 1.913]),
    "Portlandite": ("Ca(OH)2", [4.918, 2.628, 1.927]),
}

ANGLE_TOL = 0.011          # допуск сравнения углов TXT/XLSX/график, °


# --------------------------------------------------------------- синтетика
def write_raw(path, seed=42):
    """.raw в формате WinXRD: текстовый заголовок + записи <iid> (I, dwell, 2θ)."""
    rng = random.Random(seed)
    chunks = [struct.pack("<iid", 0, N, START)]
    for k in range(N):
        ang = START + k * STEP
        y = 45.0 + 5.0 * math.sin(ang * 0.9) + rng.gauss(0.0, 1.5)
        for pos, amp, width in SYNTH_PEAKS:
            y += amp * math.exp(-0.5 * ((ang - pos) / width) ** 2)
        chunks.append(struct.pack("<iid", int(max(0, round(y))), 150, ang))
    header = b"WinXRD 2.0-6  Copper  Solid State  air\n"
    Path(path).write_bytes(header + b"".join(chunks))
    return path


def two_theta_of_d(d, lam=radiation.BASE_WAVELENGTH):
    return 2.0 * math.degrees(math.asin(lam / (2.0 * d)))


def synthetic_db():
    """Маленькая база: три фазы с рефлексами, посчитанными по d."""
    items = []
    for name, (formula, ds) in REF_LINES.items():
        peaks = []
        for k, d in enumerate(ds):
            peaks.append([round(two_theta_of_d(d), 4), d,
                          100.0 if k == 0 else max(6.0, 70.0 - 12 * k)])
        items.append({
            "src": "COD", "id": "9000001", "name": name, "formula": formula,
            "elements": sorted(set(re.findall(r"[A-Z][a-z]?", formula))),
            "sg": "P 1", "peaks": peaks, "prio": 1,
        })
    return {"meta": {"version": "test", "wavelength": "CuKa 1.5418 A",
                     "count": len(items), "sources": "COD (test fixture)",
                     "from_COD": len(items), "from_JARVIS": 0},
            "phases": items}


# ------------------------------------------------------------------ утилиты
def png_info(path):
    """(ширина, высота, dpi_x) из заголовка PNG."""
    raw = Path(path).read_bytes()
    assert raw[:8] == b"\x89PNG\r\n\x1a\n", "не PNG"
    width = height = dpi = None
    pos = 8
    while pos + 8 <= len(raw):
        length = struct.unpack(">I", raw[pos:pos + 4])[0]
        kind = raw[pos + 4:pos + 8]
        body = raw[pos + 8:pos + 8 + length]
        if kind == b"IHDR":
            width, height = struct.unpack(">II", body[:8])
        elif kind == b"pHYs":
            px, _py, unit = struct.unpack(">IIB", body[:9])
            if unit == 1:
                dpi = round(px * 0.0254)
        pos += 12 + length
    return width, height, dpi


def txt_status_rows(txt):
    """Строки таблицы «ПИКИ: СТАТУСЫ» -> [(2θ, статус)]."""
    rows = []
    pattern = re.compile(r"^\s{3,}\d+\s+(\d+\.\d{2})\s+(\d+\.\d{3})\s+(\d+\.\d)\s{2}(\S.*?)\s{2,}")
    for line in txt.splitlines():
        m = pattern.match(line)
        if m:
            rows.append((float(m.group(1)), m.group(4).strip()))
    return rows


class TempCase(unittest.TestCase):
    """Общая обвязка: временная папка и синтетический .raw."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="xrd_export_"))
        cls.srcdir = cls.tmp / "src"
        cls.srcdir.mkdir()
        cls.raw = write_raw(cls.srcdir / "synthetic.raw")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(str(cls.tmp), ignore_errors=True)

    def assertAnglesClose(self, got, expected, tol=ANGLE_TOL):
        self.assertEqual(len(got), len(expected),
                         "разное число пиков: {} vs {}".format(len(got), len(expected)))
        for a, b in zip(sorted(got), sorted(expected)):
            self.assertLess(abs(a - b), tol, "2θ расходятся: {} vs {}".format(a, b))


# ------------------------------------------------------------------- тесты
class TestFullExport(TempCase):
    """Полный проход identify.run на синтетике (без полной базы)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = cls.tmp / "out"
        cls.db = synthetic_db()
        with mock.patch.object(identify.phases, "load_db", return_value=cls.db):
            cls.ans, cls.files, cls.txt = identify.run(
                str(cls.raw), theme="paper", outdir=str(cls.outdir), shift_mode="off")
        cls.names = sorted(p.name for p in cls.files)

    def test_source_dir_untouched(self):
        """При указанном outdir рядом с исходником не появляется ни одного файла."""
        self.assertEqual(sorted(os.listdir(str(self.srcdir))), ["synthetic.raw"])

    def test_output_names_compatible(self):
        """Имена файлов прежние: <stem>_phases_<theme>.png/.svg, _phases.txt, .xlsx."""
        for expected in ("synthetic_phases_paper.png", "synthetic_phases_paper.svg",
                         "synthetic_phases.txt", "synthetic.xlsx"):
            self.assertIn(expected, self.names)
        for path in self.files:
            self.assertTrue(Path(path).is_file(), path)
            self.assertEqual(Path(path).parent, self.outdir)

    def test_txt_required_content(self):
        """TXT: SHA-256, дата/режим, излучение, настройки, статусы, ограничения, FOM."""
        txt = self.txt
        self.assertIn(hashlib.sha256(self.raw.read_bytes()).hexdigest(), txt)
        self.assertIn("Дата запуска:", txt)
        self.assertIn("Режим запуска:", txt)
        self.assertIn("Излучение:", txt)
        self.assertIn("λ = 1.5418 Å", txt)
        self.assertIn("НАСТРОЙКИ ДЕТЕКЦИИ И СОПОСТАВЛЕНИЯ", txt)
        self.assertIn("допуск по 2θ", txt)
        self.assertIn("ПАРАМЕТРЫ СЪЁМКИ — подтверждены файлом", txt)
        self.assertIn("ПАРАМЕТРЫ СЪЁМКИ — предположения", txt)
        self.assertIn("ПИКИ: СТАТУСЫ", txt)
        self.assertIn("ОГРАНИЧЕНИЯ МЕТОДА", txt)
        self.assertIn("FOM — оценка совпадения", txt)
        self.assertIn("НЕ вероятность", txt)
        self.assertIn("происхождение базы", txt)
        self.assertIn("assigned", txt)
        self.assertIn("below_threshold", txt)
        # прежние разделы отчёта сохранены
        self.assertIn("Пиков найдено:", txt)
        self.assertIn("ФАЗЫ (по убыванию достоверности):", txt)
        self.assertIn("ЭЛЕМЕНТЫ", txt)

    def test_txt_encoding(self):
        """Файл отчёта — utf-8-sig."""
        raw = (self.outdir / "synthetic_phases.txt").read_bytes()
        self.assertTrue(raw.startswith(b"\xef\xbb\xbf"))
        raw.decode("utf-8-sig")

    def test_xlsx_sheets_and_chart(self):
        """Состав листов и родная диаграмма Excel сохранены."""
        from openpyxl import load_workbook

        wb = load_workbook(str(self.outdir / "synthetic.xlsx"))
        self.assertEqual(wb.sheetnames,
                         ["Данные", "График", "Сводка", "Пики", "Фазы", "Элементы"])
        self.assertEqual(len(wb["График"]._charts), 1)

    def test_xlsx_single_wavelength(self):
        """d во всех листах считается в той же длине волны, что и в radiation."""
        from openpyxl import load_workbook

        _data, _step, meta = parse(str(self.raw))
        rad = radiation.Radiation.from_settings(meta)
        self.assertAlmostEqual(rad.lam, 1.54184, places=5)
        wb = load_workbook(str(self.outdir / "synthetic.xlsx"))
        for row in list(wb["Данные"].iter_rows(min_row=2, max_row=6, values_only=True)):
            ang, d_cell = float(row[0]), float(row[1])
            expected = rad.lam / (2.0 * math.sin(math.radians(ang / 2.0)))
            self.assertAlmostEqual(d_cell, expected, places=4)
            # именно 1.5418, а не 1.5406 (иначе d уезжает на ~0.008 %)
            self.assertNotAlmostEqual(d_cell, d_spacing(ang), places=6)

        texts = {str(r[0]): r[1] for r in wb["Сводка"].iter_rows(values_only=True)
                 if r and r[0]}
        self.assertIn("Излучение", texts)
        self.assertIn("λ = 1.5418 Å", str(texts["Излучение"]))
        self.assertAlmostEqual(float(texts["Длина волны, Å"]), rad.lam, places=5)
        self.assertEqual(texts["SHA-256 исходного файла"],
                         hashlib.sha256(self.raw.read_bytes()).hexdigest())

    def test_xlsx_peaks_statuses(self):
        """Лист «Пики»: колонка статуса и все фазы-претенденты."""
        from openpyxl import load_workbook

        wb = load_workbook(str(self.outdir / "synthetic.xlsx"))
        head = [c.value for c in wb["Пики"][1]]
        for column in ("Статус", "Фазы-претенденты", "Δ2θ лучшей, °"):
            self.assertIn(column, head)
        rows = list(wb["Пики"].iter_rows(min_row=2, values_only=True))
        self.assertEqual(len(rows), len(self.ans["peaks"]))
        statuses = {r[head.index("Статус")] for r in rows}
        self.assertTrue(statuses <= {"assigned", "ambiguous", "unexplained",
                                     "below_threshold"}, statuses)
        self.assertIn("assigned", statuses)
        self.assertTrue(any(r[head.index("Фазы-претенденты")] for r in rows),
                        "у объяснённых пиков должен быть список фаз-претендентов")

    def test_xlsx_phases_columns(self):
        """Лист «Фазы»: FOM, статус, n_ref, datasource, prior, flags, ID."""
        from openpyxl import load_workbook

        wb = load_workbook(str(self.outdir / "synthetic.xlsx"))
        head = [c.value for c in wb["Фазы"][1]]
        for column in ("FOM", "Статус", "Линий в диапазоне", "Источник данных",
                       "Профиль (prior)", "Флаги", "ID"):
            self.assertIn(column, head)
        rows = list(wb["Фазы"].iter_rows(min_row=2, values_only=True))
        self.assertEqual(len(rows), len(self.ans["result"]["phases"]))
        self.assertGreaterEqual(len(rows), 2, "синтетика должна дать >=2 фазы")

    def test_consistency_txt_xlsx_graph(self):
        """Одни и те же пики (по 2θ) и фазы во всех трёх представлениях."""
        from openpyxl import load_workbook

        wb = load_workbook(str(self.outdir / "synthetic.xlsx"))
        head = [c.value for c in wb["Пики"][1]]
        status_col = head.index("Статус")
        xlsx_peaks = [float(r[1]) for r in wb["Пики"].iter_rows(min_row=2, values_only=True)
                      if r[status_col] != result_view.STATUS_BELOW]
        txt_peaks = [angle for angle, status in txt_status_rows(self.txt)
                     if status != result_view.STATUS_RU[result_view.STATUS_BELOW]]
        labels = list(plot_xrd.LAST["labels"])

        self.assertTrue(xlsx_peaks, "в синтетике должны быть объяснённые пики")
        self.assertAnglesClose(txt_peaks, xlsx_peaks)
        self.assertTrue(labels, "на графике не оказалось подписей пиков")
        self.assertLessEqual(len(labels), plot_xrd.MAX_LABELS)
        for angle in labels:
            nearest = min(xlsx_peaks, key=lambda a: abs(a - angle))
            self.assertLess(abs(nearest - angle), ANGLE_TOL,
                            "подпись {} не соответствует ни одному пику XLSX".format(angle))

        # фазы: одинаковое число и одинаковые формулы в TXT и XLSX
        xlsx_formulas = [r[2] for r in wb["Фазы"].iter_rows(min_row=2, values_only=True)]
        scores = self.ans["result"]["phases"]
        self.assertEqual(len(xlsx_formulas), len(scores))
        for formula in xlsx_formulas:
            self.assertIn(formula, self.txt)
        codes = plot_xrd.LAST["codes"]
        self.assertEqual(len(codes), len(scores))
        self.assertEqual(len(set(codes)), len(codes), "коды фаз должны быть уникальны")

    def test_phase_labels_land_on_matching_peaks(self):
        """Подпись фазы стоит на пике, где действительно есть её рефлекс.

        Это и есть регрессия P0: при сдвиге индексов подписи уезжали на чужие пики.
        """
        view = result_view.PeakView(self.ans["result"], self.ans["peaks"])
        self.assertEqual(view.n, len(self.ans["peaks"]))
        checked = 0
        for i in view.assigned + view.ambiguous:
            tt = float(view.peaks[i]["two_theta"])
            for entry in view.candidates(i):
                refs = [float(p[0]) for p in
                        view.phase(entry["phase_no"])["phase"]["peaks"]]
                nearest = min(refs, key=lambda ref: abs(ref - tt))
                self.assertLess(abs(nearest - tt), 0.30 + 1e-6,
                                "фаза {} подписана на пике {:.2f}°, ближайший рефлекс {:.2f}°"
                                .format(entry["phase_no"], tt, nearest))
                checked += 1
        self.assertGreater(checked, 0)

    def test_counts_agree_with_n_obs(self):
        """below_threshold + объяснённые + необъяснённые = все; n_obs сходится."""
        view = result_view.PeakView(self.ans["result"], self.ans["peaks"])
        counts = view.counts()
        self.assertEqual(sum(counts.values()), view.n)
        self.assertEqual(view.n - counts[result_view.STATUS_BELOW],
                         self.ans["result"]["n_obs"])


class TestPlotOutputs(TempCase):
    """PNG/SVG: dpi, векторность, чистый график без фаз, тема paper."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = cls.tmp / "plot"

    def test_png_300dpi_and_svg_vector(self):
        png, svg, _meta = plot_xrd.make(str(self.raw), theme="paper",
                                        outdir=str(self.outdir))
        width, _height, dpi = png_info(png)
        self.assertGreaterEqual(dpi, 300, "PNG должен быть не ниже 300 dpi")
        self.assertGreater(width, 3000)

        text = Path(svg).read_text(encoding="utf-8")
        self.assertIn("<svg", text)
        self.assertNotIn("<image", text, "SVG не должен содержать растровых вставок")
        self.assertNotIn("base64", text)
        self.assertGreater(text.count("<path"), 20, "SVG должен быть векторным")
        self.assertTrue("pt" in text or "mm" in text)

    def test_paper_theme_light_background(self):
        """Тема paper: белый фон, тонкая тёмная кривая; dark/light на месте."""
        theme = plot_xrd.THEMES["paper"]
        self.assertEqual(theme["bg"].lower(), "#ffffff")
        self.assertLessEqual(theme["curve_lw"], 1.0)
        self.assertIn("dark", plot_xrd.THEMES)
        self.assertIn("light", plot_xrd.THEMES)

    def test_clean_graph_without_phases(self):
        """with_phases=False — чистый график: без Bragg-полос и кодов фаз."""
        png, _svg, _meta = plot_xrd.make(str(self.raw), theme="dark",
                                         outdir=str(self.outdir), with_phases=False)
        self.assertTrue(png.name.endswith("_graph_dark.png"))
        self.assertEqual(plot_xrd.LAST["n_bragg"], 0)
        self.assertEqual(plot_xrd.LAST["codes"], [])

    def test_paper_theme_used_by_identify(self):
        """identify.run(theme='paper') пишет <stem>_phases_paper.png и SVG."""
        db = synthetic_db()
        out = self.tmp / "paper_run"
        with mock.patch.object(identify.phases, "load_db", return_value=db):
            _ans, _files, _txt = identify.run(str(self.raw), theme="paper",
                                              outdir=str(out), shift_mode="off")
        png = out / "synthetic_phases_paper.png"
        self.assertTrue(png.is_file())
        self.assertGreaterEqual(png_info(png)[2], 300)
        self.assertTrue((out / "synthetic_phases_paper.svg").is_file())


class TestLegacyStructure(TempCase):
    """Старая структура phases.explain: индексы отфильтрованного списка.

    Регрессия на P0: раньше make_graph строил idx_map по id(dict) и брал
    per_peak[i] по индексу полного списка, из-за чего подписи фаз уезжали
    на чужие пики.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        data, _step, meta = parse(str(cls.raw))
        cls.peaks = identify._pick([a for a, _ in data], [i for _, i in data],
                                   radiation.Radiation.from_settings(meta))
        # индексы отфильтрованного списка (как это делал старый explain)
        cls.obs = [i for i, p in enumerate(cls.peaks) if p["intensity"] >= 3.0]
        weak = [i for i, p in enumerate(cls.peaks) if p["intensity"] < 3.0]
        cls.weak_index = weak[0] if weak else None
        cls.legacy = cls._legacy_result(cls.peaks, cls.obs)

    @staticmethod
    def _legacy_result(peaks, obs):
        """Результат в старом формате: per_peak кортежами, индексы obs-списка."""
        obs_pos = {full: k for k, full in enumerate(obs)}
        spec = (("Quartz", "SiO2", [4.257, 3.343]),
                ("Calcite", "CaCO3", [3.035]),
                ("Vaterite", "CaCO3", [3.030]))     # делит пик 29.4° с кальцитом
        scores, per_peak = [], {}
        for no, (name, formula, ds) in enumerate(spec, 1):
            refs = [two_theta_of_d(d) for d in ds]
            matched = []
            for ri, ref in enumerate(refs):
                best = None
                for oi in obs:
                    dth = abs(peaks[oi]["two_theta"] - ref)
                    if dth <= 0.30 and (best is None or dth < best[0]):
                        best = (dth, oi)
                if best is not None:
                    matched.append((ri, obs_pos[best[1]], round(best[0], 4)))
            for _ri, oi, dth in matched:
                per_peak.setdefault(oi, []).append((no, dth))
            scores.append({
                "phase": {"formula": formula, "name": name, "src": "COD", "id": str(no),
                          "elements": sorted(set(re.findall(r"[A-Z][a-z]?", formula))),
                          "peaks": [[round(ref, 4), d, 100.0] for ref, d in zip(refs, ds)]},
                "fom": 90.0 - 10 * (no - 1), "fom_raw": 90.0, "n_matched": len(matched),
                "n_ref": len(refs), "mean_dtheta": 0.05, "matched": matched,
                "explained_peaks": sorted({oi for _ri, oi, _d in matched})})
        return {
            "phases": scores, "per_peak": per_peak,
            "unassigned": [i for i in range(len(obs)) if i not in per_peak],
            "elements": {"Si": ["SiO2"], "Ca": ["CaCO3"]},
            "n_obs": len(obs), "shift": 0.0,
        }

    def _view(self):
        return result_view.PeakView(self.legacy, self.peaks)

    def test_space_detected_as_filtered(self):
        view = self._view()
        self.assertTrue(view.legacy)
        self.assertEqual(view.space, "filtered")
        self.assertLess(len(self.obs), view.n)
        self.assertEqual(view.n - len(view.below), len(self.obs))

    def test_indices_map_to_full_list(self):
        """per_peak[obs_index] -> индекс полного списка, а не сдвинутый."""
        view = self._view()
        # пик 29.4°: индекс в полном списке != индекс в отфильтрованном
        c1 = min(range(len(self.peaks)),
                 key=lambda i: abs(self.peaks[i]["two_theta"] - 29.40))
        self.assertNotEqual(self.obs.index(c1), c1)
        self.assertEqual(set(view.phases_of(c1)), {2, 3})
        self.assertTrue(view.is_ambiguous(c1))
        # каждый ключ per_peak (индекс отфильтрованного списка) попал на свой пик
        for obs_index in self.legacy["per_peak"]:
            full = self.obs[obs_index]
            self.assertTrue(view.phases_of(full),
                            "потеряна привязка пика obs[{}] -> full[{}]".format(obs_index, full))

    def test_below_threshold_not_unexplained(self):
        view = self._view()
        if self.weak_index is None:
            self.skipTest("в синтетике нет пиков ниже порога")
        self.assertEqual(view.status_of(self.weak_index), result_view.STATUS_BELOW)
        self.assertNotIn(self.weak_index, view.unexplained)
        self.assertEqual(view.counts()[result_view.STATUS_BELOW], view.n - len(self.obs))

    def test_phase_labels_land_on_matching_peaks(self):
        """Подпись фазы стоит на пике с её рефлексом (а не на соседнем по индексу)."""
        view = self._view()
        for i in view.assigned + view.ambiguous:
            tt = float(view.peaks[i]["two_theta"])
            for entry in view.candidates(i):
                refs = [float(p[0]) for p in view.phase(entry["phase_no"])["phase"]["peaks"]]
                nearest = min(refs, key=lambda ref: abs(ref - tt))
                self.assertLess(abs(nearest - tt), 0.30 + 1e-6)

    def test_report_and_graph_survive_old_structure(self):
        """TXT/XLSX/график строятся и на старой структуре (без новых ключей)."""
        data, _step, meta = parse(str(self.raw))
        rad = radiation.Radiation.from_settings(meta)
        ans = {"peaks": self.peaks, "result": self.legacy, "meta": meta, "data": data,
               "step": STEP, "db": None, "radiation": rad,
               "sha256": identify.sha256_file(self.raw), "settings": {}}
        txt = identify.build_report(ans, self.raw, run_mode="тест", db_meta={}, n_db=0)
        self.assertIn("перекрытие", txt)
        self.assertIn("под порогом поиска", txt)

        xlsx = xlsx_out.make_xlsx(str(self.raw), out=str(self.tmp / "legacy.xlsx"),
                                  phases_result=self.legacy, peaks=self.peaks,
                                  radiation=rad)
        self.assertTrue(Path(xlsx).is_file())

        png, svg, _meta = plot_xrd.make(str(self.raw), theme="light",
                                        outdir=str(self.tmp / "legacy_plot"),
                                        with_phases=True, phases_result=self.legacy,
                                        peaks=self.peaks, radiation=rad)
        self.assertTrue(Path(png).is_file() and Path(svg).is_file())
        self.assertGreaterEqual(png_info(png)[2], 300)
        self.assertGreater(plot_xrd.LAST["n_bragg"], 0, "Bragg-полосы должны быть нарисованы")
        self.assertLessEqual(len(plot_xrd.LAST["labels"]), plot_xrd.MAX_LABELS)
        self.assertTrue(set(plot_xrd.LAST["codes"]) <= {"Q", "C", "V"},
                        plot_xrd.LAST["codes"])


class TestExplicitSpaceMarker(TempCase):
    """Метка index_space: analyze помечает пространство индексов результата."""

    def test_analyze_marks_space(self):
        db = synthetic_db()
        with mock.patch.object(identify.phases, "load_db", return_value=db):
            ans = identify.analyze(str(self.raw), shift_mode="off")
        self.assertIn(result_view.SPACE_KEY, ans["result"])
        self.assertIn(ans["result"][result_view.SPACE_KEY], ("full", "filtered"))
        view = result_view.PeakView(ans["result"], ans["peaks"])
        self.assertEqual(view.space, ans["result"][result_view.SPACE_KEY])

    def test_contract_structure_kept_as_full_space(self):
        """Результат с peak_status/ambiguous читается как полный список."""
        peaks = [{"two_theta": 10.0, "intensity": 1.0, "d": 8.8},
                 {"two_theta": 20.0, "intensity": 100.0, "d": 4.4},
                 {"two_theta": 30.0, "intensity": 50.0, "d": 3.0}]
        res = {"phases": [{"phase": {"formula": "X", "name": "X"}, "fom": 80}],
               "per_peak": {1: [{"phase_no": 1, "dtheta": 0.01, "ref_int": 100.0,
                                 "overlap": False}]},
               "unassigned": [2], "ambiguous": [],
               "peak_status": {1: "assigned", 2: "unexplained"},
               "settings": {"min_obs_int": 3.0}}
        view = result_view.PeakView(res, peaks)
        self.assertEqual(view.space, "full")
        self.assertFalse(view.legacy)
        self.assertEqual(view.phases_of(1), [1])
        self.assertEqual(view.status_of(0), result_view.STATUS_BELOW)   # I = 1 % < 3 %
        self.assertNotIn(0, view.unexplained)


class TestPublicSignatures(unittest.TestCase):
    """Публичные сигнатуры: старые аргументы на месте, новые — kwarg с дефолтом."""

    def test_identify_run(self):
        import inspect

        params = inspect.signature(identify.run).parameters
        for name in ("path", "theme", "both_themes", "tol", "min_fom", "want_xlsx",
                     "outdir", "with_phases", "shift_mode", "use_profile",
                     "db_filter", "anode", "model"):
            self.assertIn(name, params)
            if name != "path":
                self.assertIsNot(params[name].default, inspect.Parameter.empty, name)
        self.assertEqual(params["theme"].default, "dark")
        self.assertIs(params["outdir"].default, None)

    def test_identify_analyze(self):
        import inspect

        params = inspect.signature(identify.analyze).parameters
        for name in ("path", "tol", "min_fom", "shift_mode", "use_profile",
                     "db_filter", "anode", "model"):
            self.assertIn(name, params)

    def test_plot_xrd_make(self):
        import inspect

        params = inspect.signature(plot_xrd.make).parameters
        for name in ("path", "theme", "outdir", "with_phases", "phases_result",
                     "peaks", "radiation"):
            self.assertIn(name, params)
        self.assertIs(params["with_phases"].default, False)

    def test_xlsx_make(self):
        import inspect

        params = inspect.signature(xlsx_out.make_xlsx).parameters
        for name in ("path", "out", "phases_result", "peaks"):
            self.assertIn(name, params)

    def test_themes_present(self):
        self.assertEqual(set(("dark", "light", "paper")) - set(plot_xrd.THEMES), set())

    def test_result_view_codes_stable(self):
        codes = result_view.phase_codes([{"name": "Quartz low", "formula": "SiO2"},
                                         {"name": "Calcite", "formula": "CaCO3"},
                                         {"name": "Portlandite", "formula": "Ca(OH)2"}])
        self.assertEqual(codes, ["Q", "C", "P"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
