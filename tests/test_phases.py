#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Тесты tools/phases.py и tools/building_materials.py (unittest, без сети/SciPy).

База для тестов — МАЛЕНЬКАЯ СИНТЕТИЧЕСКАЯ, собирается здесь же в gzip-json.
Полная база (32k фаз) — только один тест под ``XRD_SLOW=1``.

Покрываемые требования контракта (docs/CONTRACT.md, п.4-5):
  * индексы результата — в пространстве полного списка пиков (слабый пик перед сильным);
  * два пика одной фазы могут принадлежать разным фазам (перекрытие, ambiguous);
  * рефлексы вне измеренного диапазона не штрафуют фазу;
  * отрицательный контроль: случайное совпадение одной линии не даёт кандидата;
  * пустой / почти пустой спектр;
  * автосдвиг находит внесённые 0.1° и сообщает его в результате и отчёте;
  * JARVIS-DFT помечается и штрафуется, фильтр db_filter='cod' работает;
  * профиль стройматериалов — мягкий бонус, не фильтр.
"""

import gzip
import json
import math
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import building_materials  # noqa: E402
import phases  # noqa: E402
import radiation  # noqa: E402

LAM = radiation.BASE_WAVELENGTH
TOL = 0.25


# ------------------------------------------------------------------ утилиты
def d_of(two_theta):
    return radiation.d_from_two_theta(two_theta, LAM)


def ref(two_theta, intensity):
    """Рефлекс базы в формате [2θ, d, I%] — как в phases_db.json.gz."""
    return [round(two_theta, 3), round(d_of(two_theta), 4), float(intensity)]


def make_phase(pid, formula, lines, elements=("X",), src="COD", prio=1, name=None):
    return {
        "src": src,
        "id": pid,
        "name": name or pid,
        "formula": formula,
        "elements": list(elements),
        "sg": "P 1",
        "prio": prio,
        "peaks": [ref(tt, inten) for tt, inten in lines],
    }


def make_peaks(spec, shift=0.0):
    """spec: [(2θ, I%)] -> список пиков в формате peakpick.pick()."""
    out = []
    for i, (tt, inten) in enumerate(sorted(spec, key=lambda r: r[0])):
        pos = round(tt + shift, 4)
        out.append({
            "two_theta": pos,
            "d": round(d_of(pos), 4),
            "intensity": float(inten),
            "height": float(inten),
            "prominence": float(inten),
            "fwhm": 0.08,
            "n": i + 1,
            "peak_id": i,
            "snr": 10.0,
        })
    return out


def make_db(tmpdir, phases_list, name="mini_db.json.gz"):
    """Записать синтетическую базу в gzip-json и загрузить её через load_db."""
    path = Path(tmpdir) / name
    payload = {
        "meta": {"version": "test", "wavelength": "CuKa 1.5418 A",
                 "count": len(phases_list), "note": "synthetic test db"},
        "phases": phases_list,
    }
    with gzip.open(str(path), "wt", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)
    return phases.load_db(path)


# Линии фаз мини-базы (2θ, I%). Позиции специально разнесены, кроме 20.00°
# у ALPHA/BETA — это общая линия для проверки перекрытия.
ALPHA_LINES = [(20.00, 100), (25.00, 70), (30.00, 50), (35.00, 30), (40.00, 20)]
BETA_LINES = [(20.00, 100), (28.00, 80), (33.00, 60), (38.00, 40), (45.00, 25)]
GAMMA_LINES = [(10.00, 100), (15.00, 80), (18.00, 60)]
NARROW_LINES = [(8.00, 100), (12.00, 90), (16.00, 80), (40.00, 50), (44.00, 40), (48.00, 30)]
NARROW2_LINES = [(40.00, 50), (44.00, 40), (48.00, 30)]
SOLO_LINES = [(5.00, 90), (30.00, 100)]
QUARTZ_LINES = [(26.66, 100), (20.88, 22), (36.54, 8), (39.47, 5), (50.14, 12), (45.79, 4)]
COMPETITOR_LINES = [(26.66, 100), (20.88, 22), (36.54, 8), (39.47, 5), (50.14, 12), (45.79, 4)]


def mini_phases():
    return [
        make_phase("COD 1", "Aa2Bb", ALPHA_LINES, elements=("A", "B"), name="Alpha"),
        make_phase("COD 2", "Cc3Dd", BETA_LINES, elements=("C", "D"), name="Beta"),
        make_phase("COD 3", "EeFf3", GAMMA_LINES, elements=("E", "F"), name="Gamma"),
        make_phase("COD 4", "Gg2Hh5", NARROW_LINES, elements=("G",), name="Narrow"),
        make_phase("COD 5", "Gg2Hh5b", NARROW2_LINES, elements=("G",), name="Narrow2"),
        make_phase("COD 6", "IiJj", SOLO_LINES, elements=("I",), name="Solo"),
        make_phase("COD 7", "Aa2BbJ", ALPHA_LINES, elements=("A", "B"),
                   src="JARVIS-DFT", prio=2, name="AlphaJarvis"),
        # реальный id из профиля стройматериалов (кварц) — для проверки приора
        make_phase("COD 1011097", "SiO2", QUARTZ_LINES, elements=("Si", "O"),
                   name="Quartz low"),
        # «конкурент» с тем же паттерном, но не из профиля
        make_phase("COD 8", "Kk2Ll3", COMPETITOR_LINES, elements=("K",), name="Competitor"),
    ]


class MiniDbCase(unittest.TestCase):
    """Общая синтетическая база во временной папке."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="xrd_phases_")
        cls.db = make_db(cls.tmpdir, mini_phases())

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def by_id(self, res):
        return {sc["phase"]["id"]: sc for sc in res["phases"]}

    def ids(self, res):
        return [sc["phase"]["id"] for sc in res["phases"]]

    def status_of(self, res, phase_id):
        sc = self.by_id(res).get(phase_id)
        return sc["status"] if sc else None

    def flags_of(self, res, phase_id):
        sc = self.by_id(res).get(phase_id)
        return list(sc["flags"]) if sc else []


# --------------------------------------------------------------- 1. индексы
class TestPeakIndexSpace(MiniDbCase):
    """Регрессия дефекта: слабый пик перед сильными сдвигал подписи на 1."""

    def setUp(self):
        # СЛАБЫЙ пик 10.00° (I=2% < min_obs_int=3%) стоит ПЕРЕД сильными:
        # ALPHA даёт 5 линий 20/25/30/35/40°.
        self.peaks = make_peaks([(10.00, 2.0)] + ALPHA_LINES)
        self.res = phases.explain(self.peaks, self.db, tol=TOL, min_fom=40.0,
                                  shift_mode="off")

    def test_counts(self):
        self.assertEqual(6, len(self.peaks))
        self.assertEqual(6, self.res["n_peaks"])
        self.assertEqual(5, self.res["n_obs"])
        self.assertEqual(1, self.res["n_below_threshold"])

    def test_per_peak_uses_full_list_indices(self):
        keys = sorted(self.res["per_peak"])
        self.assertEqual([1, 2, 3, 4, 5], keys)       # 0 — слабый пик, не участвовал
        self.assertEqual("below_threshold", self.res["peak_status"][0])
        self.assertNotIn(0, self.res["unassigned"])

    def test_matched_and_explained_in_full_space(self):
        sc = self.res["phases"][0]
        self.assertEqual([1, 2, 3, 4, 5], sc["explained_peaks"])
        for _ref_index, peak_index, _dtheta in sc["matched"]:
            self.assertIn(peak_index, (1, 2, 3, 4, 5))

    def test_report_labels_not_shifted(self):
        txt = phases.format_report(self.res, self.peaks, filename="synthetic.raw")
        weak = [ln for ln in txt.splitlines() if ln.strip().startswith("10.00")]
        strong = [ln for ln in txt.splitlines() if ln.strip().startswith("20.00")]
        self.assertTrue(weak, "в отчёте нет строки пика 10.00°")
        self.assertIn("ниже порога", weak[0])
        self.assertNotIn("#1", weak[0])
        self.assertTrue(strong, "в отчёте нет строки пика 20.00°")
        self.assertIn("#1", strong[0])

    def test_unexplained_peak_is_not_lost(self):
        peaks = make_peaks([(10.00, 2.0)] + ALPHA_LINES + [(48.42, 60.0)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=40.0, shift_mode="off")
        self.assertEqual([6], res["unassigned"])
        self.assertEqual("unexplained", res["peak_status"][6])
        self.assertNotIn(0, res["unassigned"])


# ------------------------------------------------------------- 3. перекрытие
class TestOverlap(MiniDbCase):
    """Пик 20.00° принадлежит и ALPHA, и BETA; ни одна фаза не исчезает."""

    def setUp(self):
        # отдельная база из двух фаз: без JARVIS-двойника, чтобы перекрытие
        # было ровно между ALPHA и BETA
        self.db2 = make_db(self.tmpdir, [make_phase("COD 1", "Aa2Bb", ALPHA_LINES),
                                         make_phase("COD 2", "Cc3Dd", BETA_LINES)],
                           name="overlap_db.json.gz")
        spec = ALPHA_LINES + [ln for ln in BETA_LINES if ln[0] != 20.00]
        self.peaks = make_peaks(spec)
        self.res = phases.explain(self.peaks, self.db2, tol=TOL, min_fom=45.0,
                                  shift_mode="off")

    def test_both_phases_found(self):
        ids = self.ids(self.res)
        self.assertIn("COD 1", ids)
        self.assertIn("COD 2", ids)

    def test_shared_peak_marked_ambiguous_and_overlap(self):
        self.assertIn(0, self.res["ambiguous"])          # 20.00° — общая линия
        entries = self.res["per_peak"][0]
        self.assertEqual(2, len(entries))
        self.assertTrue(all(e["overlap"] for e in entries))
        self.assertEqual({1, 2}, {e["phase_no"] for e in entries})
        self.assertEqual("ambiguous", self.res["peak_status"][0])

    def test_no_exclusive_assignment(self):
        # у обеих фаз в matched есть пик 20.00° (индекс 0 полного списка)
        for sc in self.res["phases"]:
            if sc["phase"]["id"] in ("COD 1", "COD 2"):
                self.assertIn(0, [p for _r, p, _d in sc["matched"]])

    def test_unique_peaks_counted(self):
        for sc in self.res["phases"]:
            if sc["phase"]["id"] in ("COD 1", "COD 2"):
                self.assertGreater(sc["n_unique_peaks"], 0)
                self.assertIn(sc["status"], ("candidate", "ambiguous"))


# --------------------------------------------------------- 2. диапазон скана
class TestNarrowRange(MiniDbCase):
    """Линии эталона вне измеренного диапазона не штрафуют FOM."""

    def test_out_of_range_lines_do_not_change_score(self):
        peaks = make_peaks([(8.00, 100), (40.00, 50), (44.00, 40), (48.00, 30)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=40.0, shift_mode="off",
                             two_theta_min=38.0, two_theta_max=50.0)
        full = self.by_id(res)["COD 4"]        # 3 сильные линии вне диапазона
        only = self.by_id(res)["COD 5"]        # те же линии, но их нет в карточке
        self.assertEqual(3, full["n_ref_out_of_range"])
        self.assertEqual(0, only["n_ref_out_of_range"])
        self.assertEqual(3, full["n_ref"])
        self.assertEqual(3, only["n_ref"])
        self.assertAlmostEqual(full["fom"], only["fom"], places=6)
        self.assertAlmostEqual(full["cov_ref"], only["cov_ref"], places=6)

    def test_peak_outside_declared_range_is_reported(self):
        peaks = make_peaks([(8.00, 100), (40.00, 50), (44.00, 40), (48.00, 30)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=40.0, shift_mode="off",
                             two_theta_min=38.0, two_theta_max=50.0)
        self.assertEqual(1, res["n_outside_range"])
        self.assertEqual("outside_range", res["peak_status"][0])
        self.assertEqual(3, res["n_obs"])
        self.assertNotIn(0, res["unassigned"])

    def test_few_diagnostic_lines_flag(self):
        # SOLO: одна линия в диапазоне 25-35° -> «недостаточно данных»
        peaks = make_peaks([(30.02, 100)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=20.0, shift_mode="off",
                             two_theta_min=25.0, two_theta_max=35.0)
        sc = self.by_id(res)["COD 6"]
        self.assertIn("FEW_DIAGNOSTIC_LINES", sc["flags"])
        self.assertEqual("weak", sc["status"])
        self.assertNotEqual("candidate", sc["status"])


# ------------------------------------------------------- отрицательный контроль
class TestNegativeControl(MiniDbCase):
    def test_single_random_line_is_not_candidate(self):
        # GAMMA имеет линии 10/15/18; в спектре только случайный пик 18.03°
        peaks = make_peaks([(18.03, 100.0), (46.10, 40.0), (52.30, 25.0), (57.00, 15.0)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=55.0, shift_mode="off")
        sc = self.by_id(res).get("COD 3")
        if sc is not None:
            self.assertEqual("weak", sc["status"])
            self.assertIn(sc["status"], ("weak", "ambiguous"))
        self.assertNotIn("candidate", [s["status"] for s in res["phases"]
                                       if s["phase"]["id"] == "COD 3"])

    def test_loose_threshold_still_marks_single_line(self):
        peaks = make_peaks([(18.03, 100.0), (46.10, 40.0), (52.30, 25.0), (57.00, 15.0)])
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=10.0, shift_mode="off")
        sc = self.by_id(res)["COD 3"]
        self.assertIn("SINGLE_LINE_MATCH", sc["flags"])
        self.assertEqual("weak", sc["status"])
        self.assertLessEqual(sc["fom"], phases.SINGLE_LINE_CAP)


# ------------------------------------------------------------ пустой спектр
class TestEmptySpectrum(MiniDbCase):
    def test_empty_list(self):
        res = phases.explain([], self.db)
        self.assertEqual([], res["phases"])
        self.assertEqual({}, res["per_peak"])
        self.assertEqual([], res["unassigned"])
        self.assertEqual([], res["ambiguous"])
        self.assertEqual({}, res["peak_status"])
        self.assertEqual(0, res["n_obs"])
        self.assertEqual(0, res["n_peaks"])
        self.assertEqual(0.0, res["shift"])
        txt = phases.format_report(res, [], filename="empty.raw")
        self.assertIn("ФАЗЫ", txt)
        self.assertIn("НЕ вероятность", txt)

    def test_all_peaks_below_threshold(self):
        peaks = make_peaks([(20.0, 1.0), (25.0, 2.0), (30.0, 0.5)])
        res = phases.explain(peaks, self.db)
        self.assertEqual([], res["phases"])
        self.assertEqual(0, res["n_obs"])
        self.assertEqual(3, res["n_below_threshold"])
        self.assertEqual(
            {0: "below_threshold", 1: "below_threshold", 2: "below_threshold"},
            res["peak_status"])
        self.assertEqual([], res["unassigned"])

    def test_single_peak_spectrum(self):
        res = phases.explain(make_peaks([(30.0, 100.0)]), self.db, min_fom=55.0,
                             shift_mode="off")
        for sc in res["phases"]:
            self.assertNotEqual("candidate", sc["status"])


# ------------------------------------------------------------------ 4. сдвиг
class TestShift(MiniDbCase):
    SHIFT = 0.10

    def shifted_peaks(self):
        spec = ALPHA_LINES + [ln for ln in BETA_LINES if ln[0] != 20.00] + NARROW2_LINES
        return make_peaks(spec, shift=self.SHIFT)

    def test_auto_finds_shift(self):
        peaks = self.shifted_peaks()
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=45.0)
        self.assertEqual("auto", res["shift_mode"])
        self.assertAlmostEqual(self.SHIFT, res["shift"], delta=0.051)
        self.assertGreaterEqual(res["shift_gain"], phases.SHIFT_GAIN)
        self.assertTrue(res["shift_detail"]["accepted"])

    def test_auto_shift_reported_in_text(self):
        peaks = self.shifted_peaks()
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=45.0)
        txt = phases.format_report(res, peaks)
        self.assertIn("Сдвиг 2θ", txt)
        self.assertIn("+0.10", txt)
        self.assertIn("auto", txt)

    def test_shift_off(self):
        peaks = self.shifted_peaks()
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=45.0, shift_mode="off")
        self.assertEqual("off", res["shift_mode"])
        self.assertEqual(0.0, res["shift"])

    def test_manual_shift(self):
        peaks = self.shifted_peaks()
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=45.0, shift=self.SHIFT)
        self.assertEqual("manual", res["shift_mode"])
        self.assertAlmostEqual(self.SHIFT, res["shift"], places=3)
        sc = self.by_id(res)["COD 1"]
        self.assertLess(sc["mean_dtheta"], 0.02)

    def test_estimate_shift_api(self):
        peaks = self.shifted_peaks()
        plain = phases.estimate_shift(peaks, self.db, tol=TOL)
        self.assertIsInstance(plain, tuple)
        self.assertEqual(2, len(plain))
        self.assertAlmostEqual(self.SHIFT, plain[0], delta=0.051)
        detail = phases.estimate_shift(peaks, self.db, tol=TOL, detail=True)
        self.assertIn("accepted", detail)
        self.assertTrue(detail["accepted"])

    def test_no_shift_on_aligned_data(self):
        spec = ALPHA_LINES + [ln for ln in BETA_LINES if ln[0] != 20.00] + NARROW2_LINES
        peaks = make_peaks(spec)
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=45.0)
        self.assertEqual(0.0, res["shift"])


# ------------------------------------------------------- 5. JARVIS vs COD
class TestJarvisVsCod(MiniDbCase):
    def setUp(self):
        self.peaks = make_peaks(ALPHA_LINES)
        self.res = phases.explain(self.peaks, self.db, tol=TOL, min_fom=45.0,
                                  shift_mode="off")

    def test_datasource_and_penalty(self):
        sc = self.by_id(self.res)
        self.assertIn("COD 1", sc)
        self.assertIn("COD 7", sc)
        self.assertEqual("COD", sc["COD 1"]["datasource"])
        self.assertEqual("JARVIS-DFT", sc["COD 7"]["datasource"])
        self.assertAlmostEqual(phases.DFT_PENALTY, sc["COD 1"]["fom"] - sc["COD 7"]["fom"],
                               places=6)
        self.assertLess(sc["COD 7"]["rank_score"], sc["COD 1"]["rank_score"])

    def test_db_filter_cod(self):
        res = phases.explain(self.peaks, self.db, tol=TOL, min_fom=45.0,
                             shift_mode="off", db_filter="cod")
        ids = self.ids(res)
        self.assertIn("COD 1", ids)
        self.assertNotIn("COD 7", ids)
        self.assertEqual("cod", res["settings"]["db_filter"])

    def test_db_filter_all_is_default(self):
        self.assertEqual("all", self.res["settings"]["db_filter"])
        self.assertEqual("all", self.res["db_meta"]["db_filter"])

    def test_invalid_filter(self):
        with self.assertRaises(ValueError):
            phases.explain(self.peaks, self.db, db_filter="pdf")


# --------------------------------------------------- 6. профиль стройматериалов
class TestProfile(MiniDbCase):
    def setUp(self):
        self.peaks = make_peaks(QUARTZ_LINES)

    def test_profile_records_exist_and_consistent(self):
        ids = [e["id"] for e in building_materials.PROFILE]
        self.assertEqual(len(ids), len(set(ids)), "дубли id в профиле")
        for entry in building_materials.PROFILE:
            self.assertTrue(entry["id"].startswith("COD "), entry)
            self.assertTrue(entry["formula"])
            self.assertIn(entry["weight"], (0.5, 1.0))
            self.assertTrue(entry["group"])

    def test_bonus_is_ranking_only(self):
        on = phases.explain(self.peaks, self.db, tol=TOL, min_fom=45.0,
                            shift_mode="off", use_profile=True)
        off = phases.explain(self.peaks, self.db, tol=TOL, min_fom=45.0,
                             shift_mode="off", use_profile=False)
        q_on = self.by_id(on)["COD 1011097"]
        q_off = self.by_id(off)["COD 1011097"]
        self.assertEqual(["building_materials"], q_on["prior"])
        self.assertEqual([], q_off["prior"])
        self.assertAlmostEqual(q_on["fom"] + phases.PRIOR_BONUS, q_on["rank_score"], places=3)
        self.assertAlmostEqual(q_off["fom"], q_off["rank_score"], places=3)
        # приор ничего не фильтрует: набор фаз тот же
        self.assertEqual(set(self.ids(on)), set(self.ids(off)))
        # и не скрывает альтернативы
        self.assertIn("COD 8", self.ids(on))

    def test_report_marks_prior_and_label(self):
        on = phases.explain(self.peaks, self.db, tol=TOL, min_fom=45.0,
                            shift_mode="off", use_profile=True)
        txt = phases.format_report(on, self.peaks)
        self.assertIn("prior: building_materials", txt)
        self.assertIn("Кварц", txt)

    def test_profile_helpers(self):
        quartz = {"id": "COD 1011097", "formula": "SiO2"}
        unknown = {"id": "COD 999", "formula": "Zz"}
        self.assertEqual(1.0, building_materials.prior_weight(quartz))
        self.assertEqual(["building_materials"], building_materials.prior_tags(quartz))
        self.assertEqual("Кварц", building_materials.label(quartz))
        self.assertEqual(0.0, building_materials.prior_weight(unknown))
        self.assertEqual([], building_materials.prior_tags(unknown))
        self.assertIsNone(building_materials.label(unknown))


# ------------------------------------------------------- контракт результата
class TestResultContract(MiniDbCase):
    def setUp(self):
        self.peaks = make_peaks([(10.00, 2.0)] + ALPHA_LINES + [(48.42, 60.0)])
        self.res = phases.explain(self.peaks, self.db, tol=TOL, min_fom=40.0,
                                  min_share=0.25)

    def test_legacy_keys(self):
        for key in ("phases", "per_peak", "unassigned", "elements", "n_obs", "shift"):
            self.assertIn(key, self.res)

    def test_new_keys(self):
        for key in ("ambiguous", "peak_status", "n_peaks", "n_below_threshold",
                    "n_outside_range", "shift_mode", "shift_gain", "radiation",
                    "settings", "db_meta", "peak_density", "p_chance"):
            self.assertIn(key, self.res)
        for key in ("tol", "min_fom", "min_obs_int", "use_profile", "db_filter",
                    "shift_mode", "two_theta_min", "two_theta_max"):
            self.assertIn(key, self.res["settings"])
        self.assertIn("wavelength_A", self.res["radiation"])
        self.assertEqual("Cu", self.res["radiation"]["anode"])

    def test_phase_entry_fields(self):
        sc = self.res["phases"][0]
        for key in ("phase", "fom", "fom_raw", "cov_ref", "cov_obs", "key_hits",
                    "key_total", "n_matched", "n_ref", "n_ref_total",
                    "n_ref_out_of_range", "mean_dtheta", "matched", "explained_peaks",
                    "share", "status", "datasource", "prior", "flags", "rank_score"):
            self.assertIn(key, sc)
        self.assertIn(sc["status"], ("candidate", "ambiguous", "weak"))
        self.assertLessEqual(sc["n_matched"], sc["n_ref"])
        self.assertLessEqual(sc["n_ref"], sc["n_ref_total"])
        self.assertEqual(sc["n_ref"] + sc["n_ref_out_of_range"], sc["n_ref_total"])

    def test_per_peak_entry_fields(self):
        for entries in self.res["per_peak"].values():
            self.assertIsInstance(entries, list)
            for entry in entries:
                self.assertEqual({"phase_no", "dtheta", "ref_int", "overlap"},
                                 set(entry))
                self.assertIsInstance(entry["overlap"], bool)
                self.assertGreaterEqual(entry["phase_no"], 1)
                self.assertLessEqual(entry["phase_no"], len(self.res["phases"]))

    def test_peak_status_partition(self):
        statuses = self.res["peak_status"]
        self.assertEqual(set(range(len(self.peaks))), set(statuses))
        counts = {}
        for value in statuses.values():
            counts[value] = counts.get(value, 0) + 1
        self.assertEqual(len(self.peaks), sum(counts.values()))
        self.assertEqual(counts.get("below_threshold", 0), self.res["n_below_threshold"])
        self.assertEqual(counts.get("unexplained", 0), len(self.res["unassigned"]))
        for pidx in self.res["ambiguous"]:
            self.assertEqual("ambiguous", statuses[pidx])

    def test_report_contains_required_warnings(self):
        txt = phases.format_report(self.res, self.peaks, filename="x.raw")
        self.assertIn("FOM — оценка совпадения, НЕ вероятность", txt)
        self.assertIn("Сдвиг 2θ", txt)
        self.assertIn("Диапазон съёмки", txt)
        self.assertIn("ниже порога", txt)
        self.assertIn("статус:", txt)
        self.assertIn("ПЕРЕКРЫТИЯ", txt)


# ------------------------------------------------------------- излучение
class TestRadiation(MiniDbCase):
    def test_ref_positions_recomputed_by_d(self):
        co_lam = radiation.wavelength("Co", "ka_avg")
        d_quartz = 3.343                      # сильнейшая линия кварца
        cu_tt = radiation.two_theta_from_d(d_quartz, LAM)
        co_tt = radiation.two_theta_from_d(d_quartz, co_lam)
        self.assertGreater(abs(co_tt - cu_tt), 3.0)      # шкалы сильно разные

        db = make_db(self.tmpdir, [make_phase("COD 9", "Rr2", [(cu_tt, 100)],
                                             name="RadTest")], name="rad_db.json.gz")
        peaks = make_peaks([(co_tt, 100)])
        res_cu = phases.explain(peaks, db, tol=TOL, min_fom=30.0, shift_mode="off")
        self.assertEqual([], res_cu["phases"])
        res_co = phases.explain(peaks, db, tol=TOL, min_fom=30.0, shift_mode="off",
                                wavelength=co_lam)
        self.assertEqual(["COD 9"], self.ids(res_co))
        self.assertAlmostEqual(co_lam, res_co["radiation"]["wavelength_A"], places=5)

    def test_radiation_object_accepted(self):
        rad = radiation.Radiation.from_settings({"anode": "Cobalt"})
        self.assertEqual("Co", rad.anode)
        peaks = make_peaks(ALPHA_LINES)
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=40.0, shift_mode="off",
                             radiation=rad)
        self.assertEqual("Co", res["radiation"]["anode"])
        self.assertAlmostEqual(rad.lam, res["radiation"]["wavelength_A"], places=5)


# --------------------------------------------------------- search_phases API
class TestSearchPhases(MiniDbCase):
    def test_returns_full_space_indices(self):
        peaks = make_peaks([(10.00, 2.0)] + ALPHA_LINES)
        res = phases.search_phases(peaks, self.db, tol=TOL, min_fom=30.0,
                                   shift_mode="off")
        self.assertTrue(res)
        for sc in res:
            for _ri, pidx, _dth in sc["matched"]:
                self.assertGreaterEqual(pidx, 1)      # индекс 0 — слабый пик
                self.assertLess(pidx, len(peaks))
            self.assertIn("status", sc)
            self.assertIn("datasource", sc)

    def test_legacy_positional_args(self):
        peaks = make_peaks(ALPHA_LINES)
        res = phases.search_phases(peaks, self.db, TOL, 30.0, 5, None, 0.0)
        self.assertTrue(res)
        self.assertIn(res[0]["phase"]["formula"], ("Aa2Bb", "Aa2BbJ"))

    def test_invalid_shift_mode(self):
        with self.assertRaises(ValueError):
            phases.explain(make_peaks(ALPHA_LINES), self.db, shift_mode="magic")


# ----------------------------------------------------------- полная база (slow)
@unittest.skipUnless(os.environ.get("XRD_SLOW"), "XRD_SLOW не задан: тест на полной базе пропущен")
class TestFullDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = phases.load_db()
        cls.by_id = {p["id"]: p for p in cls.db["phases"]}

    def lines(self, pid, imin=5.0, tmax=56.0):
        return [(p[0], p[2]) for p in self.by_id[pid]["peaks"]
                if p[2] >= imin and p[0] <= tmax]

    def test_profile_ids_exist_in_real_db(self):
        found, missing = building_materials.verify(self.db)
        self.assertEqual([], missing)
        self.assertEqual(len(building_materials.PROFILE), len(found))

    def test_three_phase_mixture_found(self):
        spec = (self.lines("COD 1011097") + self.lines("COD 4502443")
                + self.lines("COD 9000112"))
        peaks = make_peaks(spec)
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=55.0,
                             two_theta_min=4.0, two_theta_max=56.0)
        ids = [sc["phase"]["id"] for sc in res["phases"]]
        for pid in ("COD 4502443", "COD 1011097", "COD 9000112"):
            self.assertIn(pid, ids[:6], "фаза {} не найдена в топ-6: {}".format(pid, ids))
        statuses = {sc["phase"]["id"]: sc["status"] for sc in res["phases"]}
        self.assertEqual("candidate", statuses["COD 4502443"])

    def test_no_phases_on_noise_like_spectrum(self):
        # 120 «пиков» через 0.43° — совпадения с эталонами должны быть
        # распознаны как случайные (chance-level), а не как фазы
        spec = [(4.5 + 0.43 * i, 100.0 - (i % 7) * 6.0) for i in range(120)]
        peaks = make_peaks(spec)
        res = phases.explain(peaks, self.db, tol=TOL, min_fom=55.0,
                             two_theta_min=4.0, two_theta_max=56.0)
        self.assertEqual([], res["phases"])
        self.assertGreater(res["p_chance"], 0.9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
