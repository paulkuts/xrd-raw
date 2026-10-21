#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Регрессионные тесты научной валидации xrd-raw (исполнитель D).

Тесты гоняют РЕАЛЬНЫЙ конвейер (winxrd2xy.parse -> peakpick.pick ->
phases.explain) на синтетических фикстурах из tests/fixtures/ и сверяют
результат с независимой истиной из литературных d-значений (JSON рядом
с каждым .raw).

Принципы:
  * тесты не подгоняются под текущее поведение: если код нарушает
    научное требование контракта (docs/CONTRACT.md), тест помечен
    @unittest.expectedFailure с объяснением и числами;
  * отсутствие в коде ещё не реализованного поля контракта (flags,
    status, n_ref_out_of_range, per_peak-словари) — это НЕ провал науки,
    поэтому такие проверки пропускаются через skipUnless и включаются
    автоматически, когда поле появится (ветка A);
  * expectedFailure, который вдруг стал проходить, unittest показывает как
    "unexpected success" — это сигнал, что дефект исправлен и декоратор
    пора снять (см. docs/SCIENTIFIC_VALIDATION.md).

Запуск (из корня репозитория):
    python -m unittest tests.test_regression -v
    XRD_SLOW=1 python -m unittest tests.test_regression -v   # + сдвиг по полной базе
"""

import copy
import json
import os
import re
import shutil
import sys
import tempfile
import time
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TOOLS = ROOT / "tools"
FIX = HERE / "fixtures"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(FIX))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp"))

for _stream in (sys.stdout, sys.stderr):     # кириллица в cp866/cp1251 консоли
    try:
        _stream.reconfigure(errors="replace")
    except Exception:                        # noqa: BLE001
        pass

import peakpick  # noqa: E402
import phases  # noqa: E402
import winxrd2xy  # noqa: E402

TOL = 0.25          # допуск по 2θ, как TOL_DEFAULT в phases.py
MIN_FOM = 55.0      # порог «кандидата» по умолчанию

FIXTURES = [
    "quartz_pure", "mix_qcg", "overlap_calcite_gypsum", "noise_slope",
    "shift_plus010", "narrow_20_35", "wrongwl_co_as_cu", "noise_only",
    "dft_weak",
]

_CACHE = {}


# ------------------------------------------------------------------ утилиты
def fixture(name):
    """Разобрать фикстуру и прогнать конвейер один раз на процесс."""
    if name in _CACHE:
        return _CACHE[name]
    truth = json.loads((FIX / (name + ".json")).read_text(encoding="utf-8"))
    data, step, meta = winxrd2xy.parse(str(FIX / (name + ".raw")))
    x = [a for a, _ in data]
    y = [i for _, i in data]
    t0 = time.time()
    peaks = peakpick.pick(x, y)
    t1 = time.time()
    db = phases.load_db()
    res = phases.explain(peaks, db, tol=TOL, min_fom=MIN_FOM)
    t2 = time.time()
    out = {"name": name, "truth": truth, "data": data, "step": step,
           "meta": meta, "peaks": peaks, "res": res, "db": db,
           "t_pick": t1 - t0, "t_explain": t2 - t1}
    _CACHE[name] = out
    return out


def truth_lines(truth_phase, min_int=0.0):
    """Позиции 2θ истинных линий фазы в диапазоне скана."""
    if "lines_in_range" in truth_phase:
        lines = truth_phase["lines_in_range"]
        pts = [(l["two_theta"], l["I_percent"]) for l in lines]
    else:
        pts = [(l["two_theta"], l["I_percent"])
               for l in truth_phase.get("all_lines_d_I", [])]
    # dft_weak хранит линии как пары (d, I) без 2θ
    if not pts and "all_lines_d_I" in truth_phase:
        pts = [(a, b) for a, b in truth_phase["all_lines_d_I"]]
    return [tt for tt, ii in pts if ii >= min_int]


def formula_elements(formula):
    """Химические элементы из формулы: CaSO4*2H2O -> {Ca, S, O, H}."""
    return set(re.findall(r"[A-Z][a-z]?", formula or ""))


def same_phase(cand_phase, truth_phase, tol=TOL, min_frac=0.6, min_int=10.0):
    """Совпадает ли запись базы с истинной фазой.

    Два независимых условия:
      1) набор отражений — не менее 60 % истинных линий (I >= 10 %) имеют
         рефлекс базы в пределах tol. Сравнение по d/2θ, а не по названию:
         в базе гипс называется и CaSO4*2H2O, и CaH4SO6, кварц —
         "Quartz low"/"SiO2";
      2) химия — набор элементов совпадает. Без этого условия любая
         случайная многолинейная структура «покрывает» 2-3 линии кварца:
         на фикстуре wrongwl Co->Cu запись DyGe3 давала 2 попадания из 3.
    """
    el_t = formula_elements(truth_phase.get("formula", ""))
    el_c = set(cand_phase.get("elements") or [])
    if el_t and el_c and el_t != el_c:
        return False
    lines = truth_lines(truth_phase, min_int=min_int)
    if len(lines) < 2:
        lines = truth_lines(truth_phase)
    if not lines:
        return False
    refs = [p[0] for p in cand_phase.get("peaks", [])]
    hits = sum(1 for tt in lines if any(abs(tt - r) <= tol for r in refs))
    return hits / float(len(lines)) >= min_frac


def found(name):
    """-> {ключ истинной фазы: (ранг 1..N, FOM)} для найденных истинных фаз."""
    fx = fixture(name)
    out = {}
    for key, tp in enumerate(fx["truth"]["truth"]["present"]):
        tkey = tp.get("key") or tp.get("name") or str(key)
        for rank, sc in enumerate(fx["res"]["phases"], 1):
            if same_phase(sc["phase"], tp):
                out[tkey] = (rank, sc["fom"])
                break
    return out


def foreign(name):
    """Фазы результата, не совпавшие ни с одной истинной фазой."""
    fx = fixture(name)
    bad = []
    for rank, sc in enumerate(fx["res"]["phases"], 1):
        if not any(same_phase(sc["phase"], tp)
                   for tp in fx["truth"]["truth"]["present"]):
            bad.append((rank, (sc["phase"].get("name") or "")[:34],
                        sc["phase"].get("formula", ""), round(sc["fom"], 1),
                        sc["phase"].get("src", "")))
    return bad


def peak_rows(name):
    """Пики с интенсивностью >= MIN_OBS_INT (участвуют в поиске фаз)."""
    fx = fixture(name)
    return [p for p in fx["peaks"] if p["intensity"] >= phases.MIN_OBS_INT]


def unexplained(name):
    """Пики, которые конвейер не объяснил (по per_peak результата)."""
    fx = fixture(name)
    peaks = fx["peaks"]
    res = fx["res"]
    per = res.get("per_peak") or {}
    return [i for i, p in enumerate(peaks)
            if p["intensity"] >= phases.MIN_OBS_INT and i not in per]


def describe(name):
    """Диагностическая строка: что именно вернул конвейер."""
    fx = fixture(name)
    res = fx["res"]
    rows = ["{}: пиков {}, учтено {}, не объяснено {}".format(
        name, len(fx["peaks"]), res.get("n_obs"), len(unexplained(name)))]
    for rank, sc in enumerate(res.get("phases", []), 1):
        ph = sc["phase"]
        rows.append("  #{} {:<34s} {:<14s} FOM {:5.1f} (raw {:5.1f}) "
                    "линий {}/{} key {}/{} src {}".format(
                        rank, (ph.get("name") or "")[:34],
                        (ph.get("formula") or "")[:14], sc["fom"],
                        sc.get("fom_raw", float("nan")), sc.get("n_matched", -1),
                        sc.get("n_ref", -1), sc.get("key_hits", -1),
                        sc.get("key_total", -1), ph.get("src", "")))
    return "\n".join(rows)


def small_db(formulas=("SiO2", "CaCO3", "CaH4SO6", "Ca(OH)2")):
    """Мини-база для быстрых проверок сдвига (полная база — 65 с на estimate_shift).

    Это урезанное ПРОСТРАНСТВО ПОИСКА, а не источник истины: истина фикстур
    задана литературными d и лежит в JSON.
    """
    db = phases.load_db()
    keep = []
    for p in db["phases"]:
        if p.get("prio") != 1:
            continue
        if p.get("formula") not in formulas:
            continue
        name = (p.get("name") or "")
        if p.get("formula") == "SiO2" and name != "Quartz low":
            continue
        keep.append(p)
    return {"meta": db.get("meta", {}), "phases": keep}


def caps(name):
    """Какие поля расширенного контракта реально есть в результате."""
    res = fixture(name)["res"]
    any_sc = res["phases"][0] if res["phases"] else {}
    return {
        "flags": "flags" in any_sc,
        "status": "status" in any_sc,
        "n_ref_out_of_range": "n_ref_out_of_range" in any_sc,
        "peak_status": "peak_status" in res,
        "shift_mode": "shift_mode" in res,
        "per_peak_dict": bool(res.get("per_peak")) and isinstance(
            list(res["per_peak"].values())[0], list) and bool(
            res["per_peak"]) and isinstance(
            list(res["per_peak"].values())[0][0], dict),
    }


# --------------------------------------------------------- целостность фикстур
class TestFixtureIntegrity(unittest.TestCase):
    """Фикстуры должны быть маленькими, детерминированными и без потерь данных."""

    def test_all_fixtures_listed_and_present(self):
        for name in FIXTURES:
            self.assertTrue((FIX / (name + ".raw")).is_file(), name)
            self.assertTrue((FIX / (name + ".json")).is_file(), name)

    def test_size_budget(self):
        for name in FIXTURES:
            size = (FIX / (name + ".raw")).stat().st_size
            self.assertLess(size, 200 * 1024,
                            "{}: {} байт > 200 КБ".format(name, size))

    def test_parser_reads_truth_geometry(self):
        for name in FIXTURES:
            fx = fixture(name)
            tt = fx["truth"]["two_theta"]
            self.assertEqual(len(fx["data"]), tt["points"], name)
            self.assertAlmostEqual(fx["data"][0][0], tt["start"], places=6)
            self.assertAlmostEqual(fx["data"][-1][0], tt["stop"], places=6)
            self.assertAlmostEqual(fx["step"], tt["step"], places=6)
            self.assertEqual(fx["meta"].get("anode"), "Copper", name)

    def test_no_data_loss_against_generator(self):
        """Числа в .raw должны совпасть с тем, что строит генератор (шаг/точки/уровни)."""
        import make_fixtures

        for name in ("quartz_pure", "narrow_20_35"):
            fx = fixture(name)
            spec = make_fixtures.SCENARIOS[name][0]()
            self.assertEqual(len(fx["data"]), len(spec["y"]))
            diffs = [abs(a - b) for (_, a), b in zip(fx["data"], spec["y"])]
            self.assertEqual(max(diffs), 0, "{}: данные искажены".format(name))

    def test_generator_is_deterministic(self):
        import make_fixtures

        tmp = Path(tempfile.mkdtemp(prefix="xrd_fix_"))
        try:
            for name in ("quartz_pure", "noise_only"):
                raw, _doc = make_fixtures.build_one(name, tmp)
                committed = (FIX / (name + ".raw")).read_bytes()
                self.assertEqual(raw.read_bytes(), committed,
                                 "{}: файл не воспроизводится побитово".format(name))
        finally:
            shutil.rmtree(str(tmp), ignore_errors=True)

    def test_truth_json_is_self_describing(self):
        for name in FIXTURES:
            doc = fixture(name)["truth"]
            self.assertEqual(doc["kind"], "synthetic")
            self.assertIn("MIT", doc["license"])
            self.assertIn("generator", doc)
            self.assertIn("seed", doc)
            self.assertIn("wavelength_A", doc)
            self.assertIn("truth", doc)


# --------------------------------------------------------------- одна фаза
class TestSinglePhase(unittest.TestCase):
    def test_quartz_is_first_and_clean(self):
        res = fixture("quartz_pure")["res"]
        self.assertTrue(res["phases"], "кварц не найден:\n" + describe("quartz_pure"))
        fx = fixture("quartz_pure")
        self.assertTrue(same_phase(res["phases"][0]["phase"],
                                   fx["truth"]["truth"]["present"][0]),
                        "первая фаза не кварц:\n" + describe("quartz_pure"))
        self.assertGreaterEqual(res["phases"][0]["fom"], 70.0)

    @unittest.expectedFailure
    def test_quartz_result_has_no_foreign_phases(self):
        """Кроме кварца в результате не должно быть ничего.

        НЕ РЕШЕНО (научная проблема, а не дефект скоринга): BeF2 (α-кварцевый
        тип) и LiSmAlF6 дают набор линий, совпадающий с кварцем в пределах
        допуска 0.25° — по позициям PXRD они неразличимы, отличаются только
        относительные интенсивности (у BeF2 среднее расхождение в 2.1 раза,
        FOM 59.7). Пока нет ни химического ввода (какие элементы есть в
        образце), ни профиля/текстуры, отсечь изоструктурный вариант нельзя,
        не отсекая заодно настоящие фазы в смесях. Сейчас такие карточки
        остаются кандидатами с пометкой NOT_SIGNIFICANT; правильно — понижать
        их до «альтернативная интерпретация» (см. docs/VALIDATION.md).
        """
        self.assertEqual(foreign("quartz_pure"), [],
                         "ложные фазы:\n" + describe("quartz_pure"))

    def test_quartz_strongest_line_is_right(self):
        """Контрольная линия кварца 3.343 Å (26.66° 2θ) должна быть среди пиков."""
        fx = fixture("quartz_pure")
        best = max(fx["peaks"], key=lambda p: p["intensity"])
        self.assertAlmostEqual(best["two_theta"], 26.66, delta=0.05)
        self.assertAlmostEqual(best["d"], 3.343, delta=0.006)


# --------------------------------------------------------------- смеси
class TestMixture(unittest.TestCase):
    def test_all_three_phases_found(self):
        """Смесь кварц+кальцит+гипс: все три фазы обязаны быть в результате.

        ТЕКУЩЕЕ СОСТОЯНИЕ (baseline): не найдена НИ ОДНА. Первая итерация
        жадного explain() выбирает случайную запись базы (V2CuO6, JARVIS-DFT,
        FOM 64.0), которая «объясняет» 7 пиков, включая сильнейшие 26.67°
        (кварц) и 29.43° (кальцит). После этого настоящие фазы уже не набирают
        min_share=0.25, потому что их пики заняты. Контракт §5 требует
        совместной привязки без монополии одной фазы на пик.
        """
        got = found("mix_qcg")
        missing = [tp.get("key") for tp in
                   fixture("mix_qcg")["truth"]["truth"]["present"]
                   if (tp.get("key") or tp.get("name")) not in got]
        self.assertEqual(missing, [],
                         "не найдены фазы {}:\n{}".format(missing,
                                                          describe("mix_qcg")))

    @unittest.expectedFailure
    def test_mixture_has_no_junk_above_true_phases(self):
        """Ложные кандидаты не должны стоять выше настоящих фаз."""
        fx = fixture("mix_qcg")
        true_ranks = [r for r, _ in found("mix_qcg").values()]
        bad = foreign("mix_qcg")
        self.assertTrue(bad, "ожидались ложные кандидаты — проверка потеряла смысл")
        best_true = min(true_ranks) if true_ranks else 99
        self.assertEqual([b for b in bad if b[0] < best_true], [],
                         "ложные фазы выше настоящих:\n" + describe("mix_qcg"))


# --------------------------------------------------------------- перекрытие
class TestOverlap(unittest.TestCase):
    """Кальцит 3.035 Å (29.40°) и гипс 3.06 Å (29.17°) — общий допуск 0.25°."""

    def test_both_phases_survive_overlap(self):
        got = found("overlap_calcite_gypsum")
        self.assertIn("calcite", got, describe("overlap_calcite_gypsum"))
        self.assertIn("gypsum", got, describe("overlap_calcite_gypsum"))

    def test_true_phases_rank_above_junk(self):
        got = found("overlap_calcite_gypsum")
        worst_true = max(r for r, _ in got.values())
        self.assertTrue(worst_true <= 2,
                        "настоящая фаза ниже 2-го места:\n" +
                        describe("overlap_calcite_gypsum"))

    def test_overlapping_peak_listed_for_both_phases(self):
        if not caps("overlap_calcite_gypsum")["per_peak_dict"]:
            self.skipTest("per_peak ещё не в формате контракта §4 (ветка A)")
        fx = fixture("overlap_calcite_gypsum")
        res = fx["res"]
        # пики в области перекрытия 29.0-29.6°
        idx = [i for i, p in enumerate(fx["peaks"])
               if 29.0 <= p["two_theta"] <= 29.6]
        self.assertTrue(idx)
        claimants = set()
        for i in idx:
            for ent in res["per_peak"].get(i, []):
                claimants.add(ent["phase_no"])
        self.assertGreaterEqual(len(claimants), 2,
                                "общий пик приписан только одной фазе: "
                                "{}".format({round(fx['peaks'][i]['two_theta'], 2):
                                             res['per_peak'].get(i) for i in idx}))


# --------------------------------------------------------- шум и фон
class TestNoiseAndBackground(unittest.TestCase):
    def test_no_phase_on_pure_noise(self):
        res = fixture("noise_only")["res"]
        self.assertEqual(res["phases"], [],
                         "фазы найдены на чистом шуме:\n" + describe("noise_only"))

    def test_false_peaks_on_noise_are_few(self):
        """На чистом шуме (σ=12, фон 190–211 имп.) пиков быть почти не должно.

        ТЕКУЩЕЕ СОСТОЯНИЕ: peakpick.pick() возвращает 120 пиков — ровно
        max_peaks (упор в потолок). Порог prominence привязан к максимуму
        сигнала (2 % от максимума), а не к шуму (MAD σ), как требует
        контракт §6. Любой порог «от максимума» на шумовом спектре даёт
        десятки ложных пиков.
        """
        n = len(fixture("noise_only")["peaks"])
        self.assertLessEqual(n, 10, "ложных пиков {} (порог 10)".format(n))

    def test_weak_signal_on_steep_background(self):
        """Кварц 260 имп. на фоне 1160->60 имп. должен быть найден.

        Проверено после снижения порога пик-пикинга до 4.3σ (шум по MAD):
        кварц находится, ложных пиков на наклонном фоне 3–4.
        """
        got = found("noise_slope")
        self.assertIn("quartz", got,
                      "кварц не найден на наклонном фоне:\n" +
                      describe("noise_slope"))

    def test_background_does_not_create_extra_peaks(self):
        """Число пиков при слабом сигнале не должно втрое превышать число линий."""
        fx = fixture("noise_slope")
        n_lines = len([l for l in
                       fx["truth"]["truth"]["present"][0]["lines_in_range"]
                       if l["I_percent"] >= 5.0])
        n_peaks = len(fx["peaks"])
        self.assertLessEqual(n_peaks, 2 * n_lines,
                             "пиков {} при {} линиях с I>=5%".format(n_peaks, n_lines))


# --------------------------------------------------------------- сдвиг 2θ
class TestShift(unittest.TestCase):
    """Сдвиг проверяется на мини-базе: estimate_shift по полной базе ~65 с."""

    @classmethod
    def setUpClass(cls):
        cls.db = small_db()

    def test_shift_detected(self):
        fx = fixture("shift_plus010")
        sh, gain = phases.estimate_shift(fx["peaks"], self.db)
        self.assertAlmostEqual(sh, 0.10, delta=0.10,
                               msg="сдвиг +0.10° не обнаружен (получено {})".format(sh))
        self.assertGreater(sh, 0.0, "знак сдвига неверен")
        self.assertGreater(gain, 0.0)

    def test_shift_not_invented(self):
        fx = fixture("quartz_pure")
        sh, _gain = phases.estimate_shift(fx["peaks"], self.db)
        self.assertLess(abs(sh), 0.03,
                        "сдвиг выдуман на несдвинутом образце: {}".format(sh))

    def test_shift_improves_identification(self):
        fx = fixture("shift_plus010")
        sh, _ = phases.estimate_shift(fx["peaks"], self.db)
        res = phases.explain(fx["peaks"], self.db, tol=TOL, min_fom=MIN_FOM,
                             shift=sh)
        keys = set()
        for tp in fx["truth"]["truth"]["present"]:
            for sc in res["phases"]:
                if same_phase(sc["phase"], tp):
                    keys.add(tp["key"])
        self.assertEqual(keys, {"quartz", "calcite"},
                         "после поправки сдвига фазы не найдены:\n{}".format(
                             [(round(s["fom"], 1),
                               (s["phase"].get("name") or "")[:24])
                              for s in res["phases"]]))

    @unittest.skipUnless(os.environ.get("XRD_SLOW"),
                         "полная база: estimate_shift ~65 с (XRD_SLOW=1)")
    def test_shift_on_full_database(self):
        fx = fixture("shift_plus010")
        sh, _gain = phases.estimate_shift(fx["peaks"], phases.load_db())
        self.assertAlmostEqual(sh, 0.10, delta=0.10)


# --------------------------------------------------------- узкий диапазон
class TestNarrowRange(unittest.TestCase):
    """20–35°: сильные линии гипса (7.63, 4.28 Å) и портландита (4.90 Å) вне скана."""

    @unittest.expectedFailure
    def test_gypsum_found_despite_out_of_range_lines(self):
        """Отсутствие линий вне диапазона не должно убивать фазу.

        ТЕКУЩЕЕ СОСТОЯНИЕ: гипс в результате отсутствует; выбираются
        Sr6Mg23/YbFeO3. Скоринг берёт 12 сильнейших рефлексов эталона ВСЕГДА,
        включая 7.63 Å (11.64°) и 4.28 Å (20.7°, за нижней границей), и
        штрафует за них как за пропущенные ключевые линии.
        """
        got = found("narrow_20_35")
        self.assertIn("gypsum", got,
                      "гипс не найден:\n" + describe("narrow_20_35"))

    def test_single_line_phase_is_not_confirmed(self):
        """Портландит в 20–35° имеет одну линию (2.628 Å) — уверенного кандидата быть не должно."""
        fx = fixture("narrow_20_35")
        tp = [t for t in fx["truth"]["truth"]["present"]
              if t["key"] == "portlandite"][0]
        self.assertEqual(len(tp["lines_in_range"]), 1)
        res = fx["res"]
        cand = None
        for rank, sc in enumerate(res["phases"], 1):
            if same_phase(sc["phase"], tp):
                cand = (rank, sc)
                break
        if cand is None:
            return          # фаза не выдана — это и есть корректное поведение
        rank, sc = cand
        # если фаза всё же выдана, она обязана быть помечена как «мало данных»
        flags = sc.get("flags") or []
        status = sc.get("status")
        self.assertTrue(status in ("weak", "insufficient") or
                        "FEW_DIAGNOSTIC_LINES" in flags or
                        (status not in ("candidate",) and rank > 1),
                        "портландит по одной линии подан как уверенный кандидат "
                        "(ранг {}, FOM {:.1f}, status={}, flags={})".format(
                            rank, sc["fom"], status, flags))

    def test_out_of_range_lines_do_not_penalize_score(self):
        """Штраф за линии вне диапазона скана: контракт §5.3 запрещает.

        ТЕКУЩЕЕ СОСТОЯНИЕ: FOM гипса на фикстуре 20–35° = 55.7, а если у той же
        записи базы удалить рефлексы вне диапазона — 83.9. Разница 28.2 балла
        возникает только из-за линий, которые физически не могли наблюдаться.
        """
        fx = fixture("narrow_20_35")
        pk = peak_rows("narrow_20_35")
        gyp = None
        for p in phases.load_db()["phases"]:
            if p.get("formula") == "CaH4SO6" and p.get("prio") == 1:
                gyp = p
                break
        self.assertIsNotNone(gyp, "в базе нет гипса CaH4SO6")
        lo = fx["truth"]["two_theta"]["start"]
        hi = fx["truth"]["two_theta"]["stop"]
        trunc = copy.deepcopy(gyp)
        trunc["peaks"] = [r for r in gyp["peaks"] if lo - TOL <= r[0] <= hi + TOL]
        a = phases._score_phase(gyp, pk, TOL)
        b = phases._score_phase(trunc, pk, TOL)
        self.assertIsNotNone(a)
        self.assertIsNotNone(b)
        self.assertLess(abs(a["fom"] - b["fom"]), 5.0,
                        "линии вне диапазона штрафуют FOM на {:.1f} балла "
                        "(полный {:.1f} vs только в диапазоне {:.1f})".format(
                            b["fom"] - a["fom"], a["fom"], b["fom"]))


# ------------------------------------------------------- неверная длина волны
class TestWavelength(unittest.TestCase):
    """Снято в Co Kα, а метаданные прибора говорят Copper."""

    def test_wrong_wavelength_breaks_identification(self):
        """Ожидаемое поведение: истинные фазы не находятся (и это надо сообщать).

        Тест фиксирует ОТСУТСТВИЕ ложной уверенности: при неверной длине волны
        конвейер не должен выдавать кварц/кальцит как найденные. Если однажды
        появится пересчёт по излучению (radiation.py), тест начнёт падать —
        тогда его нужно переписать на «фазы найдены после пересчёта по d».
        """
        got = found("wrongwl_co_as_cu")
        self.assertEqual(got, {},
                         "при неверной длине волны «найдены» истинные фазы: "
                         "{}\n{}".format(got, describe("wrongwl_co_as_cu")))

    def test_wrong_wavelength_result_is_not_confident(self):
        """И даже ложные кандидаты не должны выглядеть как уверенное определение."""
        fx = fixture("wrongwl_co_as_cu")
        top_fom = fx["res"]["phases"][0]["fom"] if fx["res"]["phases"] else 0.0
        # Верхний порог «уверенности» — 85 баллов FOM у лучшего кандидата.
        self.assertLess(top_fom, 85.0,
                        "ложный кандидат с FOM {:.1f}:\n{}".format(
                            top_fom, describe("wrongwl_co_as_cu")))

    def test_metadata_wavelength_conflict_is_visible(self):
        """Фикстура обязана нести конфликт метаданных и фактического излучения."""
        doc = fixture("wrongwl_co_as_cu")["truth"]
        self.assertEqual(doc["anode_metadata"], "Copper")
        self.assertEqual(doc["anode_actual"], "Co")
        self.assertAlmostEqual(doc["wavelength_A"], 1.78901, places=5)


# ------------------------------------------------------------- DFT-кандидат
class TestDFTCandidate(unittest.TestCase):
    """Слабый набор из 4 линий расчётной фазы (JARVIS-DFT NbN, JVASP-15025)."""

    def test_dft_phase_never_presented_as_experimental(self):
        fx = fixture("dft_weak")
        for sc in fx["res"]["phases"]:
            src = sc["phase"].get("src", "")
            if src != "COD":
                self.assertLessEqual(sc["fom"], sc.get("fom_raw", sc["fom"]) - 7.9,
                                     "DFT-фаза без штрафа: {}".format(sc["phase"]))

    @unittest.expectedFailure
    def test_weak_four_line_pattern_has_no_confident_candidate(self):
        """4 линии (27/11/100/27 %) — недостаточно для уверенного кандидата.

        ТЕКУЩЕЕ СОСТОЯНИЕ: выдаются два кандидата COD (K3BiSe3 FOM 76.4,
        Bi2SeO2 67.3), не имеющих отношения к истинной фазе; статуса
        "weak"/флагов FEW_DIAGNOSTIC_LINES в структуре результата нет.
        """
        fx = fixture("dft_weak")
        res = fx["res"]
        if not res["phases"]:
            return
        statuses = set(sc.get("status") for sc in res["phases"])
        has_status = caps("dft_weak")["status"]
        if has_status:
            self.assertTrue(statuses <= {"weak"},
                            "уверенные кандидаты на 4 линиях: {}".format(statuses))
        else:
            self.assertEqual(res["phases"], [],
                             "уверенные кандидаты без статуса:\n" +
                             describe("dft_weak"))


# ------------------------------------------------- политика скоринга (общая)
class TestScoringPolicy(unittest.TestCase):
    def test_dft_penalty_applied_everywhere(self):
        for name in FIXTURES:
            for sc in fixture(name)["res"]["phases"]:
                if sc["phase"].get("src") == "JARVIS-DFT":
                    self.assertLessEqual(
                        sc["fom"], sc.get("fom_raw", sc["fom"]) - 7.9,
                        "{}: DFT-фаза без штрафа".format(name))

    def test_fom_is_bounded(self):
        for name in FIXTURES:
            for sc in fixture(name)["res"]["phases"]:
                self.assertGreaterEqual(sc["fom"], 0.0)
                self.assertLessEqual(sc["fom"], 100.0 + 1e-9)

    def test_unassigned_indices_point_into_full_peak_list(self):
        """Индексы результата — в пространстве полного списка пиков (контракт §3)."""
        for name in FIXTURES:
            fx = fixture(name)
            n = len(fx["peaks"])
            for i in fx["res"].get("unassigned", []):
                self.assertTrue(0 <= i < n, "{}: индекс {} вне списка".format(name, i))
            for i in fx["res"].get("per_peak", {}):
                self.assertTrue(0 <= i < n, "{}: per_peak {} вне списка".format(name, i))

    def test_report_renders_without_crash(self):
        for name in FIXTURES:
            fx = fixture(name)
            txt = phases.format_report(fx["res"], fx["peaks"], filename=name)
            self.assertIn("ФАЗЫ", txt)
            self.assertIn("ПИКИ", txt)


# ---------------------------------------------- пробелы контракта (ветка A)
class TestContractFields(unittest.TestCase):
    """Проверки полей, которых в текущем коде ещё нет (их делает ветка A).

    Пропускаются, пока поля отсутствуют; включаются автоматически, когда
    структура результата станет соответствовать docs/CONTRACT.md §4.
    """

    def test_peaks_have_peak_id_and_snr(self):
        peaks = fixture("quartz_pure")["peaks"]
        self.assertTrue(peaks)
        missing = [k for k in ("peak_id", "snr", "wavelength")
                   if k not in peaks[0]]
        self.assertEqual(missing, [], "в пиках нет полей {}".format(missing))

    def test_explain_exposes_shift_mode(self):
        res = fixture("quartz_pure")["res"]
        self.assertIn("shift_mode", res)

    def test_peak_status_covers_all_peaks(self):
        res = fixture("quartz_pure")["res"]
        self.assertIn("peak_status", res)
        self.assertEqual(len(res["peak_status"]),
                         len(fixture("quartz_pure")["peaks"]))

    def test_n_ref_limited_to_scan_range(self):
        """n_ref — только рефлексы в диапазоне скана (контракт §5.3)."""
        res = fixture("narrow_20_35")["res"]
        if not res["phases"]:
            self.skipTest("на фикстуре нет фаз — проверять нечего")
        for sc in res["phases"]:
            if sc["phase"].get("formula") != "CaH4SO6":
                continue
            self.assertIn("n_ref_out_of_range", sc)
            self.assertGreater(sc["n_ref_out_of_range"], 0)
            self.assertLessEqual(sc["n_ref"], 8)

    def test_few_diagnostic_lines_flag(self):
        """Фаза с одной линией в диапазоне обязана получить флаг."""
        res = fixture("narrow_20_35")["res"]
        flagged = [sc for sc in res["phases"]
                   if "FEW_DIAGNOSTIC_LINES" in (sc.get("flags") or [])]
        self.assertTrue(flagged, "нет ни одного флага FEW_DIAGNOSTIC_LINES")


def _install_contract_skips():
    """Пропустить проверки полей, которых ещё нет в текущем коде."""
    c = caps("quartz_pure")
    gates = {
        "test_peaks_have_peak_id_and_snr": None,
        "test_explain_exposes_shift_mode": c["shift_mode"],
        "test_peak_status_covers_all_peaks": c["peak_status"],
        "test_n_ref_limited_to_scan_range": c["n_ref_out_of_range"],
        "test_few_diagnostic_lines_flag": c["flags"],
    }
    for name, ok in gates.items():
        fn = getattr(TestContractFields, name)
        reason = ("поле контракта ещё не реализовано (ветка A); "
                  "проверка включится автоматически")
        if ok is None:
            peaks = fixture("quartz_pure")["peaks"]
            ok = bool(peaks) and all(k in peaks[0] for k in
                                     ("peak_id", "snr", "wavelength"))
        if not ok:
            setattr(TestContractFields, name, unittest.skip(reason)(fn))


_install_contract_skips()


if __name__ == "__main__":
    unittest.main(verbosity=2)
