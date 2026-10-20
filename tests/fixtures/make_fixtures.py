#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Генератор синтетических фикстур для научной валидации xrd-raw.

Все дифрактограммы в tests/fixtures/*.raw — СИНТЕТИКА, созданная в этом
проекте (MIT). Это не измерения и не выдержки из картотеки PDF ICDD.

Источник истины — литературные значения межплоскостных расстояний d (Å) и
относительных интенсивностей I (%) для фаз, взятые из открытых справочников:
  * кварц SiO2        — PDF 46-1045 (общеизвестный набор);
  * кальцит CaCO3     — PDF 05-0586;
  * гипс CaSO4*2H2O   — PDF 33-0311;
  * портландит Ca(OH)2— PDF 44-1481;
  * эттрингит         — Handbook of Mineralogy (используется в docs как контроль).

Позиции пиков считаются из d по закону Вульфа—Брэга, а НЕ берутся из
tools/phases_db.json.gz. Единственное исключение — фикстура dft_weak:
для фаз, существующих только в расчётной базе JARVIS-DFT, независимого
литературного набора линий нет, поэтому её рефлексы взяты из базы и это
явно помечено в truth-JSON (line_source) и в README.

Профиль пика — псевдо-Фойгт (eta=0.5), FWHM растёт с углом, фон —
константа + наклон + (опционально) экспоненциальный спад, шум — гауссов
с фиксированным seed. Всё детерминировано: один и тот же seed даёт
побайтово тот же файл.

Использование:
    python tests/fixtures/make_fixtures.py            # перезаписать все .raw
    python tests/fixtures/make_fixtures.py --only quartz_pure
    python tests/fixtures/make_fixtures.py --xy       # ещё и .xy (не коммитятся)

Формат .raw — WinXRD (Thermo ARL X'TRA), см. build/make_test_raw.py:
текстовые строки метаданных, затем запись-заголовок (0, N, start_2theta)
и N записей по 16 байт: int32 intensity | int32 dwell | double two_theta.
"""

import argparse
import json
import math
import random
import struct
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent            # корень репозитория
DB_PATH = ROOT / "tools" / "phases_db.json.gz"

# Длина волны, в которой посчитана эталонная база (Cu Kα средневзвешенная).
CU_KA = 1.54184
CO_KA = 1.78901
MO_KA = 0.71073

FIXTURE_DATE = "2026-10-20"

# --------------------------------------------------------------- литература
# (d, Å; I, %). Порядок — по d. Значения округлены до точности справочника.
LITERATURE = {
    "quartz": {
        "name": "quartz (alpha-SiO2)",
        "formula": "SiO2",
        "reference": "PDF 46-1045 (literature d/I set)",
        "lines": [
            (4.257, 22.0), (3.343, 100.0), (2.457, 8.0), (2.282, 8.0),
            (2.237, 4.0), (2.127, 6.0), (1.980, 4.0), (1.818, 14.0),
            (1.802, 4.0), (1.671, 5.0), (1.659, 4.0), (1.608, 2.0),
            (1.541, 4.0), (1.453, 3.0), (1.418, 3.0), (1.375, 6.0),
        ],
    },
    "calcite": {
        "name": "calcite (CaCO3)",
        "formula": "CaCO3",
        "reference": "PDF 05-0586 (literature d/I set)",
        "lines": [
            (3.035, 100.0), (2.495, 14.0), (2.285, 18.0), (2.095, 17.0),
            (1.927, 8.0), (1.875, 17.0), (1.604, 8.0), (1.525, 6.0),
            (1.510, 4.0), (1.473, 4.0), (1.440, 4.0), (1.422, 3.0),
            (1.355, 3.0), (1.297, 3.0),
        ],
    },
    "gypsum": {
        "name": "gypsum (CaSO4*2H2O)",
        "formula": "CaSO4*2H2O",
        "reference": "PDF 33-0311 (literature d/I set)",
        "lines": [
            (7.63, 100.0), (4.28, 60.0), (3.80, 25.0), (3.06, 90.0),
            (2.87, 25.0), (2.69, 15.0), (2.60, 10.0), (2.50, 8.0),
            (2.44, 6.0), (2.32, 8.0), (2.22, 10.0), (2.13, 6.0),
            (2.07, 10.0), (1.90, 8.0),
        ],
    },
    "portlandite": {
        "name": "portlandite (Ca(OH)2)",
        "formula": "Ca(OH)2",
        "reference": "PDF 44-1481 (literature d/I set); strongest line 2.63 A",
        "lines": [
            (4.90, 74.0), (2.628, 100.0), (1.927, 41.0), (1.796, 30.0),
            (1.687, 15.0), (1.573, 8.0), (1.482, 6.0), (1.414, 5.0),
            (1.361, 4.0), (1.312, 3.0),
        ],
    },
    "ettringite": {
        "name": "ettringite",
        "formula": "Ca6Al2(SO4)3(OH)12*26H2O",
        "reference": "Handbook of Mineralogy (9.73/5.61/3.88/2.564/2.209 A)",
        "lines": [
            (9.73, 100.0), (5.61, 80.0), (4.69, 35.0), (3.88, 50.0),
            (2.773, 40.0), (2.564, 45.0), (2.209, 45.0),
        ],
    },
}


# ------------------------------------------------------------------ физика
def two_theta_from_d(d, lam=CU_KA):
    """2θ (град) по d (Å) и длине волны (Å)."""
    x = lam / (2.0 * d)
    if x >= 1.0:
        return None
    return 2.0 * math.degrees(math.asin(x))


def pseudo_voigt(dx, fwhm, eta=0.5):
    """Псевдо-Фойгт с максимумом 1 в центре; dx — отступ от центра, град."""
    if fwhm <= 0:
        return 0.0
    g = math.exp(-4.0 * math.log(2.0) * (dx / fwhm) ** 2)
    lor = 1.0 / (1.0 + 4.0 * (dx / fwhm) ** 2)
    return eta * lor + (1.0 - eta) * g


def fwhm_at(two_theta, fwhm0=0.090, fwhm_slope=0.0012):
    """Ширина пика прибора: линейно растёт с углом (0.11° при 20°, 0.16° при 60°)."""
    return fwhm0 + fwhm_slope * two_theta


def build_pattern(components, start, stop, step, lam=CU_KA, shift=0.0,
                  noise_sigma=6.0, bg_offset=30.0, bg_slope=0.0,
                  bg_exp=None, seed=1, eta=0.5):
    """Собрать синтетическую дифрактограмму.

    components — список dict: {"lines": [(d, I%)], "scale": максимум сильнейшей
    линии в импульсах, "name": ...}. Возвращает (x, y, info).
    """
    n = int(round((stop - start) / step)) + 1
    x = [start + k * step for k in range(n)]
    y = [0.0] * n

    # фон: константа + наклон + экспоненциальный спад (аморфное рассеяние)
    for k in range(n):
        v = bg_offset + bg_slope * (x[k] - start)
        if bg_exp:
            amp, decay = bg_exp
            v += amp * math.exp(-decay * (x[k] - start))
        y[k] += v

    placed = []
    for comp in components:
        lines = comp["lines"]
        imax = max(i for _, i in lines) or 1.0
        scale = float(comp.get("scale", 1000.0))
        for d, rel in lines:
            tt = two_theta_from_d(d, lam)
            if tt is None or not (start - 1.0 <= tt + shift <= stop + 1.0):
                continue
            pos = tt + shift
            height = scale * rel / imax
            fw = fwhm_at(pos)
            k0 = int(round((pos - start) / step))
            half = max(4, int(round(3.0 * fw / step)))
            for k in range(max(0, k0 - half), min(n, k0 + half + 1)):
                y[k] += height * pseudo_voigt(x[k] - pos, fw, eta=eta)
            placed.append({"phase": comp.get("name", ""), "d": d, "rel": rel,
                           "two_theta": round(pos, 3), "height": round(height, 2)})

    rnd = random.Random(seed)
    out = []
    for v in y:
        v = v + rnd.gauss(0.0, noise_sigma)
        out.append(int(round(max(0.0, v))))
    info = {"components": placed, "seed": seed}
    return x, out, info


# --------------------------------------------------------------- WinXRD .raw
def write_raw(path, start, step, intensities, dwell=100, anode="Copper",
              comment="synthetic"):
    """Записать .raw в формате WinXRD (16-байтных записей)."""
    head = ("{}\r\nThermo ARL X'TRA\r\nWinXRD 2.0-6\r\n{}\r\n"
            "Solid State Thermo\r\nair\r\nOperator: {}\r\n").format(
        path.stem if hasattr(path, "stem") else "synthetic",
        anode, comment)
    data = bytearray(head.encode("latin-1"))
    data += struct.pack("<iid", 0, len(intensities), start)
    for k, v in enumerate(intensities):
        data += struct.pack("<iid", int(v), dwell, start + k * step)
    Path(path).write_bytes(bytes(data))
    return Path(path)


def write_xy(path, x, y):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for a, v in zip(x, y):
            f.write("{:8.2f}{:12d}\n".format(a, v))
    return Path(path)


# ------------------------------------------------------------------ сценарии
def _tt(d, lam):
    """2θ по d; ошибка, если отражение недостижимо (не должно случаться в фикстурах)."""
    tt = two_theta_from_d(d, lam)
    if tt is None:
        raise ValueError("d={} Å недостижимо при λ={} Å".format(d, lam))
    return tt


def _phase(key, scale):
    lit = LITERATURE[key]
    return {"name": lit["name"], "formula": lit["formula"],
            "lines": lit["lines"], "scale": scale, "key": key}


def _truth_phase(key, scale, start, stop, lam, shift):
    lit = LITERATURE[key]
    lines = []
    for d, rel in lit["lines"]:
        tt = two_theta_from_d(d, lam)
        if tt is not None and start - 0.2 <= tt + shift <= stop + 0.2:
            lines.append({"d": d, "I_percent": rel,
                          "two_theta": round(tt + shift, 3)})
    return {"key": key, "name": lit["name"], "formula": lit["formula"],
            "reference": lit["reference"], "scale_counts": scale,
            "lines_in_range": lines}


def load_jarvis_lines(phase_id):
    """Рефлексы фазы JARVIS-DFT из базы (единственный случай использования базы).

    Возвращает (name, formula, src, [(d, I%), ...]) в пространстве Cu Kα.
    """
    import gzip

    with gzip.open(str(DB_PATH), "rt", encoding="utf-8") as f:
        db = json.load(f)
    for p in db["phases"]:
        if p.get("id") == phase_id:
            lines = [(round(r[1], 3), round(r[2], 1)) for r in p["peaks"]]
            return p.get("name", ""), p["formula"], p.get("src", ""), lines
    raise KeyError("нет фазы {} в базе".format(phase_id))


SCENARIOS = {}


def scenario(name, **kw):
    def deco(fn):
        SCENARIOS[name] = (fn, kw)
        return fn
    return deco


@scenario("quartz_pure", title="Одиночная фаза: кварц",
          expectation="кварц найден первым, ложных кандидатов выше порога нет")
def _quartz_pure():
    start, stop, step, lam, shift, seed = 5.0, 60.0, 0.02, CU_KA, 0.0, 20261020
    comps = [_phase("quartz", 1500.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=6.0, bg_offset=28.0, bg_slope=-0.10,
                               seed=seed)
    truth = {
        "present": [_truth_phase("quartz", 1500.0, start, stop, lam, shift)],
        "absent": ["calcite", "gypsum", "portlandite"],
        "expect_no_confident_candidate": False,
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 6.0},
                background={"offset": 28.0, "slope_per_deg": -0.10})


@scenario("mix_qcg", title="Смесь: кварц + кальцит + гипс",
          expectation="все три фазы в результате, кварц/кальцит/гипс в правильном порядке по FOM")
def _mix_qcg():
    start, stop, step, lam, shift, seed = 5.0, 60.0, 0.02, CU_KA, 0.0, 20261021
    comps = [_phase("quartz", 900.0), _phase("calcite", 1600.0),
             _phase("gypsum", 750.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=8.0, bg_offset=30.0, bg_slope=-0.08,
                               seed=seed)
    truth = {
        "present": [_truth_phase("calcite", 1600.0, start, stop, lam, shift),
                    _truth_phase("quartz", 900.0, start, stop, lam, shift),
                    _truth_phase("gypsum", 750.0, start, stop, lam, shift)],
        "absent": ["portlandite", "ettringite"],
        "expect_no_confident_candidate": False,
        "ranking_note": "кальцит сильнейший по I, но FOM зависит от покрытия; "
                        "проверяется состав множества и отсутствие мусора",
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 8.0},
                background={"offset": 30.0, "slope_per_deg": -0.08})


@scenario("overlap_calcite_gypsum",
          title="Перекрытие: кальцит 3.035 + гипс 3.06 (Δ2θ ≈ 0.23°)",
          expectation="обе фазы остаются в результате; общий пик отмечен у обеих")
def _overlap():
    start, stop, step, lam, shift, seed = 8.0, 60.0, 0.02, CU_KA, 0.0, 20261022
    comps = [_phase("calcite", 1600.0), _phase("gypsum", 1400.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=6.0, bg_offset=26.0, bg_slope=-0.05,
                               seed=seed)
    truth = {
        "present": [_truth_phase("calcite", 1600.0, start, stop, lam, shift),
                    _truth_phase("gypsum", 1400.0, start, stop, lam, shift)],
        "absent": ["quartz", "portlandite"],
        "expect_no_confident_candidate": False,
        "overlap_2theta": [round(_tt(3.035, lam), 2), round(_tt(3.06, lam), 2)],
        "note": "линии 3.035 (кальцит, 100%) и 3.06 Å (гипс, 90%) попадают "
                "в один допуск 0.25° — классическое перекрытие",
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 6.0},
                background={"offset": 26.0, "slope_per_deg": -0.05})


@scenario("noise_slope", title="Слабый сигнал на сильном наклонном фоне",
          expectation="кварц найден; ложных пиков мало (порог по шуму, не по максимуму)")
def _noise_slope():
    start, stop, step, lam, shift, seed = 5.0, 60.0, 0.02, CU_KA, 0.0, 20261023
    comps = [_phase("quartz", 260.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=10.0, bg_offset=60.0,
                               bg_exp=(1100.0, 0.055), seed=seed)
    truth = {
        "present": [_truth_phase("quartz", 260.0, start, stop, lam, shift)],
        "absent": ["calcite", "gypsum"],
        "expect_no_confident_candidate": False,
        "note": "фон 1160 -> 60 имп., сигнал кварца 260 имп. при 3.343 Å",
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 10.0},
                background={"offset": 60.0, "exp_amp": 1100.0, "exp_decay": 0.055})


@scenario("shift_plus010", title="Систематический сдвиг 2θ = +0.10°",
          expectation="сдвиг обнаруживается (~+0.10°), фазы остаются кварц+кальцит")
def _shift():
    start, stop, step, lam, shift, seed = 10.0, 60.0, 0.02, CU_KA, 0.10, 20261024
    comps = [_phase("quartz", 1200.0), _phase("calcite", 1500.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=6.0, bg_offset=28.0, bg_slope=-0.06,
                               seed=seed)
    truth = {
        "present": [_truth_phase("calcite", 1500.0, start, stop, lam, shift),
                    _truth_phase("quartz", 1200.0, start, stop, lam, shift)],
        "absent": ["gypsum", "portlandite"],
        "expect_no_confident_candidate": False,
        "shift_deg": shift,
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 6.0},
                background={"offset": 28.0, "slope_per_deg": -0.06})


@scenario("narrow_20_35", title="Узкий диапазон 20–35°: гипс + портландит",
          expectation="сильные линии вне диапазона не штрафуют; при 1 диагностической "
                      "линии фаза помечается «недостаточно данных», а не candidate")
def _narrow():
    start, stop, step, lam, shift, seed = 20.0, 35.0, 0.02, CU_KA, 0.0, 20261025
    comps = [_phase("gypsum", 1200.0), _phase("portlandite", 1200.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=5.0, bg_offset=25.0, bg_slope=0.0,
                               seed=seed)
    truth = {
        "present": [_truth_phase("gypsum", 1200.0, start, stop, lam, shift),
                    _truth_phase("portlandite", 1200.0, start, stop, lam, shift)],
        "absent": ["quartz", "calcite"],
        "expect_no_confident_candidate": False,
        "out_of_range_strong_lines": {
            "gypsum": [{"d": 7.63, "I_percent": 100.0},
                       {"d": 4.28, "I_percent": 60.0}],
            "portlandite": [{"d": 4.90, "I_percent": 74.0},
                            {"d": 1.927, "I_percent": 41.0}],
        },
        "note": "портландит имеет в диапазоне единственную линию 2.628 Å (100%)",
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 5.0},
                background={"offset": 25.0, "slope_per_deg": 0.0})


@scenario("wrongwl_co_as_cu",
          title="Неверная длина волны: снято в Co Kα, в метаданных «Copper»",
          expectation="идентификация ломается (ожидаемое поведение, фиксируется тестом)")
def _wrongwl():
    start, stop, step, lam, shift, seed = 8.0, 70.0, 0.02, CO_KA, 0.0, 20261026
    comps = [_phase("quartz", 1300.0), _phase("calcite", 1500.0)]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=6.0, bg_offset=28.0, bg_slope=-0.05,
                               seed=seed)
    truth = {
        "present": [_truth_phase("quartz", 1300.0, start, stop, lam, shift),
                    _truth_phase("calcite", 1500.0, start, stop, lam, shift)],
        "absent": [],
        "expect_no_confident_candidate": True,
        "wavelength_actual_A": CO_KA,
        "anode_actual": "Co",
        "anode_metadata": "Copper",
        "note": "позиции посчитаны для Co Kα (1.78901 Å), а метаданные прибора "
                "говорят Copper: конвейер обязан либо требовать подтверждения "
                "излучения, либо давать неверный/пустой результат — это и "
                "фиксируется как ожидаемое поведение",
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 6.0},
                background={"offset": 28.0, "slope_per_deg": -0.05})


@scenario("noise_only", title="Только шум + один посторонний пик",
          expectation="уверенных кандидатов нет; посторонний пик остаётся необъяснённым")
def _noise_only():
    start, stop, step, seed = 5.0, 60.0, 0.02, 20261027
    comps = []
    x, y, info = build_pattern(comps, start, stop, step, noise_sigma=12.0,
                               bg_offset=190.0, bg_slope=0.35, seed=seed)
    # один посторонний пик, которого нет в наборе фаз (d ~ 4.13 Å)
    rnd = random.Random(seed + 1)
    pos = 21.50
    fw = fwhm_at(pos)
    for k, xk in enumerate(x):
        y[k] = int(round(y[k] + 150.0 * pseudo_voigt(xk - pos, fw)))
    truth = {
        "present": [],
        "absent": ["quartz", "calcite", "gypsum", "portlandite"],
        "expect_no_confident_candidate": True,
        "foreign_peak": {"two_theta": pos, "d": round(CU_KA / (2 * math.sin(
            math.radians(pos / 2))), 3), "height_counts": 150},
        "note": "шум sigma=12 при фоне 190-211 имп.; один пик на 21.5°, "
                "не принадлежащий ни одной фазе набора",
    }
    return dict(start=start, stop=stop, step=step, lam=CU_KA, shift=0.0, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 12.0},
                background={"offset": 190.0, "slope_per_deg": 0.35})


@scenario("dft_weak", title="Слабый набор линий расчётной фазы (JARVIS-DFT NbN)",
          expectation="кандидат из DFT не подаётся как подтверждённый; штраф и пометка источника")
def _dft():
    name, formula, src, lines = load_jarvis_lines("JVASP-15025")
    start, stop, step, lam, shift, seed = 15.0, 60.0, 0.02, CU_KA, 0.0, 20261028
    comps = [{"name": name, "formula": formula, "lines": lines, "scale": 700.0}]
    x, y, info = build_pattern(comps, start, stop, step, lam=lam, shift=shift,
                               noise_sigma=8.0, bg_offset=35.0, bg_slope=-0.05,
                               seed=seed)
    truth = {
        "present": [{"key": "dft", "name": name, "formula": formula,
                     "line_source": "phases_db.json.gz::JARVIS-DFT {}".format(
                         "JVASP-15025"),
                     "note": "единственная фикстура, где линии взяты из базы: "
                             "для DFT-фаз нет независимого литературного набора",
                     "all_lines_d_I": lines}],
        "absent": ["quartz", "calcite", "gypsum"],
        "expect_no_confident_candidate": True,
        "datasource_expected": src,
    }
    return dict(start=start, stop=stop, step=step, lam=lam, shift=shift, seed=seed,
                x=x, y=y, info=info, anode="Copper", truth=truth,
                noise={"model": "gaussian_abs", "sigma": 8.0},
                background={"offset": 35.0, "slope_per_deg": -0.05})


# ------------------------------------------------------------------- запись
def build_one(name, outdir, want_xy=False):
    fn, _kw = SCENARIOS[name]
    spec = fn()
    raw = outdir / "{}.raw".format(name)
    write_raw(raw, spec["start"], spec["step"], spec["y"],
              anode=spec["anode"], comment="synthetic fixture (MIT)")
    if want_xy:
        write_xy(outdir / "{}.xy".format(name), spec["x"], spec["y"])

    doc = {
        "fixture": name,
        "file": raw.name,
        "kind": "synthetic",
        "license": "MIT — создано в этом проекте; не является измерением",
        "generator": "tests/fixtures/make_fixtures.py",
        "created": FIXTURE_DATE,
        "seed": spec["seed"],
        "two_theta": {"start": spec["start"], "stop": spec["stop"],
                      "step": spec["step"], "points": len(spec["y"])},
        "wavelength_A": spec["lam"],
        "anode_actual": "Co" if abs(spec["lam"] - CO_KA) < 1e-6 else "Cu",
        "anode_metadata": spec["anode"],
        "shift_deg": spec["shift"],
        "profile": {"shape": "pseudo-Voigt", "eta": 0.5,
                    "fwhm_deg": "0.090 + 0.0012 * 2theta"},
        "noise": spec["noise"],
        "background": spec["background"],
        "peaks_placed": spec["info"]["components"],
        "truth": spec["truth"],
        "bytes": raw.stat().st_size,
    }
    # состав: из truth (у dft_weak он другой формы)
    comp = []
    for p in spec["truth"]["present"]:
        comp.append({"name": p.get("name"), "formula": p.get("formula"),
                     "source": p.get("line_source", "literature d/I")})
    doc["composition"] = comp
    (outdir / "{}.json".format(name)).write_text(
        json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return raw, doc


def main(argv=None):
    ap = argparse.ArgumentParser(description="Собрать синтетические фикстуры .raw")
    ap.add_argument("--outdir", default=str(HERE))
    ap.add_argument("--only", default=None, help="имя одного сценария")
    ap.add_argument("--xy", action="store_true", help="писать ещё и .xy")
    a = ap.parse_args(argv)
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = [a.only] if a.only else list(SCENARIOS)
    for name in names:
        if name not in SCENARIOS:
            raise SystemExit("нет сценария {}. Есть: {}".format(
                name, ", ".join(SCENARIOS)))
        raw, doc = build_one(name, outdir, want_xy=a.xy)
        tt = doc["two_theta"]
        print("{:<22s} {:>6d} байт  {:>5.1f}-{:>5.1f}°  шаг {}  фаз {}  файл {}".format(
            name, doc["bytes"], tt["start"], tt["stop"], tt["step"],
            len(doc["truth"]["present"]), raw.name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
