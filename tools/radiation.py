#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Единый источник параметров излучения для всего проекта.

Раньше длина волны была зашита в разных модулях по-разному:
peakpick.py брал Cu Kα 1.5418, winxrd2xy.py и xlsx_out.py — Cu Kα1 1.5406,
а эталонная база посчитана на Cu Kα 1.5418 (pymatgen XRDCalculator("CuKa")).
Из-за этого d в CSV и d в пиках расходились, а сравнение с эталонами шло
по разным шкалам.

Теперь длина волны задаётся здесь один раз и передаётся во все модули.
Эталонные линии базы хранятся как [2θ, d, I%]; при длине волны, отличной от
базовой, позиции эталонов пересчитываются ПО d (закон Вульфа-Брэгга), а не
сравниваются как 2θ другого излучения.

Никаких внешних зависимостей — только стандартный Python.
"""

import math

# --------------------------------------------------------------- длины волн, Å
# ka1, ka2 — дублет Kα; ka_avg — средневзвешенная (веса 2:1), kbeta — Kβ.
ANODES = {
    "Cu": {"ka1": 1.54056, "ka2": 1.54439, "ka_avg": 1.54184, "kbeta": 1.39222},
    "Co": {"ka1": 1.78897, "ka2": 1.79285, "ka_avg": 1.78901, "kbeta": 1.62079},
    "Mo": {"ka1": 0.70930, "ka2": 0.71359, "ka_avg": 0.71073, "kbeta": 0.63229},
    "Cr": {"ka1": 2.28970, "ka2": 2.29361, "ka_avg": 2.29092, "kbeta": 2.08487},
    "Fe": {"ka1": 1.93604, "ka2": 1.93998, "ka_avg": 1.93728, "kbeta": 1.75661},
}

# Метаданные прибора (WinXRD пишет материал анода словом)
ANODE_ALIASES = {
    "copper": "Cu", "cu": "Cu",
    "cobalt": "Co", "co": "Co",
    "molybdenum": "Mo", "mo": "Mo",
    "chromium": "Cr", "cr": "Cr",
    "iron": "Fe", "fe": "Fe",
}

# Длина волны, в которой посчитана эталонная база phases_db.json.gz
# (build_db2.py: WAVELENGTH = 1.5418; pymatgen XRDCalculator(wavelength="CuKa")).
BASE_WAVELENGTH = 1.54184
BASE_LABEL = "Cu Kα средневзвешенная, 1.5418 Å"

MODELS = {
    "ka_avg": "Kα1+Kα2 средневзвешенная (2:1)",
    "ka1": "Kα1 (монохроматор / только α1)",
}

DEFAULT_ANODE = "Cu"
DEFAULT_MODEL = "ka_avg"


def xray_number_mode():
    """Режим округления не нужен, но имя оставлено для совместимости вызовов."""
    return 0


def wavelength(anode=DEFAULT_ANODE, model=DEFAULT_MODEL):
    """Длина волны (Å) для анода и спектральной модели."""
    a = ANODES.get(anode) or ANODES[DEFAULT_ANODE]
    key = model if model in ("ka1", "ka_avg") else DEFAULT_MODEL
    return float(a[key])


def anode_from_meta(meta):
    """Анод из метаданных прибора: 'Cu' или None, если метаданные не дают ответа."""
    if not meta:
        return None
    raw = str(meta.get("anode") or "").strip().lower()
    if not raw:
        return None
    return ANODE_ALIASES.get(raw)


def d_from_two_theta(two_theta, lam=None):
    """Межплоскостное расстояние d (Å) по 2θ (град)."""
    lam = BASE_WAVELENGTH if lam is None else lam
    th = math.radians(two_theta / 2.0)
    s = math.sin(th)
    if s <= 1e-9:
        return float("inf")
    return lam / (2.0 * s)


def two_theta_from_d(d, lam=None):
    """Угол 2θ (град) по d (Å); nan, если отражение недостижимо."""
    lam = BASE_WAVELENGTH if lam is None else lam
    if d is None or d <= 0:
        return float("nan")
    x = lam / (2.0 * d)
    if x >= 1.0:
        return float("nan")
    return 2.0 * math.degrees(math.asin(x))


class Radiation(object):
    """Параметры излучения для одного запуска.

    Поля:
      anode   — 'Cu' | 'Co' | 'Mo' | 'Cr' | 'Fe'
      model   — 'ka_avg' | 'ka1'
      lam     — длина волны, Å
      source  — 'metadata' | 'user' | 'default'
      confirmed — True, если анод и модель заданы пользователем/метаданными явно
      notes   — список строк-пояснений для отчёта
    """

    def __init__(self, anode=DEFAULT_ANODE, model=DEFAULT_MODEL, source="default",
                 confirmed=False, notes=None):
        self.anode = anode
        self.model = model if model in MODELS else DEFAULT_MODEL
        self.lam = wavelength(self.anode, self.model)
        self.source = source
        self.confirmed = bool(confirmed)
        self.notes = list(notes or [])

    # ------------------------------------------------------------------ сборка
    @classmethod
    def from_settings(cls, meta=None, anode=None, model=None):
        """Собрать конфигурацию из метаданных прибора и явных настроек пользователя.

        Приоритет: явная настройка пользователя > метаданные прибора > Cu по умолчанию.
        Отсутствие надёжных метаданных не маскируется: в notes попадает
        предупреждение, а confirmed остаётся False.
        """
        notes = []
        meta_anode = anode_from_meta(meta)
        src = "default"
        conf = False
        if anode:
            a = ANODE_ALIASES.get(str(anode).strip().lower(), anode if anode in ANODES else None)
            if a not in ANODES:
                notes.append("неизвестный анод '{}' — принято Cu (предположение)".format(anode))
                a = DEFAULT_ANODE
                src = "default"
            else:
                src = "user"
                conf = True
        elif meta_anode:
            a = meta_anode
            src = "metadata"
            conf = True
        else:
            a = DEFAULT_ANODE
            notes.append("анод не определён по метаданным — принято Cu (ПРЕДПОЛОЖЕНИЕ, не подтверждено)")
        if model:
            m = str(model).strip()
            if m not in MODELS:
                notes.append("неизвестная модель '{}' — принята средневзвешенная Kα".format(model))
                m = DEFAULT_MODEL
            else:
                conf = conf and True
        else:
            m = DEFAULT_MODEL
            notes.append("спектральная модель не задана — принята {} (по умолчанию)".format(MODELS[m]))
        return cls(a, m, source=src, confirmed=conf, notes=notes)

    # ------------------------------------------------------------------ вывод
    def describe(self, short=False):
        """Строка для отчёта/подписи графика."""
        label = "Kα1" if self.model == "ka1" else "Kα"
        src = {
            "user": "задано пользователем",
            "metadata": "из метаданных прибора",
            "default": "ПРЕДПОЛОЖЕНИЕ (метаданные не дали ответа)",
        }.get(self.source, self.source)
        if short:
            return "{} {} λ={:.4f} Å".format(self.anode, label, self.lam)
        return "{} {} ({}), λ = {:.4f} Å; источник: {}".format(
            self.anode, label, MODELS[self.model], self.lam, src)

    def matches_base(self, tol=5e-4):
        """Совпадает ли длина волны с той, в которой посчитана база."""
        return abs(self.lam - BASE_WAVELENGTH) <= tol

    def ref_two_theta(self, d, ref_two_theta=None):
        """Позиция эталонного рефлекса в текущем излучении.

        Пересчёт всегда по d — сравнивать 2θ разных излучений нельзя.
        Если d отсутствует (старые записи), используется ref_two_theta как есть.
        """
        if d:
            tt = two_theta_from_d(d, self.lam)
            if tt == tt:  # не nan
                return tt
        return ref_two_theta

    def d(self, two_theta):
        return d_from_two_theta(two_theta, self.lam)

    def as_dict(self):
        return {
            "anode": self.anode,
            "model": self.model,
            "wavelength_A": round(self.lam, 5),
            "source": self.source,
            "confirmed": self.confirmed,
            "notes": list(self.notes),
            "base_wavelength_A": BASE_WAVELENGTH,
        }


def default(meta=None):
    """Короткий путь: конфигурация по умолчанию с учётом метаданных прибора."""
    return Radiation.from_settings(meta)


if __name__ == "__main__":
    for a in ANODES:
        for m in ("ka_avg", "ka1"):
            r = Radiation(a, m)
            print("{:3s} {:7s} λ = {:.5f} Å | {}".format(a, m, r.lam, r.describe()))
    import sys
    sys.stdout.write("\n")
    demo = Radiation.from_settings({"anode": "Copper"})
    print("demo:", demo.describe())
    print("2θ кальцита d=3.035 Å:", round(float(demo.ref_two_theta(3.035)), 3))
    print("as_dict:", demo.as_dict())
