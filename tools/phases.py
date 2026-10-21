#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Идентификация кристаллических фаз по пикам дифрактограммы.

Работает офлайн: база эталонных порошковых паттернов лежит рядом
(phases_db.json.gz). Логика — как в поиске по картотеке PDF: сравнение
d-позиций и относительных интенсивностей, отбор по сильнейшим линиям
(Hanawalt-подобный) и совместная привязка пиков к фазам.

Что важно знать про результат
-----------------------------
* ``FOM`` — оценка совпадения (0..100), **не вероятность** и не количество
  вещества. Количественный состав по PXRD в этой программе не определяется.
* Пики, у которых I < ``min_obs_int``, в поиске **не участвуют**: они получают
  статус ``below_threshold`` и не считаются «необъяснёнными».
* Один пик может принадлежать нескольким фазам (перекрытие рефлексов). Никакого
  взаимоисключающего назначения нет: ``per_peak[peak_index]`` — список всех
  претендентов с флагом ``overlap``, ``ambiguous`` — пики с >= 2 фазами.
* Все индексы в результате — в пространстве **полного переданного списка**
  ``peaks`` (внутренняя фильтрация идёт через карту ``obs_index -> peak_index``).

Определения счётчиков эталонных линий (важно для отчёта)
--------------------------------------------------------
* ``n_ref_total``        — сколько рефлексов вообще есть в карточке эталона;
* ``n_ref_out_of_range`` — сколько из них лежит ВНЕ измеренного диапазона
                           ``[2θmin - tol, 2θmax + tol]`` (после поправки сдвига);
* ``n_ref``              — сколько сильнейших рефлексов В ДИАПАЗОНЕ реально
                           участвовало в расчёте (<= ``max_ref_lines``).
  Покрытие ``cov_ref`` и штрафы считаются ТОЛЬКО по линиям в диапазоне: нельзя
  штрафовать фазу за линии, которые прибор физически не мог измерить.

Сдвиг нуля 2θ
--------------
``shift_mode`` = ``'off'`` | ``'manual'`` | ``'auto'`` (по умолчанию ``'auto'``).
Числовое ``shift=`` означает ``manual``. Автосдвиг ищется в узком окне
(<= ``SHIFT_SPAN``, шаг ``SHIFT_STEP``) и принимается только при выигрыше
>= ``SHIFT_GAIN``; это диагностика юстировки (смещение образца, нуль шкалы),
а не подгонка структуры. Применённая поправка всегда попадает в результат
и в отчёт.

Публичные сигнатуры (совместимость)
-----------------------------------
``load_db(path=None)``, ``db_info(db)``, ``estimate_shift(peaks, db, tol, span,
step, top_n)``, ``search_phases(peaks, db, tol, min_fom, top, elements_only,
shift)``, ``explain(peaks, db, tol, min_fom, max_phases, min_share, shift)``,
``format_report(res, peaks, filename, max_cand)`` — позиционные аргументы
сохранены, всё новое добавлено kwarg-ами с дефолтами.
"""

import gzip
import json
import math
from pathlib import Path

try:                                            # единый источник излучения
    import radiation as _radiation
except Exception:                               # pragma: no cover - страховка
    _radiation = None

try:                                            # профиль стройматериалов (приор)
    import building_materials as _profile
except Exception:                               # pragma: no cover - страховка
    _profile = None

TOL_DEFAULT = 0.25          # допуск по 2θ, градусы
MIN_OBS_INT = 3.0           # пики слабее (% от макс.) не участвуют в поиске фаз
MAX_PHASES = 8              # сколько фаз максимум объясняют образец
DFT_PENALTY = 8.0           # штраф фазам из DFT-базы (нет экспериментального эталона)
MAX_REF_LINES = 12          # сколько сильнейших рефлексов эталона берём в расчёт
STRONG_REF_INT = 40.0       # «сильная» линия эталона, % — за пропуск штрафуем
MISS_PENALTY_PER_PCT = 0.30  # штраф за каждую пропущенную сильную линию
MAX_MISS_PENALTY = 25.0     # потолок штрафа за пропущенные линии
SINGLE_LINE_CAP = 40.0      # потолок FOM, если совпала ровно одна линия
INT_LIMIT_LOG = 0.477       # log10(3): среднее расхождение относительных интенсивностей
                            # в 3 раза — противоречие картине, а не подтверждение
INT_FLOOR = 0.55            # минимум множителя интенсивностей (полное расхождение)
WEAK_ANCHOR = 0.15          # сильный необъяснённый пик (доля от макс.), для которого
                            # сообщаем фазу «недостаточно данных»
MAX_WEAK_REPORTED = 2       # сколько таких фаз максимум дописываем в результат
MAX_WEAK_CANDIDATES = 3     # сколько «слабых кандидатов» показывать, когда порог не пройден
WEAK_CANDIDATE_MIN_FOM = 35.0   # ниже этого не показываем даже как слабого кандидата
SAME_PATTERN_TOL = 0.30     # две карточки с совпадающими линиями — одна фаза
SAME_PATTERN_FRAC = 0.75
CHANCE_PENALTY_SCALE = 10.0  # штраф за log10(ожидаемого числа ложных фаз)
MAX_CHANCE_PENALTY = 20.0   # потолок этого штрафа
DIAG_REL_INT = 10.0         # линия диагностична при I >= 10 % от сильнейшей в диапазоне
DIAG_MIN_INT = 5.0          # ... и не слабее 5 %
MAX_PER_FORMULA = 2         # не больше двух карточек одной формулы в результате
DUP_FOM_TOL = 3.0           # почти идентичные карточки (|ΔFOM| <= 3) не дублируем
PRIOR_BONUS = 3.0           # мягкий бонус профилю стройматериалов (только ранжирование)
SHIFT_SPAN = 0.20           # окно автосдвига, градусы
SHIFT_STEP = 0.05
SHIFT_GAIN = 0.08           # минимальный выигрыш, чтобы принять автосдвиг (8 %)
SHIFT_MAX_CANDIDATES = 40000  # предохранитель числа фаз при оценке сдвига
SHIFT_TOP_PHASES = 100      # второй проход оценки сдвига: только лучшие фазы
SHIFT_MIN_COHERENCE = 0.75  # доля совпавшего веса в одном бине, чтобы доверять фазе
SHIFT_MIN_MATCH = 0.25      # минимум совпавшего веса (сильная линия ~1.0)
SHIFT_MIN_Z = 3.5           # значимость пика гистограммы над фоном (в сигмах)
SHIFT_MIN_VOTES = 1.0       # минимум «фаз-голосов» в лучшем бине
SHIFT_MAX_CHANCE = 0.6     # при такой доле случайных совпадений сдвиг не оцениваем
SHIFT_MIN_LINES = 3         # минимум сильных линий в диапазоне, чтобы фазе доверять
SHIFT_MIN_REF_INT = 5.0     # линия участвует в голосовании по сдвигу при I >= 5 %
BUCKET_STEP = 0.1           # шаг углового индекса
RI_LIMIT = 128              # кодирование (phase_index, ref_index) в одно число
CACHE = {}

_FLAGS_TEXT = {
    "FEW_DIAGNOSTIC_LINES": "недостаточно данных: меньше 2 диагностических линий подтверждено",
    "SINGLE_LINE_MATCH": "совпала одна линия — совпадение может быть случайным",
    "INTENSITY_CONFLICT": "относительные интенсивности противоречат эталону — скорее всего "
                          "другая фаза или сильная текстура",
    "NO_STRONG_LINE": "не совпала ни одна сильная линия эталона",
    "RANGE_LIMITED": "диапазон съёмки уже эталона — часть линий вне измерения",
    "LOW_SHARE": "малая доля объяснённой интенсивности образца",
    "CHANCE_LEVEL_MATCH": "совпадение на уровне случайных (высокая плотность пиков — "
                          "проверьте пик-пикинг и шум)",
    "NOT_SIGNIFICANT": "ожидаемо среди случайных фаз базы (мало линий, слабая статистика)",
}


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


# ------------------------------------------------------------------- излучение
class _SimpleRadiation(object):
    """Минимальная замена tools/radiation.py, если модуль недоступен.

    Нужна только для пересчёта позиции эталона по d (закон Вульфа-Брэгга).
    """

    BASE_WAVELENGTH = 1.54184

    def __init__(self, lam=None, source="default"):
        self.anode = "Cu"
        self.model = "ka_avg"
        self.lam = float(lam) if lam else self.BASE_WAVELENGTH
        self.source = source
        self.confirmed = False

    def ref_two_theta(self, d, ref_two_theta=None):
        if d:
            x = self.lam / (2.0 * d)
            if x < 1.0:
                return 2.0 * math.degrees(math.asin(x))
        return ref_two_theta

    def matches_base(self, tol=5e-4):
        return abs(self.lam - self.BASE_WAVELENGTH) <= tol

    def as_dict(self):
        return {"anode": self.anode, "model": self.model,
                "wavelength_A": round(self.lam, 5), "source": self.source,
                "confirmed": False, "base_wavelength_A": self.BASE_WAVELENGTH}


def _make_radiation(radiation=None, wavelength=None):
    """Собрать объект излучения: объект, длина волны или дефолт Cu Kα."""
    if radiation is not None:
        if hasattr(radiation, "ref_two_theta"):
            return radiation
        if isinstance(radiation, (int, float)):
            return _SimpleRadiation(radiation, source="wavelength")
        raise TypeError("radiation: ожидался объект radiation.Radiation или длина волны")
    if wavelength is not None:
        return _SimpleRadiation(wavelength, source="wavelength")
    if _radiation is not None:
        return _radiation.Radiation.from_settings(None)
    return _SimpleRadiation()


def _lam_key(rad):
    """None, если позиции эталонов можно брать из базы как есть."""
    if rad is None:
        return None
    try:
        if rad.matches_base():
            return None
        return round(float(rad.lam), 5)
    except Exception:                                  # pragma: no cover
        return None


_MAPPED_CACHE = {}


def _ref_positions(ph, rad):
    """Рефлексы эталона [(2θ_ref, I%), ...] в текущем излучении (пересчёт по d)."""
    peaks = ph.get("peaks") or ()
    lam_key = _lam_key(rad)
    if lam_key is None:
        return [(float(p[0]), float(p[2])) for p in peaks]
    key = (ph.get("id"), lam_key, len(peaks))
    cached = _MAPPED_CACHE.get(key)
    if cached is not None:
        return cached
    out = []
    for p in peaks:
        d = p[1] if len(p) > 1 else None
        tt = rad.ref_two_theta(d, p[0])
        if tt is None or tt != tt:                     # nan / нет данных
            continue
        out.append((float(tt), float(p[2])))
    if len(_MAPPED_CACHE) > 50000:
        _MAPPED_CACHE.clear()
    _MAPPED_CACHE[key] = out
    return out


# ------------------------------------------------------------- угловые индексы
def _bucket(two_theta, step=BUCKET_STEP):
    return int(math.floor(float(two_theta) / step))


def _angle_index(db, rad=None, step=BUCKET_STEP):
    """Индекс: бакет 2θ -> закодированные (phase_index, ref_index).

    Кэшируется внутри самого ``db`` (без ``id(dict)``), ключ включает шаг и
    длину волны: для не-Cu излучения позиции эталонов пересчитываются по d.
    """
    cache = db.setdefault("_cache", {})
    lam_key = _lam_key(rad)
    key = ("angle_index", step, lam_key)
    idx = cache.get(key)
    if idx is not None:
        return idx, step
    phases = db.get("phases", ())
    idx = {}
    for pi, ph in enumerate(phases):
        refs = ph.get("peaks") or ()
        if len(refs) >= RI_LIMIT:                      # страховка кодирования
            continue
        for ri, ref in enumerate(refs):
            tt = float(ref[0])
            if lam_key is not None:
                tt = rad.ref_two_theta(ref[1] if len(ref) > 1 else None, ref[0])
                if tt is None or tt != tt:
                    continue
            idx.setdefault(_bucket(tt, step), []).append(pi * RI_LIMIT + ri)
    cache[key] = idx
    return idx, step


def _obs_index(obs_peaks, step=BUCKET_STEP):
    """Индекс: бакет 2θ -> индексы наблюдаемых пиков (для быстрого поиска)."""
    idx = {}
    for oi, p in enumerate(obs_peaks):
        idx.setdefault(_bucket(p["two_theta"], step), []).append(oi)
    return idx


def _near(pos, tol, index, step=BUCKET_STEP):
    """Все значения из индекса в окне pos ± tol."""
    out = []
    for b in range(_bucket(pos - tol, step), _bucket(pos + tol, step) + 1):
        out.extend(index.get(b, ()))
    return out


def _candidate_phases(obs_peaks, index, step, tol, shift=0.0):
    """Фазы, у которых есть рефлекс рядом с наблюдаемым пиком.

    Сравнение идёт с позицией эталона, сдвинутой на ``shift``:
    ``obs ≈ ref + shift``  =>  ищем эталон около ``obs - shift``.
    """
    cand = {}
    shift = float(shift or 0.0)
    for oi, op in enumerate(obs_peaks):
        pos = float(op["two_theta"]) - shift
        for code in _near(pos, tol, index, step):
            pi, _ri = divmod(code, RI_LIMIT)
            cand.setdefault(pi, set()).add(oi)
    return cand


def _datasource(ph):
    """'COD' или 'JARVIS-DFT' — по полю src (fallback: prio)."""
    src = str(ph.get("src") or "").strip()
    if src:
        return "COD" if src.upper() == "COD" else "JARVIS-DFT"
    return "COD" if int(ph.get("prio", 1) or 1) == 1 else "JARVIS-DFT"


# --------------------------------------------------------------------- скоринг
def _scan_bounds(obs_peaks, two_theta_min, two_theta_max, tol):
    """Границы диапазона поиска по эталону: [min - tol, max + tol]."""
    tmin, tmax = two_theta_min, two_theta_max
    if tmin is None or tmax is None:
        if obs_peaks:
            tts = [float(p["two_theta"]) for p in obs_peaks]
            if tmin is None:
                tmin = min(tts)
            if tmax is None:
                tmax = max(tts)
    if tmin is None or tmax is None:
        return None, None
    tmin, tmax = float(tmin), float(tmax)
    if tmax < tmin:
        tmin, tmax = tmax, tmin
    return tmin - tol, tmax + tol


def _score_phase(ph, obs_peaks, tol=TOL_DEFAULT, allowed_obs=None,
                 unmatched_penalty=True, max_ref_lines=MAX_REF_LINES, shift=0.0,
                 two_theta_min=None, two_theta_max=None, radiation=None,
                 obs_index=None, n_tested=0):
    """Сопоставить рефлексы фазы с наблюдаемыми пиками.

    В расчёт берутся сильнейшие рефлексы, но ТОЛЬКО те, что попадают в
    измеренный диапазон ``[2θmin - tol, 2θmax + tol]`` (после поправки сдвига):
    линии вне диапазона прибор не измерял, и штрафовать за них нельзя.
    Внутри одной фазы пик используется не более одного раза, но разные фазы
    могут претендовать на один и тот же пик (перекрытие) — это разрешено.

    Возвращает dict метрик (``matched`` — в индексах ``obs_peaks``) или None,
    если фаза не совпала ни с одной линией / ни с одной из трёх главных.
    """
    refs_all = _ref_positions(ph, radiation)
    if not refs_all:
        return None
    lo, hi = _scan_bounds(obs_peaks, two_theta_min, two_theta_max, tol)
    shift = float(shift or 0.0)

    in_range = []
    for ri, (rtt, rint) in enumerate(refs_all):
        pos = rtt + shift
        if lo is not None and (pos < lo or pos > hi):
            continue
        in_range.append((pos, rint, ri))
    n_ref_total = len(refs_all)
    n_out = n_ref_total - len(in_range)
    if not in_range:
        return None
    # Слабые линии карточки (< DIAG_MIN_INT %) не диагностичны: в плотной базе они
    # «находятся» у чужого образца случайно, а их вес раздувает покрытие. Берём в
    # расчёт только диагностичные линии; если таких меньше двух — отступаем к
    # сильнейшим в диапазоне (фаза всё равно получит FEW_DIAGNOSTIC_LINES).
    diag = [r for r in in_range if r[1] >= DIAG_MIN_INT]
    if len(diag) < 2:
        diag = in_range
    n_diag = len(diag)
    diag.sort(key=lambda r: (-r[1], r[0]))
    refs = diag[:max(1, int(max_ref_lines))]

    if obs_index is None:
        obs_index = _obs_index(obs_peaks)
    allowed = None if allowed_obs is None else set(allowed_obs)

    used = set()
    matched = []                                    # (ref_index, obs_index, dtheta)
    for pos, _rint, ri in refs:
        best = None
        for oi in _near(pos, tol, obs_index):
            if oi in used or (allowed is not None and oi not in allowed):
                continue
            dth = abs(float(obs_peaks[oi]["two_theta"]) - pos)
            if dth <= tol and (best is None or dth < best[0]):
                best = (dth, oi)
        if best is not None:
            used.add(best[1])
            matched.append((ri, best[1], best[0]))
    if not matched:
        return None

    ref_int = {ri: rint for _pos, rint, ri in refs}
    matched_ri = {ri for ri, _oi, _d in matched}
    key_idx = [ri for _pos, _rint, ri in refs[:3]]
    key_hits = sum(1 for ri in key_idx if ri in matched_ri)
    if key_hits == 0:                               # ни одной из трёх главных линий
        return None

    sum_ref = sum(max(rint, 0.5) for _pos, rint, _ri in refs) or 1.0
    cov_ref = sum(max(ref_int[ri], 0.5) for ri in matched_ri) / sum_ref

    obs_pool = allowed if allowed is not None else range(len(obs_peaks))
    obs_total = sum(max(float(obs_peaks[oi]["intensity"]), 0.5) for oi in obs_pool) or 1.0
    cov_obs = sum(max(float(obs_peaks[oi]["intensity"]), 0.5)
                  for _ri, oi, _d in matched) / obs_total

    mean_dtheta = sum(d for _ri, _oi, d in matched) / len(matched)

    imax_in = max(rint for _pos, rint, _ri in in_range) or 1.0
    diag_level = max(DIAG_MIN_INT, DIAG_REL_INT / 100.0 * imax_in)
    diag_ids = {ri for _pos, rint, ri in in_range if rint >= diag_level}
    n_diag = len(diag_ids)
    # Подтверждено диагностических линий — именно это, а не «сколько их в
    # диапазоне вообще», говорит, на чём держится идентификация. Фаза, у которой
    # совпала одна сильная линия (или одна и две слабые), — «недостаточно данных».
    n_diag_hit = len(diag_ids & matched_ri)

    # ---- согласие относительных интенсивностей ------------------------------
    # Позиции совпали — это ещё половина дела: у фазы должна совпадать и
    # картина интенсивностей. Сравниваем отношения I_эталона/I_наблюдения,
    # нормируя на медиану (масштаб фазы в смеси неизвестен), и считаем среднее
    # |log10| отклонение. Расхождение в 3 раза и больше — уже не подтверждение.
    q = []
    imax_obs_all = max(max(float(p["intensity"]), 0.5) for p in obs_peaks) or 1.0
    for ri, oi, _d in matched:
        r_rel = max(ref_int[ri], 0.5) / imax_in
        o_rel = max(float(obs_peaks[oi]["intensity"]), 0.5) / imax_obs_all
        if r_rel > 0 and o_rel > 0:
            q.append(r_rel / o_rel)
    if len(q) >= 2:
        scale = sorted(q)[len(q) // 2] or 1.0
        int_dev = sum(abs(math.log10(max(v / scale, 1e-6))) for v in q) / len(q)
    else:
        int_dev = 0.0
    int_agr = max(0.0, 1.0 - int_dev / INT_LIMIT_LOG)

    # ---- поправка на случайные совпадения -----------------------------------
    # Если пиков много (шумовой спектр), любая эталонная линия «найдёт» пик в
    # пределах tol. Считаем ожидаемую долю покрытия при случайном совпадении:
    # p_chance = 2*tol*плотность_пиков (пиков на градус). Покрытие сверх
    # случайного уровня и есть доказательство фазы.
    span = (hi - lo) if (lo is not None and hi is not None) else 0.0
    n_pool = len(allowed) if allowed is not None else len(obs_peaks)
    density = (n_pool / span) if span > 0 else 0.0
    p_chance = min(1.0, 2.0 * tol * density)
    if p_chance >= 1.0:
        cov_excess = 0.0
    else:
        cov_excess = max(0.0, (cov_ref - p_chance) / (1.0 - p_chance))

    pool_int = [max(float(obs_peaks[oi]["intensity"]), 0.5) for oi in obs_pool]
    imax_obs = max(pool_int) if pool_int else 1.0
    strong_cut = 0.15 * imax_obs
    strong_obs = [oi for oi in obs_pool
                  if max(float(obs_peaks[oi]["intensity"]), 0.5) >= strong_cut]
    matched_obs = {oi for _ri, oi, _d in matched}
    n_strong_hit = sum(1 for oi in strong_obs if oi in matched_obs)
    obs_strong = (n_strong_hit / float(len(strong_obs))) if strong_obs else 0.0
    obs_term = min(1.0, obs_strong / 0.5)

    pos_score = max(0.0, 1.0 - mean_dtheta / tol) ** 2 if tol > 0 else 1.0
    # Покрытие по 2-3 линиям — не покрытие: 2 из 2 совпадений даёт любая редкая
    # карточка. Вес покрытия растёт с числом подтверждённых линий.
    cov_evidence = cov_excess * min(1.0, len(matched) / 4.0)
    fom_raw = (45.0 * cov_evidence
               + 15.0 * (key_hits / float(len(key_idx)))
               + 25.0 * pos_score
               + 15.0 * obs_term)
    # Интенсивности — второй, независимый признак фазы: если картина
    # относительных интенсивностей противоречит наблюдению, FOM не может
    # оставаться высоким, даже когда все позиции совпали.
    fom_raw *= (INT_FLOOR + (1.0 - INT_FLOOR) * int_agr)

    flags = []
    if n_diag_hit < 2:
        flags.append("FEW_DIAGNOSTIC_LINES")
    if int_dev > INT_LIMIT_LOG:
        flags.append("INTENSITY_CONFLICT")
    if len(matched) == 1:
        flags.append("SINGLE_LINE_MATCH")
    if not any(ref_int[ri] >= STRONG_REF_INT for ri in matched_ri):
        flags.append("NO_STRONG_LINE")
    if n_out and len(in_range) < 3:
        flags.append("RANGE_LIMITED")
    if cov_excess < 0.15 and p_chance >= 0.5:
        flags.append("CHANCE_LEVEL_MATCH")

    # ---- поправка на множественную проверку ---------------------------------
    # Мы проверяем десятки тысяч фаз: если у фазы в диапазоне всего 2-3 линии,
    # «полное совпадение» ожидаемо встретится у сотен чужих карточек. Считаем
    # ожидаемое число таких ложных фаз и штрафуем за него.
    expected_false = 0.0
    if n_tested and p_chance > 0:
        expected_false = float(n_tested) * (p_chance ** len(matched))
    log_expected = math.log10(max(expected_false, 1e-300))
    chance_pen = min(MAX_CHANCE_PENALTY, CHANCE_PENALTY_SCALE * max(0.0, log_expected))
    if expected_false >= 1.0:
        flags.append("NOT_SIGNIFICANT")

    fom = fom_raw - chance_pen
    if unmatched_penalty:
        missed = sum(MISS_PENALTY_PER_PCT * rint
                     for _pos, rint, ri in refs
                     if ri not in matched_ri and rint >= STRONG_REF_INT)
        fom -= min(MAX_MISS_PENALTY, missed)

    source = _datasource(ph)
    if source != "COD":
        fom -= DFT_PENALTY
    if len(matched) == 1:                           # никакой ложной уверенности
        fom = min(fom, SINGLE_LINE_CAP)

    return {
        "fom": round(max(0.0, min(100.0, fom)), 3),
        "fom_raw": round(max(0.0, min(100.0, fom_raw)), 3),
        "cov_ref": round(cov_ref, 4),
        "cov_obs": round(cov_obs, 4),
        "cov_excess": round(cov_excess, 4),
        "p_chance": round(p_chance, 4),
        "obs_strong": round(obs_strong, 4),
        "int_agr": round(int_agr, 4),
        "int_dev": round(int_dev, 4),
        "expected_false": round(expected_false, 6),
        "chance_penalty": round(chance_pen, 3),
        "key_hits": key_hits,
        "key_total": len(key_idx),
        "n_matched": len(matched),
        "n_ref": len(refs),
        "n_ref_total": n_ref_total,
        "n_ref_out_of_range": n_out,
        "n_diag": n_diag,
        "n_diag_hit": n_diag_hit,
        "mean_dtheta": round(mean_dtheta, 4),
        "matched": matched,
        "datasource": source,
        "flags": flags,
    }


# ------------------------------------------------------------------------ сдвиг
def _resolve_shift(shift, shift_mode):
    """Нормализовать запрос сдвига: ('off'|'manual'|'auto', значение)."""
    if shift_mode is not None:
        mode = str(shift_mode).strip().lower()
        val = float(shift) if isinstance(shift, (int, float)) else 0.0
    elif shift is None:
        mode, val = "auto", 0.0
    elif isinstance(shift, str):
        s = shift.strip().lower()
        if s == "auto":
            mode, val = "auto", 0.0
        elif s in ("off", "none", "no", "нет"):
            mode, val = "off", 0.0
        else:
            raise ValueError("shift: неизвестный режим '{}'".format(shift))
    else:
        mode, val = "manual", float(shift)
    if mode not in ("off", "manual", "auto"):
        raise ValueError("shift_mode: ожидалось 'off'|'manual'|'auto', получено '{}'".format(mode))
    if mode == "off":
        val = 0.0
    return mode, val


def estimate_shift(peaks, db, tol=TOL_DEFAULT, span=SHIFT_SPAN, step=SHIFT_STEP,
                   restrict=None,
                   top_n=3, min_obs_int=MIN_OBS_INT, two_theta_min=None,
                   two_theta_max=None, radiation=None, db_filter="all",
                   detail=False):
    """Подобрать систематический сдвиг нуля по 2θ (юстировка, плоскость образца).

    Метод — голосование по остаткам (по сути кросс-корреляция остатков, а не
    подгонка FOM): для каждой пары «сильная эталонная линия — наблюдаемый пик»
    считается остаток ``obs - ref``; остатки складываются в гистограмму с шагом
    ``step`` в окне ±``span``, вес = (I_эталона/100)·(I_пика/I_макс). Сдвиг берётся
    как центр тяжести лучшего бина, принимается только если он набирает
    >= SHIFT_GAIN (8 %) больше нулевого бина и заметно выше фона остальных
    бинов. Иначе применяется 0.0.

    Почему не по FOM: сумма FOM лучших фаз растёт и от случайных совпадений,
    и оптимум уезжает (проверено: на сдвиге +0.10° FOM-перебор давал −0.20°).
    Голосование устойчивее и в ~10 раз дешевле.

    Возвращает ``(shift, gain)`` как раньше; при ``detail=True`` — dict с полями
    ``shift``, ``gain``, ``accepted``, ``base`` (вес нулевого бина), ``best``,
    ``candidates``, ``tried``, ``mode``.
    """
    db_filter = _check_db_filter(db_filter)
    phases = db.get("phases", ())
    rad = _make_radiation(radiation)
    obs, _mapping = _prepare_obs(peaks, min_obs_int)
    index, idx_step = _angle_index(db, rad)

    def result(sh, gain, base, best, accepted, ncand, tried, reason=""):
        det = {"shift": round(sh, 3), "gain": round(gain, 4), "base": round(base, 4),
               "best": round(best, 4), "accepted": bool(accepted), "reason": reason,
               "candidates": ncand, "tried": tried, "mode": "auto"}
        if detail:
            return det
        return det["shift"], det["gain"] if accepted else 0.0

    if not obs or step <= 0 or span <= 0:
        return result(0.0, 0.0, 0.0, 0.0, False, 0, 0)

    tmin, tmax = two_theta_min, two_theta_max
    if tmin is None:
        tmin = min(float(p["two_theta"]) for p in obs)
    if tmax is None:
        tmax = max(float(p["two_theta"]) for p in obs)
    lo, hi = _scan_bounds(obs, tmin, tmax, tol)

    # В шумовом спектре (пиков много) сдвиг оценивать нельзя: любая карточка
    # «совпадает» при любом сдвиге, и оценка становится случайной.
    span_deg = (hi - lo) if (lo is not None and hi is not None) else 0.0
    density = (len(obs) / span_deg) if span_deg > 0 else 0.0
    if min(1.0, 2.0 * tol * density) > SHIFT_MAX_CHANCE:
        return result(0.0, 0.0, 0.0, 0.0, False, 0, 0, reason="dense_peaks")

    imax = max(float(p["intensity"]) for p in obs) or 1.0
    strong_cut = max(min_obs_int, 0.10 * imax)
    strong = [p for p in obs if float(p["intensity"]) >= strong_cut]
    cand = _candidate_phases(strong, index, idx_step, tol + span + step)
    if db_filter == "cod":
        cand = {pi: v for pi, v in cand.items()
                if str(phases[pi].get("src", "COD")).upper() == "COD"}
    if not cand:
        return result(0.0, 0.0, 0.0, 0.0, False, 0, 0)
    if restrict is not None:
        # второй проход: голосуют только фазы, прошедшие скоринг, — иначе
        # 32k карточек топят гистограмму в фоне случайных совпадений
        cand = {pi: cand[pi] for pi in cand if pi in restrict}
    if len(cand) > SHIFT_MAX_CANDIDATES:
        # Отбирать по числу «зацепленных» пиков нельзя: в плотном спектре
        # выигрывают карточки с максимумом линий, а реальные фазы (7-8 линий)
        # вылетают. Отсекаем только заведомо бедные карточки.
        keep = [pi for pi in cand
                if len(phases[pi].get("peaks") or ()) >= SHIFT_MIN_LINES]
        if len(keep) > SHIFT_MAX_CANDIDATES:
            keep = keep[:SHIFT_MAX_CANDIDATES]
        cand = {pi: cand[pi] for pi in keep}

    obs_index = _obs_index(obs)
    obs_int = [max(float(p["intensity"]), 0.5) for p in obs]
    obs_max = max(obs_int) or 1.0

    nbins = int(round(2.0 * span / step)) + 1
    votes = [0.0] * nbins
    wsum = [0.0] * nbins                       # сумма resid*weight — для центра бина
    zero_bin = nbins // 2

    # Голосуют только «когерентные» фазы: у которых сильные линии дружно
    # указывают на ОДИН и тот же сдвиг. Мусорные карточки, чьи линии случайно
    # цепляются за разные пики, размазывают голоса по всем бинам и в подсчёт
    # не попадают — иначе фон случайных совпадений топит реальный сдвиг.
    for pi in cand:
        ph = phases[pi]
        local = [0.0] * nbins
        local_sum = [0.0] * nbins
        n_lines = 0
        weights = []                           # веса сильных линий в диапазоне
        voted = 0.0                            # вес линий, нашедших партнёра
        for rtt, rint in _ref_positions(ph, rad):
            if rint < SHIFT_MIN_REF_INT or rtt < lo or rtt > hi:
                continue
            n_lines += 1
            w = rint / 100.0
            weights.append(w)
            hit = False
            for oi in _near(rtt, span + tol, obs_index):
                resid = float(obs[oi]["two_theta"]) - rtt
                if abs(resid) > span:
                    continue
                k = int(round((resid + span) / step))
                if k < 0 or k >= nbins:
                    continue
                hit = True
                local[k] += w              # вес голоса = интенсивность эталона
                local_sum[k] += w * resid
            if hit:
                voted += w
        if n_lines < SHIFT_MIN_LINES or sum(local) <= 0:
            continue
        # Когерентность = какая доля веса «голосующих» линий карточки собралась
        # в ОДНОМ бине. Считаем от линий, реально нашедших партнёра: линия,
        # которой в образце нет, не должна ни подтверждать, ни опровергать
        # юстировку. Иначе карточка с точной половиной линий отбраковывалась.
        ref_weight = max(sum(weights), 1e-9)
        if max(local) < SHIFT_MIN_COHERENCE * max(voted, 1e-9):
            continue
        if max(local) < SHIFT_MIN_MATCH:
            continue
        # нормируем вклад фазы на её собственный вес: одна сильная карточка с
        # 12 линиями не должна перевешивать три реальные фазы
        norm = 1.0 / ref_weight
        for k in range(nbins):
            votes[k] += local[k] * norm
            wsum[k] += local_sum[k] * norm

    if sum(votes) <= 0:
        return result(0.0, 0.0, 0.0, 0.0, False, len(cand), nbins)

    # фон случайных совпадений — МЕДИАНА по всем бинам: среднее завышено
    # самим пиком и нулевым бином, из-за чего значимость занижалась.
    ordered = sorted(votes)
    background = ordered[len(ordered) // 2] if ordered else 0.0
    best_k = max(range(nbins), key=lambda k: votes[k])
    base, best = votes[zero_bin], votes[best_k]

    accepted = False
    shift = 0.0
    if best > 0:
        if best_k == zero_bin:
            accepted = True
        else:
            gain_ok = best >= base * (1.0 + SHIFT_GAIN)
            # значимость пика гистограммы над фоном: пуассоновский z-тест.
            # Фиксированный «×2 над фоном» не работает в плотных спектрах
            # (фон сам велик), а z-тест масштабируется с фоном.
            z = ((best - background) / math.sqrt(background + 0.25)
                 if best > background else 0.0)
            peak_ok = bool(best >= max(SHIFT_MIN_VOTES, 2.0 * background)
                           or z >= SHIFT_MIN_Z)
            accepted = bool(gain_ok and peak_ok)
    if accepted:
        shift = wsum[best_k] / votes[best_k]
        shift = max(-span, min(span, shift))
    # Выигрыш: относительное превышение над нулевым бином. Если в нулевом бине
    # голосов нет вовсе (весь спектр сдвинут), выигрыш считаем полным (1.0), а не
    # нулевым — иначе отчёт показывал «выигрыш 0 %» при найденном сдвиге.
    if base > 0:
        gain = min(best / base - 1.0, 1.0)     # потолок 100 %: при почти нулевом
    else:                                      # нулевом бине отношение теряет смысл
        gain = 1.0 if best > 0 else 0.0
    if not accepted:
        shift = 0.0
    return result(shift, gain, base, best, accepted, len(cand), nbins)


# --------------------------------------------------------- подготовка наблюдений
def _prepare_obs(peaks, min_obs_int=MIN_OBS_INT):
    """Отделить пики, участвующие в поиске: (obs, карта obs_index -> peak_index)."""
    obs, mapping = [], []
    try:
        threshold = float(min_obs_int)
    except (TypeError, ValueError):
        threshold = MIN_OBS_INT
    for i, p in enumerate(peaks or ()):
        try:
            inten = float(p.get("intensity", 0.0))
        except (TypeError, ValueError):
            inten = 0.0
        if inten >= threshold:
            obs.append(p)
            mapping.append(i)
    return obs, mapping


def _check_db_filter(db_filter):
    value = str(db_filter or "all").strip().lower()
    if value not in ("cod", "all"):
        raise ValueError("db_filter: ожидалось 'cod'|'all', получено '{}'".format(db_filter))
    return value


def _prior_info(ph, use_profile):
    """(список пометок prior, вес профиля 0..1)."""
    if not use_profile or _profile is None:
        return [], 0.0
    tags = _profile.prior_tags(ph)
    if not tags:
        return [], 0.0
    return tags, _profile.prior_weight(ph)


# ------------------------------------------------------------------ конвейер
def _evaluate(peaks, db, tol, elements_only, shift, shift_mode, two_theta_min,
              two_theta_max, radiation, db_filter, use_profile, min_obs_int,
              max_ref_lines, profile_bonus=PRIOR_BONUS, auto_shift=None):
    """Общий конвейер: диапазон, сдвиг, кандидаты, скоринг (без выбора фаз).

    Возвращает dict: ``scored`` (отсортировано по рангу), ``obs``, ``map_obs``,
    ``context``. Все индексы в ``scored[*]["matched"]`` — в пространстве
    полного списка ``peaks``.
    """
    db_filter = _check_db_filter(db_filter)
    rad = _make_radiation(radiation)
    phases = db.get("phases", ())
    mode, shift_val = _resolve_shift(shift, shift_mode)
    obs_all, map_all = _prepare_obs(peaks, min_obs_int)
    n_below = len(peaks or ()) - len(obs_all)

    tmin, tmax = two_theta_min, two_theta_max
    declared = two_theta_min is not None or two_theta_max is not None
    if obs_all:
        obs_tt = [float(p["two_theta"]) for p in obs_all]
        if tmin is None:
            tmin = min(obs_tt)
        if tmax is None:
            tmax = max(obs_tt)

    outside = []
    if obs_all and declared and tmin is not None and tmax is not None:
        keep = [k for k, p in enumerate(obs_all)
                if float(tmin) - 1e-9 <= float(p["two_theta"]) <= float(tmax) + 1e-9]
        keep_set = set(keep)
        outside = [map_all[k] for k in range(len(obs_all)) if k not in keep_set]
        obs = [obs_all[k] for k in keep]
        mapping = [map_all[k] for k in keep]
    else:
        obs, mapping = obs_all, map_all

    if obs:
        lo, hi = _scan_bounds(obs, tmin, tmax, tol)
    else:
        lo = hi = None

    gain = 0.0
    shift_detail = None
    if auto_shift is not None:
        # сдвиг найден на втором проходе (по лучшим фазам) — уже применён
        shift_val = float(auto_shift[0])
        shift_detail = auto_shift[1]
        gain = float(shift_detail.get("gain") or 0.0)
    elif mode == "auto":
        if obs:
            shift_detail = estimate_shift(obs, db, tol=tol, min_obs_int=min_obs_int,
                                          two_theta_min=tmin, two_theta_max=tmax,
                                          radiation=rad, db_filter=db_filter,
                                          detail=True)
            shift_val = shift_detail["shift"]
            gain = shift_detail["gain"] if shift_detail["accepted"] else 0.0
        else:
            shift_val, gain = 0.0, 0.0
    context = {
        "radiation": rad,
        "shift": float(shift_val or 0.0),
        "shift_mode": mode,
        "shift_gain": float(gain or 0.0),
        "shift_detail": shift_detail,
        "two_theta_min": None if tmin is None else round(float(tmin), 4),
        "two_theta_max": None if tmax is None else round(float(tmax), 4),
        "n_peaks": len(peaks or ()),
        "n_obs": len(obs),
        "n_below_threshold": n_below,
        "n_outside_range": len(outside),
        "outside_range": outside,
        "db_filter": db_filter,
        "use_profile": bool(use_profile),
        "n_candidates": 0,
        "peak_density": 0.0,
        "p_chance": 0.0,
    }
    if obs and lo is not None and hi is not None and hi > lo:
        density = len(obs) / (hi - lo)
        context["peak_density"] = round(density, 3)
        context["p_chance"] = round(min(1.0, 2.0 * tol * density), 4)

    scored = []
    if not obs:
        return {"scored": scored, "obs": obs, "map_obs": mapping, "context": context}

    index, idx_step = _angle_index(db, rad)
    cand = _candidate_phases(obs, index, idx_step,
                             tol + abs(context["shift"]), shift=context["shift"])
    context["n_candidates"] = len(cand)
    obs_index = _obs_index(obs)
    obs_int = [max(float(p["intensity"]), 0.5) for p in obs]
    obs_total = sum(obs_int) or 1.0

    for pi in cand:
        ph = phases[pi]
        if db_filter == "cod" and str(ph.get("src", "COD")).upper() != "COD":
            continue
        if elements_only and not set(elements_only) & set(ph.get("elements") or ()):
            continue
        sc = _score_phase(ph, obs, tol, shift=context["shift"],
                          max_ref_lines=max_ref_lines,
                          two_theta_min=tmin, two_theta_max=tmax,
                          radiation=rad, obs_index=obs_index, n_tested=len(cand))
        if sc is None:
            continue
        prior, weight = _prior_info(ph, use_profile)
        explained_obs = {}
        for ri, oi, dth in sc["matched"]:
            explained_obs[oi] = (ri, dth)
        sc["matched"] = [(ri, mapping[oi], dth) for ri, oi, dth in sc["matched"]]
        sc["explained_peaks"] = sorted(mapping[oi] for oi in explained_obs)
        sc["share"] = round(sum(obs_int[oi] for oi in explained_obs) / obs_total, 4)
        sc["phase"] = ph
        sc["phase_index"] = pi
        sc["prior"] = prior
        sc["rank_score"] = round(sc["fom"] + float(profile_bonus) * weight, 3)
        scored.append(sc)

    scored.sort(key=lambda r: (-r["rank_score"], -r["share"], -r["fom"]))

    # --- сдвиг, второй проход -------------------------------------------------
    # Гистограмма по всей базе в больших библиотеках тонет в фоне случайных
    # карточек: 32k карточек дают десятки «голосов» в каждый бин. Если первый
    # проход сдвиг не объявил, повторяем голосование ТОЛЬКО по лучшим фазам
    # (они уже прошли скоринг) и при значимом результате пересчитываем метрики.
    if (mode == "auto" and auto_shift is None and scored
            and not (shift_detail and shift_detail["accepted"])):
        top_pi = {sc["phase_index"] for sc in scored[:SHIFT_TOP_PHASES]}
        again = estimate_shift(obs, db, tol=tol, min_obs_int=min_obs_int,
                               two_theta_min=tmin, two_theta_max=tmax,
                               radiation=rad, db_filter=db_filter,
                               detail=True, restrict=top_pi)
        if again.get("accepted"):
            return _evaluate(peaks, db, tol, elements_only, again["shift"], "manual",
                             two_theta_min, two_theta_max, rad, db_filter,
                             use_profile, min_obs_int, max_ref_lines,
                             profile_bonus, auto_shift=(again["shift"], again))

    if auto_shift is not None:
        context["shift_mode"] = "auto"
        context["shift_gain"] = float(shift_detail.get("gain") or 0.0)
        context["shift_detail"] = shift_detail

    return {"scored": scored, "obs": obs, "map_obs": mapping, "context": context}


def _pattern_lines(ph, rad=None):
    """Сильные линии карточки (I >= 20 % от максимума) в текущем излучении."""
    refs = _ref_positions(ph, rad)
    if not refs:
        return []
    imax = max(r[1] for r in refs) or 1.0
    return [tt for tt, inten in refs if inten >= 0.2 * imax]


def _heavy_composition(ph):
    """Состав без водорода: CaH4SO6 и CaSO6 — одна и та же фаза для PXRD."""
    formula = str(ph.get("formula") or "")
    out = []
    i = 0
    while i < len(formula):
        ch = formula[i]
        if ch.isupper():
            j = i + 1
            while j < len(formula) and formula[j].islower():
                j += 1
            el = formula[i:j]
            k = j
            while k < len(formula) and (formula[k].isdigit() or formula[k] == '.'):
                k += 1
            num = formula[j:k]
            if el != "H":
                out.append(el + num)
            i = k
        else:
            # «хвост» без элемента (Gg2Hh5b): часть обозначения, не теряем
            j = i
            while j < len(formula) and not formula[j].isupper():
                j += 1
            out.append(formula[i:j])
            i = j
    return "".join(sorted(out))


def _same_pattern(ph1, ph2, rad=None, tol=SAME_PATTERN_TOL, frac=SAME_PATTERN_FRAC):
    """Дубль описания одной фазы: тот же состав без H и та же порошковая картина.

    В базе гипс описан и как CaH4SO6, и как CaSO6; показывать оба как разные
    кандидаты нельзя — это удваивает «подтверждение» и вытесняет настоящие фазы
    из результата. Требуем совпадения и состава (без водорода — рентген его не
    видит), и сильных линий: иначе изоструктурные, но разные фазы (кварц/BeF2)
    склеились бы в одну.
    """
    if _heavy_composition(ph1) != _heavy_composition(ph2):
        return False
    a, b = _pattern_lines(ph1, rad), _pattern_lines(ph2, rad)
    if len(a) < 3 or len(b) < 3:
        return False
    if len(a) > len(b):
        a, b = b, a
    hit = sum(1 for tt in a if any(abs(tt - t) <= tol for t in b))
    return hit / float(len(a)) >= frac


def _select(scored, min_fom, max_phases, radiation=None):
    """Отобрать фазы: порог FOM, без дублей карточек одной формулы и одной картины."""
    out = []
    per_formula = {}
    limit = max(1, int(max_phases))
    for sc in scored:
        if sc["fom"] < min_fom:
            continue
        formula = sc["phase"].get("formula") or ""
        kept = per_formula.setdefault(formula, [])
        if len(kept) >= MAX_PER_FORMULA:
            continue
        if any(abs(sc["fom"] - prev) <= DUP_FOM_TOL for prev in kept):
            continue
        twin = next((other for other in out
                     if _same_pattern(sc["phase"], other["phase"], radiation)), None)
        if twin is not None:
            # одна фаза описана дважды: оставляем более полную карточку
            if sc["n_matched"] > twin["n_matched"]:
                out[out.index(twin)] = sc
                per_formula[twin["phase"].get("formula") or ""] = [
                    v for v in per_formula[twin["phase"].get("formula") or ""] if v != twin["fom"]]
                kept = per_formula.setdefault(formula, [])
                kept.append(sc["fom"])
            continue
        kept.append(sc["fom"])
        out.append(sc)
        if len(out) >= limit:
            break
    return out


def _weak_claimants(chosen, scored, obs, mapping, min_fom):
    """Фазы «недостаточно данных» для сильных необъяснённых пиков.

    Если сильный пик образца не объяснён ни одним кандидатом, а какая-то фаза
    совпала с ним своей единственной линией — об этом надо сообщить (status
    ``weak``, флаг FEW_DIAGNOSTIC_LINES), иначе сильнейший пик молча остаётся
    «не объяснён». Такие фазы НЕ считаются подтверждёнными и не попадают в
    элементы; их FOM остаётся низким. Ограничение — MAX_WEAK_REPORTED.
    """
    if not chosen:
        return []
    claimed = set()
    for sc in chosen:
        claimed.update(pi for _ri, pi, _d in sc["matched"])
    imax = max((float(o["intensity"]) for o in obs), default=0.0) or 1.0
    orphans = [pi for k, pi in enumerate(mapping)
               if pi not in claimed and float(obs[k]["intensity"]) >= WEAK_ANCHOR * imax]
    if not orphans:
        return []
    out = []
    taken = {sc["phase"].get("id") for sc in chosen}
    peak_int = {}
    for k, o in enumerate(obs):
        peak_int[mapping[k]] = max(peak_int.get(mapping[k], 0.0), float(o["intensity"]))
    for pidx in sorted(orphans, key=lambda p: -peak_int.get(p, 0.0)):
        for sc in scored:
            if sc["fom"] >= min_fom or "FEW_DIAGNOSTIC_LINES" not in sc["flags"]:
                continue
            if sc["phase"].get("id") in taken:
                continue
            if pidx not in [pi for _ri, pi, _d in sc["matched"]]:
                continue
            if any(_same_pattern(sc["phase"], other["phase"], None) for other in out):
                continue
            out.append(sc)
            taken.add(sc["phase"].get("id"))
            break
        if len(out) >= MAX_WEAK_REPORTED:
            break
    return out


def _phase_status(sc, min_fom, min_share):
    """Статус фазы: candidate | ambiguous | weak."""
    flags = sc.get("flags") or []
    if ("FEW_DIAGNOSTIC_LINES" in flags or "SINGLE_LINE_MATCH" in flags
            or "INTENSITY_CONFLICT" in flags):
        return "weak"
    if sc.get("n_unique_peaks", 0) <= 0:
        return "ambiguous"
    if sc["fom"] < min_fom or sc.get("n_matched", 0) < 2:
        return "weak"
    if min_share and sc.get("share", 0.0) < min_share:
        return "weak"
    return "candidate"


def _settings(context, tol, min_fom, min_obs_int, use_profile, max_phases,
              min_share, max_ref_lines):
    return {
        "tol": tol,
        "min_fom": min_fom,
        "min_obs_int": min_obs_int,
        "max_phases": max_phases,
        "min_share": min_share,
        "max_ref_lines": max_ref_lines,
        "use_profile": bool(use_profile),
        "db_filter": context["db_filter"],
        "shift_mode": context["shift_mode"],
        "two_theta_min": context["two_theta_min"],
        "two_theta_max": context["two_theta_max"],
    }


def search_phases(peaks, db, tol=TOL_DEFAULT, min_fom=45.0, top=30,
                  elements_only=None, shift=None, shift_mode=None,
                  two_theta_min=None, two_theta_max=None, radiation=None,
                  wavelength=None, db_filter="all", use_profile=False,
                  min_obs_int=MIN_OBS_INT, max_ref_lines=MAX_REF_LINES):
    """Ранжированный список фаз-кандидатов по всему набору пиков.

    ``shift``: ``None`` -> авто, число -> ручной сдвиг, ``'off'`` -> без сдвига;
    ``shift_mode`` перекрывает это явно. Индексы в ``matched`` — в пространстве
    полного списка ``peaks``.
    """
    rad = radiation if radiation is not None else _make_radiation(None, wavelength)
    ev = _evaluate(peaks, db, tol, elements_only, shift, shift_mode,
                   two_theta_min, two_theta_max, rad, db_filter, use_profile,
                   min_obs_int, max_ref_lines)
    chosen = _select(ev["scored"], min_fom, top, ev["context"]["radiation"])
    for sc in chosen:
        sc.setdefault("n_unique_peaks", sc["n_matched"])
        sc["status"] = _phase_status(sc, min_fom, 0.0)
    return chosen


def explain(peaks, db, tol=TOL_DEFAULT, min_fom=55.0, max_phases=MAX_PHASES,
            min_share=0.25, shift=None, shift_mode=None, two_theta_min=None,
            two_theta_max=None, radiation=None, wavelength=None, db_filter="all",
            use_profile=False, min_obs_int=MIN_OBS_INT,
            max_ref_lines=MAX_REF_LINES, profile_bonus=PRIOR_BONUS):
    """Совместно объяснить набор пиков набором фаз (без взаимоисключения).

    Каждая фаза привязывается независимо: пик может принадлежать нескольким
    фазам. Фаза не исчезает из результата только потому, что её линии «заняла»
    другая фаза — она получает статус ``ambiguous``.

    ``min_share`` больше не отбрасывает фазы (это теряло малые фазы в смесях),
    а понижает статус до ``weak`` при малой доле объяснённой интенсивности.
    """
    rad = radiation if radiation is not None else _make_radiation(None, wavelength)
    ev = _evaluate(peaks, db, tol, None, shift, shift_mode, two_theta_min,
                   two_theta_max, rad, db_filter, use_profile,
                   min_obs_int, max_ref_lines, profile_bonus)
    scored = ev["scored"]
    obs, mapping, context = ev["obs"], ev["map_obs"], ev["context"]
    chosen = _select(scored, min_fom, max_phases, context["radiation"])
    chosen += _weak_claimants(chosen, scored, obs, mapping, min_fom)
    # Кандидаты ниже порога: нужны, чтобы отчёт не выглядел «пустым» и было видно,
    # ЧТО именно не дотянуло до порога. Это НЕ идентификация: все помечены
    # status="below_threshold" и их флаги показываются в отчёте.
    chosen_keys = {(id(sc), sc["phase"].get("id")) for sc in chosen}
    weak_candidates = []
    for sc in scored:
        if len(weak_candidates) >= MAX_WEAK_CANDIDATES:
            break
        if (id(sc), sc["phase"].get("id")) in chosen_keys:
            continue
        if sc["fom"] < WEAK_CANDIDATE_MIN_FOM:
            continue
        if any(abs(sc["fom"] - w["fom"]) < 1e-9
               and sc["phase"].get("formula") == w["phase"].get("formula")
               for w in weak_candidates):
            continue
        sc = dict(sc)
        sc["status"] = "below_threshold"
        sc.setdefault("flags", [])
        weak_candidates.append(sc)
    n_peaks = len(peaks or ())
    searched = set(mapping)

    # --- совместная привязка: пик может принадлежать нескольким фазам ---------
    per_peak = {}
    for no, sc in enumerate(chosen, 1):
        for ri, pidx, dth in sc["matched"]:
            refs = sc["phase"].get("peaks") or ()
            rint = float(refs[ri][2]) if ri < len(refs) else 0.0
            per_peak.setdefault(pidx, []).append(PeakAssignment(
                phase_no=no,
                dtheta=round(float(dth), 3),
                ref_int=rint,
                overlap=False,
            ))
    for entries in per_peak.values():
        multi = len({e["phase_no"] for e in entries}) >= 2
        for entry in entries:
            entry["overlap"] = multi

    ambiguous = sorted({p for p, lst in per_peak.items()
                        if len({e["phase_no"] for e in lst}) >= 2})

    for no, sc in enumerate(chosen, 1):
        unique = 0
        for pidx in sc["explained_peaks"]:
            owners = {e["phase_no"] for e in per_peak.get(pidx, ())}
            if owners == {no}:
                unique += 1
        sc["n_unique_peaks"] = unique
        sc["status"] = _phase_status(sc, min_fom, min_share)
        if sc["status"] == "weak" and min_share and sc["share"] < min_share:
            if "LOW_SHARE" not in sc["flags"]:
                sc["flags"].append("LOW_SHARE")

    # --- статусы пиков и необъяснённые ---------------------------------------
    outside = set(context.get("outside_range") or ())
    peak_status = {}
    for i in range(n_peaks):
        if i in outside:
            peak_status[i] = "outside_range"
        elif i not in searched:
            peak_status[i] = "below_threshold"
        elif i in ambiguous:
            peak_status[i] = "ambiguous"
        elif i in per_peak:
            peak_status[i] = "assigned"
        else:
            peak_status[i] = "unexplained"
    unassigned = sorted(i for i in searched if i not in per_peak)

    elements = {}
    for sc in chosen:
        if sc["status"] == "weak":
            continue
        for el in sc["phase"].get("elements") or ():
            elements.setdefault(el, set()).add(sc["phase"].get("formula") or "")

    db_meta = dict(ev["context"].get("db_meta") or db.get("meta") or {})
    db_meta["n_phases_in_db"] = len(db.get("phases", ()))
    db_meta["db_filter"] = context["db_filter"]

    return {
        "phases": chosen,
        "weak_candidates": weak_candidates,
        "per_peak": per_peak,
        "unassigned": unassigned,
        "ambiguous": ambiguous,
        "peak_status": peak_status,
        "elements": {k: sorted(v) for k, v in elements.items()},
        "n_obs": context["n_obs"],
        "n_peaks": context["n_peaks"],
        "n_below_threshold": context["n_below_threshold"],
        "n_outside_range": context["n_outside_range"],
        "shift": context["shift"],
        "shift_mode": context["shift_mode"],
        "shift_gain": context["shift_gain"],
        "shift_detail": context["shift_detail"],
        "two_theta_min": context["two_theta_min"],
        "two_theta_max": context["two_theta_max"],
        "radiation": context["radiation"].as_dict() if hasattr(context["radiation"], "as_dict") else {},
        "settings": _settings(context, tol, min_fom, min_obs_int, use_profile,
                              max_phases, min_share, max_ref_lines),
        "db_meta": db_meta,
        "n_candidates": context["n_candidates"],
        "peak_density": context["peak_density"],
        "p_chance": context["p_chance"],
    }


# ------------------------------------------------------- записи о перекрытии пиков
class PeakAssignment(dict):
    """Кто претендует на пик: ``{"phase_no", "dtheta", "ref_int", "overlap"}``.

    Наследник dict, но дополнительно ведёт себя как старый кортеж
    ``(phase_no, dtheta)``: старые потребители (identify.py, xlsx_out.py)
    делали ``sorted(ass)[0][0]``, и такой код продолжает работать без правок.
    """

    __slots__ = ()

    def _pair(self):
        return (int(self.get("phase_no", 0)), float(self.get("dtheta", 0.0)))

    def __lt__(self, other):
        if isinstance(other, PeakAssignment):
            return self._pair() < other._pair()
        return self._pair() < tuple(other)[:2]

    def __getitem__(self, key):
        if isinstance(key, int) and not isinstance(key, bool):
            return self._pair()[key]
        return dict.__getitem__(self, key)


# --------------------------------------------------------------------- отчёт
def _radiation_text(rad):
    anode = rad.get("anode") or "?"
    model = "Kα1" if rad.get("model") == "ka1" else "Kα"
    lam = rad.get("wavelength_A") or 0.0
    src = {"user": "задано пользователем", "metadata": "из метаданных прибора",
           "default": "ПРЕДПОЛОЖЕНИЕ (метаданные не дали ответа)",
           "wavelength": "задана длина волны"}.get(rad.get("source"), rad.get("source") or "?")
    return "{} {} λ = {:.5f} Å; источник: {}".format(anode, model, lam, src)


def _plural(n, one, few, many):
    n = int(n) % 100
    if 11 <= n <= 19:
        return many
    n = n % 10
    if n == 1:
        return one
    if 2 <= n <= 4:
        return few
    return many


def _shift_text(res):
    shift = float(res.get("shift") or 0.0)
    mode = res.get("shift_mode") or "manual"
    gain = float(res.get("shift_gain") or 0.0)
    detail = res.get("shift_detail") or {}
    if mode == "off":
        return "Сдвиг 2θ: не искался (shift_mode=off) — поправка юстировки не применялась"
    if mode == "manual":
        if abs(shift) < 1e-9:
            return "Сдвиг 2θ: 0.00° (режим manual) — поправка не задана"
        return ("Сдвиг 2θ: {:+.2f}° (задан вручную, режим manual); "
                "поправка нуля шкалы/высоты образца".format(shift))
    raw_gain = float(detail.get("gain", gain) or 0.0)
    if detail.get("reason") == "dense_peaks":
        return ("Сдвиг 2θ: не оценивался (режим auto) — пиков слишком много "
                "(плотность {:.2f} на градус), оценка юстировки ненадёжна; "
                "сначала проверьте пик-пикинг".format(
                    float(res.get("peak_density") or 0.0)))
    if abs(shift) > 1e-9:
        return ("Сдвиг 2θ: {:+.2f}° (режим auto, выигрыш {:.1f}% при пороге {:.0f}%); "
                "это диагностика юстировки, а не подгонка структуры").format(
                    shift, raw_gain * 100.0, SHIFT_GAIN * 100.0)
    return ("Сдвиг 2θ: не применялся (режим auto; лучший вариант {:+.2f}°, выигрыш {:.1f}% "
            "< порога {:.0f}%) — систематического смещения не обнаружено").format(
                float(detail.get("shift", 0.0) or 0.0), raw_gain * 100.0, SHIFT_GAIN * 100.0)


def _phase_name(ph):
    if _profile is not None:
        label = _profile.label(ph)
        if label:
            return label
    name = (ph.get("name") or "").strip()
    if name and name.replace("-", "").replace("_", "").replace(".", "").isalnum():
        return name
    return ph.get("formula") or "?"


def _assignment(entry):
    """Нормализовать запись per_peak (dict нового формата или старый кортеж)."""
    if isinstance(entry, dict):
        return (int(entry.get("phase_no", 0)), float(entry.get("dtheta", 0.0)),
                float(entry.get("ref_int", 0.0)), bool(entry.get("overlap")))
    return (int(entry[0]), float(entry[1]), 0.0, False)


def format_report(res, peaks, filename="", max_cand=3):
    """Текстовый отчёт для консоли и .txt (структура результата ``explain``)."""
    res = res or {}
    peaks = peaks or []
    settings = res.get("settings") or {}
    tol = float(settings.get("tol", TOL_DEFAULT))
    min_obs = float(settings.get("min_obs_int", MIN_OBS_INT))
    min_fom = float(settings.get("min_fom", 55.0))
    phases = res.get("phases") or []
    per_peak = res.get("per_peak") or {}
    peak_status = res.get("peak_status") or {}
    unassigned = list(res.get("unassigned") or [])
    ambiguous = list(res.get("ambiguous") or [])

    lines = []
    if filename:
        lines.append("Файл: {}".format(filename))
    rad = res.get("radiation") or {}
    if rad:
        lines.append("Излучение: {}".format(_radiation_text(rad)))
    tmin, tmax = res.get("two_theta_min"), res.get("two_theta_max")
    if tmin is not None and tmax is not None:
        lines.append("Диапазон съёмки: {:.2f}–{:.2f}° 2θ; сравнение с эталоном: "
                     "{:.2f}–{:.2f}° 2θ".format(tmin, tmax, tmin - tol, tmax + tol))
    n_peaks = int(res.get("n_peaks", len(peaks)))
    n_obs = int(res.get("n_obs", n_peaks))
    n_outside = int(res.get("n_outside_range", 0))
    lines.append("Пиков найдено: {}; участвовало в поиске: {} (I ≥ {:.1f}%); "
                 "ниже порога: {}{}".format(
                     n_peaks, n_obs, min_obs,
                     int(res.get("n_below_threshold", n_peaks - n_obs)),
                     "; вне объявленного диапазона: {}".format(n_outside) if n_outside else ""))
    density = float(res.get("peak_density") or 0.0)
    p_chance = float(res.get("p_chance") or 0.0)
    if density > 0:
        lines.append("Плотность пиков: {:.2f} на градус 2θ; ожидаемая доля случайных "
                     "совпадений при tol {:.2f}°: {:.0f}%".format(density, tol,
                                                                 100.0 * p_chance))
    if p_chance >= 0.5:
        lines.append("  ВНИМАНИЕ: пиков слишком много для такого диапазона — совпадения "
                     "с эталоном в основном случайны. Проверьте пик-пикинг (шум/фон), "
                     "прежде чем доверять фазам.")
    lines.append(_shift_text(res))

    db_meta = res.get("db_meta") or {}
    n_db = db_meta.get("n_phases_in_db") or db_meta.get("count")
    if n_db:
        lines.append("База эталонов: {} {}, фильтр '{}'; профиль строительных "
                     "материалов: {}".format(
                         n_db, _plural(n_db, "фаза", "фазы", "фаз"),
                         settings.get("db_filter", "all"),
                         "включён (мягкий приор)" if settings.get("use_profile") else "выключен"))
    lines.append("FOM — оценка совпадения, НЕ вероятность; количественный состав "
                 "по PXRD не определяется.")
    lines.append("")

    counts = {"candidate": 0, "ambiguous": 0, "weak": 0}
    for sc in phases:
        counts[sc.get("status", "candidate")] = counts.get(sc.get("status", "candidate"), 0) + 1
    lines.append("ФАЗЫ (ранжирование; статусы: candidate {}, ambiguous {}, weak {}):".format(
        counts.get("candidate", 0), counts.get("ambiguous", 0), counts.get("weak", 0)))
    if not phases:
        lines.append("  ничего не идентифицировано — слабый сигнал или мало пиков; "
                     "ослабьте min_fom / проверьте юстировку и базу")
    for i, sc in enumerate(phases, 1):
        ph = sc.get("phase") or {}
        status = sc.get("status") or "candidate"
        lines.append("  {}. {}  [{}]  {}  {}".format(
            i, _phase_name(ph), ph.get("formula", "?"),
            sc.get("datasource") or _datasource(ph), ph.get("id", "")))
        lines.append("     FOM {:5.1f} (до штрафов {:5.1f}), ранг {:5.1f}; статус: {}".format(
            float(sc.get("fom", 0.0)), float(sc.get("fom_raw", 0.0)),
            float(sc.get("rank_score", sc.get("fom", 0.0))), status))
        lines.append("     линии: совпало {}/{} в диапазоне (в карточке всего {}, "
                     "вне диапазона {}); ключевых {}/{}".format(
                         sc.get("n_matched", 0), sc.get("n_ref", 0),
                         sc.get("n_ref_total", 0), sc.get("n_ref_out_of_range", 0),
                         sc.get("key_hits", 0), sc.get("key_total", 0)))
        if sc.get("cov_ref") is not None:
            lines.append("     покрытие эталона {:.2f} (случайный уровень {:.2f} → значимое "
                         "{:.2f}); сильных пиков образца объяснено {:.0%}".format(
                             float(sc.get("cov_ref", 0.0)), float(sc.get("p_chance", 0.0)),
                             float(sc.get("cov_excess", 0.0)),
                             float(sc.get("obs_strong", 0.0))))
        if sc.get("expected_false") is not None:
            lines.append("     ожидаемое число случайных фаз базы с таким совпадением: "
                         "{:.3g}; штраф за это {:+.1f}".format(
                             float(sc.get("expected_false", 0.0)),
                             -float(sc.get("chance_penalty", 0.0))))
        lines.append("     пиков объяснено {} (только этой фазой {}), доля "
                     "интенсивности {:.2f}, средн. Δ2θ {:.3f}°".format(
                         len(sc.get("explained_peaks") or []),
                         sc.get("n_unique_peaks", 0), float(sc.get("share", 0.0)),
                         float(sc.get("mean_dtheta", 0.0))))
        extra = []
        elements = ph.get("elements") or []
        if elements:
            extra.append("элементы: " + ", ".join(elements))
        if sc.get("prior"):
            extra.append("prior: " + ", ".join(sc["prior"]))
        if extra:
            lines.append("     " + "; ".join(extra))
        flags = sc.get("flags") or []
        if flags:
            lines.append("     пометки: " + "; ".join(_FLAGS_TEXT.get(f, f) for f in flags))
    lines.append("")

    if res.get("elements"):
        lines.append("ЭЛЕМЕНТЫ (расчётно, только по фазам со статусом candidate/ambiguous; "
                     "XRD определяет фазы, не элементы):")
        for el in sorted(res["elements"]):
            lines.append("  {:3s} — {}".format(el, ", ".join(res["elements"][el][:4])))
    else:
        lines.append("ЭЛЕМЕНТЫ: не определены")
    lines.append("")

    if ambiguous:
        lines.append("ПЕРЕКРЫТИЯ (пик принадлежит >= 2 фазам; ни одна фаза не «заняла» пик монопольно):")
        for pidx in ambiguous[:24]:
            owners = sorted({e["phase_no"] for e in per_peak.get(pidx, ())})
            tt = peaks[pidx]["two_theta"] if pidx < len(peaks) else float("nan")
            names = []
            for no in owners:
                sc = phases[no - 1] if 0 < no <= len(phases) else None
                names.append("#{} {}".format(no, _phase_name(sc["phase"]) if sc else "?"))
            lines.append("  {:6.2f}° — {}".format(tt, ", ".join(names)))
        lines.append("")

    lines.append("ПИКИ -> ФАЗА:")
    for i, p in enumerate(peaks):
        status = peak_status.get(i)
        ass = per_peak.get(i)
        if ass:
            pairs = sorted((_assignment(e) for e in ass), key=lambda e: (e[0], -e[2]))[:max_cand]
            tag = ", ".join("#{} {}".format(no, _phase_name(phases[no - 1]["phase"])
                                            if 0 < no <= len(phases) else "?")
                            for no, _dth, _ri, _ov in pairs)
            overlapped = any(ov for _no, _dth, _ri, ov in (_assignment(e) for e in ass))
            if overlapped:
                tag += "  (перекрытие)"
        elif status == "below_threshold":
            tag = "· ниже порога I < {:.1f}%".format(min_obs)
        elif status == "outside_range":
            tag = "· вне объявленного диапазона съёмки"
        elif status == "unexplained" or i in unassigned:
            tag = "— не объяснён"
        else:
            tag = "—"
        fwhm = "  FWHM {:.3f}°".format(p["fwhm"]) if p.get("fwhm") else ""
        lines.append("  {:6.2f}°  d={:6.3f} Å  I={:5.1f}%{}   {}".format(
            p["two_theta"], p["d"], p["intensity"], fwhm, tag))
    if unassigned:
        lines.append("  не объяснены пики: {}".format(
            ", ".join("{:.2f}°".format(peaks[i]["two_theta"]) for i in unassigned[:12])))
    lines.append("")
    lines.append("Ограничения: результат — предварительная идентификация кандидатов по "
                 "расчётным эталонам (не картотека PDF/ICDD), без независимой валидации "
                 "(эталон, РФА, ТГ). FOM не является вероятностью.")
    return "\n".join(lines)
