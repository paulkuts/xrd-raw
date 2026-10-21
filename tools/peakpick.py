#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Поиск пиков на дифрактограмме (только numpy, без scipy).

Вход: два списка/массива — 2θ (град) и интенсивность.
Выход: список пиков с положением 2θ, d (Å), интенсивностью, I отн. (%) и FWHM.

Пороги адаптивные. Раньше prominence/height считались только как доля от
максимума сигнала, поэтому на слабом шумном образце шум принимался за пики
(реальный файл 261007115241.raw: 120 «пиков» при ~14 значимых). Теперь:

  σ  — оценка шума по MAD разностей соседних точек (устойчива к пикам и фону);
  порог prominence = max(min_prominence * max(сигнал), prominence_snr * σ);
  порог высоты     = max(min_height     * max(сигнал), height_snr     * σ).

Старые kwarg-и и их дефолты сохранены: при малом шуме поведение прежнее,
при большом — порог поднимается вместе с шумом. `adaptive=False` возвращает
чисто «долевые» пороги.

Детекция идёт по копии сигнала: сглаживание и вычитание локального фона
(скользящий минимум) выполняются внутри функции, исходные данные не меняются.

Длина волны — из единого источника tools/radiation.py (по умолчанию Cu Kα
средневзвешенная 1.54184 Å, та же, в которой посчитана база эталонов).
"""

import math

import numpy as np

try:  # обычный запуск: tools/ лежит в sys.path
    import radiation as _rad
    from radiation import BASE_WAVELENGTH, d_from_two_theta as _rad_d
    from radiation import two_theta_from_d as _rad_tt
except ImportError:  # запуск как пакет tools.peakpick
    from . import radiation as _rad
    from .radiation import BASE_WAVELENGTH, d_from_two_theta as _rad_d
    from .radiation import two_theta_from_d as _rad_tt

# Совместимость: раньше константы жили здесь. Источник истины — radiation.py.
CU_KA1 = _rad.ANODES["Cu"]["ka1"]          # 1.54056 Å, Cu Kα1
CU_KA = BASE_WAVELENGTH                    # 1.54184 Å, средневзвешенная Cu Kα

# Пороги по шуму (в единицах σ). Пик должен быть виден над шумом:
# prominence >= 4.3σ и высота над локальным фоном >= 3σ. Порог 5σ/4σ срезал
# реальные слабые линии: на образце 261007115241.raw не находились 39.41°,
# 43.2°, 47.2°, 48.6° (3.2–3.5σ). После снижения до 4.3σ/3σ найдено 10 линий
# вместо 7, на чистом шуме по-прежнему 1 пик, на наклонном фоне — 3–4,
# синтетический тест «фон без пиков» даёт 0. Ниже 4.0σ появляется ложный пик
# на чистом гауссовом фоне (проверено), поэтому 4.3σ — нижняя безопасная граница.
PROMINENCE_SNR = 4.3
HEIGHT_SNR = 3.0
BASELINE_WINDOW = 1.2   # °2θ — окно скользящего минимума (локальный фон)
MERGE_DEG = 0.05        # °2θ — ближе этого пики считаются одним


def _resolve_wavelength(wavelength=None, radiation=None):
    """λ (Å): объект Radiation > явный wavelength > база (Cu Kα средневзвешенная)."""
    if radiation is not None:
        lam = getattr(radiation, "lam", radiation)
        try:
            lam = float(lam)
        except (TypeError, ValueError):
            lam = None
        if lam and lam > 0:
            return lam
    if wavelength is not None:
        return float(wavelength)
    return float(BASE_WAVELENGTH)


def d_from_2theta(two_theta, wavelength=None, radiation=None):
    """Межплоскостное расстояние по углу 2θ (град)."""
    return _rad_d(float(two_theta), _resolve_wavelength(wavelength, radiation))


def two_theta_from_d(d, wavelength=None, radiation=None):
    """Угол 2θ (град) по межплоскостному расстоянию."""
    return _rad_tt(float(d), _resolve_wavelength(wavelength, radiation))


# ------------------------------------------------------------------- обработка
def _savgol_coeffs(window, poly):
    """Коэффициенты сглаживающего фильтра Савицкого–Голея (по центру)."""
    half = window // 2
    x = np.arange(-half, half + 1, dtype=float)
    A = np.vander(x, poly + 1, increasing=True)
    pinv = np.linalg.pinv(A)
    return pinv[0]


def smooth(y, window=7, poly=3):
    """Сглаживание Савицкого–Голея с сохранением краёв (копия, не in-place)."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < window or window < 3:
        return y.copy()
    if window % 2 == 0:
        window += 1
    c = _savgol_coeffs(window, min(poly, window - 1))
    half = window // 2
    out = y.copy()
    out[half:n - half] = np.convolve(y, c[::-1], mode="same")[half:n - half]
    return out


def noise_sigma(y):
    """Оценка σ шума по разностям соседних точек (MAD, устойчиво к пикам).

    Для белого шума σ разности = σ·√2, медиана |разность| = 0.6745·σ·√2,
    отсюда σ = 1.4826·median(|Δ|)/√2. Пики и фон занимают малую долю точек,
    поэтому медиана их почти не замечает.
    """
    y = np.asarray(y, dtype=float)
    if y.size < 3:
        return 0.0
    d = np.diff(y)
    d = d[np.isfinite(d)]
    if d.size < 2:
        return 0.0
    med = float(np.median(d))
    mad = float(np.median(np.abs(d - med)))
    sigma = 1.4826 * mad / math.sqrt(2.0)
    if sigma <= 0.0:  # квантованные/округлённые данные: MAD вырождается в 0
        sigma = float(np.std(d)) / math.sqrt(2.0)
    return max(0.0, float(sigma))


estimate_noise = noise_sigma


def _rolling_min(y, window):
    """Скользящий минимум (локальный фон) через strided-представление numpy."""
    y = np.asarray(y, dtype=float)
    n = y.size
    w = int(window)
    if n == 0:
        return y.copy()
    if w < 3:
        return y.copy()
    if w >= n:
        return np.full(n, float(np.min(y)))
    if w % 2 == 0:
        w += 1
    half = w // 2
    padded = np.pad(y, (half, half), mode="edge")
    view = np.lib.stride_tricks.sliding_window_view(padded, w)
    return view.min(axis=1)


def _prominences(y):
    """Prominence каждого локального максимума (упрощённый, но корректный)."""
    n = len(y)
    prom = np.zeros(n)
    peaks = [i for i in range(1, n - 1) if y[i] >= y[i - 1] and y[i] > y[i + 1]]
    if n > 1 and y[0] > y[1]:
        peaks.insert(0, 0)
    if n > 1 and y[n - 1] >= y[n - 2]:
        peaks.append(n - 1)
    for i in peaks:
        lmin = y[i]
        j = i - 1
        while j >= 0:
            if y[j] > y[i]:
                break
            lmin = min(lmin, y[j])
            j -= 1
        rmin = y[i]
        j = i + 1
        while j < n:
            if y[j] > y[i]:
                break
            rmin = min(rmin, y[j])
            j += 1
        base = max(lmin, rmin)
        prom[i] = y[i] - base
    return prom, peaks


def _refine_position(x, y, i):
    """Параболическое уточнение положения максимума по трём точкам."""
    if i <= 0 or i >= len(y) - 1:
        return x[i]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    denom = (y0 - 2 * y1 + y2)
    if abs(denom) < 1e-12:
        return x[i]
    delta = 0.5 * (y0 - y2) / denom
    dx = x[1] - x[0] if len(x) > 1 else 0.0
    return x[i] + delta * dx


def _fwhm(x, y, i, half):
    """Ширина на половине высоты (в градусах 2θ) с линейной интерполяцией."""
    n = len(y)
    target = y[i] - half
    left = None
    for j in range(i, 0, -1):
        if y[j] <= target:
            y_hi, y_lo = y[j], y[j + 1]
            if abs(y_hi - y_lo) < 1e-12:
                left = x[j]
            else:
                frac = (target - y_hi) / (y_lo - y_hi)
                left = x[j] + frac * (x[j + 1] - x[j])
            break
    right = None
    for j in range(i, n - 1):
        if y[j] <= target:
            y_hi, y_lo = y[j], y[j - 1]
            if abs(y_hi - y_lo) < 1e-12:
                right = x[j]
            else:
                frac = (target - y_hi) / (y_lo - y_hi)
                right = x[j] - frac * (x[j] - x[j - 1])
            break
    if left is None or right is None or right <= left:
        return None
    return right - left


def pick(x, y, min_prominence=0.02, min_height=0.01, smooth_window=7,
         wavelength=None, max_peaks=120, radiation=None,
         prominence_snr=PROMINENCE_SNR, height_snr=HEIGHT_SNR,
         baseline_window=BASELINE_WINDOW, merge_deg=MERGE_DEG,
         adaptive=True, sigma_noise=None):
    """Найти пики.

    min_prominence — нижний порог prominence в долях от максимума сигнала;
    min_height     — нижний порог высоты в долях от максимума сигнала;
    smooth_window  — окно сглаживания (точек), 0/1 — без сглаживания;
    wavelength     — λ, Å; radiation — объект radiation.Radiation (приоритетнее);
    max_peaks      — сколько пиков вернуть максимум;
    prominence_snr, height_snr — пороги в σ шума (адаптивная часть);
    baseline_window — окно локального фона, °2θ;
    merge_deg      — пики ближе этого по 2θ объединяются;
    adaptive       — False: только «долевые» пороги (поведение до 2.1);
    sigma_noise    — задать σ шума вручную (иначе оценка по MAD разностей).

    Возвращает список dict по возрастанию 2θ: two_theta, d, height,
    prominence, fwhm, intensity, snr, wavelength, n (1-based), peak_id (0-based).
    Исходные x/y не изменяются.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        raise ValueError("длины x и y не совпадают: {} и {}".format(x.size, y.size))
    if x.size < 7:
        return []
    if not np.all(np.isfinite(y)):
        return []

    lam = _resolve_wavelength(wavelength, radiation)

    step = float(np.median(np.diff(x))) if x.size > 1 else 0.0
    if step <= 0:
        step = 0.02

    # --- копии для детекции: исходный сигнал не трогаем
    ys = smooth(y, smooth_window) if smooth_window and smooth_window > 2 else y.copy()

    win = int(round(float(baseline_window) / step)) if step > 0 else 31
    win = max(5, win | 1)
    # Локальный фон — скользящий минимум, но не ниже минимума исходной кривой.
    # Иначе одиночный выброс (космика, край скана) даёт «звон» Савицкого–Голея:
    # провал в сглаженной кривой, фон уходит в него, и рядом вырастают ложные
    # пики. Медианный фильтр здесь не годится — на слабом образце настоящий
    # пик шириной 1 канал (кварц 26.66°) неотличим от выброса, и медиана его
    # стирает.
    def _base(curve):
        rm = _rolling_min(curve, win)
        return smooth(rm, smooth_window) if smooth_window and smooth_window > 2 else rm

    baseline = np.maximum(_base(ys), _base(y))
    z = np.clip(ys - baseline, 0.0, None)

    ymax = float(np.max(z))
    if ymax <= 0.0:
        return []

    sigma = float(sigma_noise) if sigma_noise is not None else noise_sigma(y)
    prom_floor = float(min_prominence) * ymax
    height_floor = float(min_height) * ymax
    if adaptive:
        prom_floor = max(prom_floor, float(prominence_snr) * sigma)
        height_floor = max(height_floor, float(height_snr) * sigma)

    prom, idxs = _prominences(z)
    peaks = []
    for i in idxs:
        if i <= 0 or i >= z.size - 1:      # краевые точки не считаем пиками
            continue
        if prom[i] < prom_floor:
            continue
        height = float(z[i])
        if height < height_floor:
            continue
        fw = _fwhm(x, z, i, height / 2.0)
        if fw is not None:
            fw = float(fw)
            if fw <= 0:
                fw = None
        tt = _refine_position(x, z, i)
        peaks.append({
            "two_theta": float(tt),
            "d": d_from_2theta(tt, lam),
            "height": height,
            "prominence": float(prom[i]),
            "fwhm": fw,
            "snr": (float(prom[i]) / sigma) if sigma > 0 else None,
            "wavelength": lam,
        })

    if not peaks:
        return []

    # убрать дубли ближе merge_deg
    peaks.sort(key=lambda p: p["two_theta"])
    merged = []
    for p in peaks:
        if merged and abs(p["two_theta"] - merged[-1]["two_theta"]) < merge_deg:
            if p["height"] > merged[-1]["height"]:
                merged[-1] = p
            continue
        merged.append(p)

    imax = max(p["height"] for p in merged)
    for p in merged:
        p["intensity"] = 100.0 * p["height"] / imax
    if max_peaks and len(merged) > max_peaks:
        merged.sort(key=lambda p: -p["intensity"])
        merged = merged[:int(max_peaks)]
    merged.sort(key=lambda p: p["two_theta"])
    for k, p in enumerate(merged):
        p["peak_id"] = k
        p["n"] = k + 1
    return merged
