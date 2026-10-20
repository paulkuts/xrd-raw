#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Поиск пиков на дифрактограмме (только numpy, без scipy).

Вход: два списка/массива — 2θ (град) и интенсивность.
Выход: список пиков с положением 2θ, d (Å), интенсивностью, I отн. (%) и FWHM.
"""

import math

import numpy as np

CU_KA1 = 1.54056   # Å, Cu Kα1
CU_KA = 1.5418     # Å, средневзвешенная Cu Kα


def d_from_2theta(two_theta, wavelength=CU_KA):
    """Межплоскостное расстояние по углу 2θ (град)."""
    th = math.radians(two_theta / 2.0)
    s = math.sin(th)
    if s <= 1e-9:
        return float("inf")
    return wavelength / (2.0 * s)


def two_theta_from_d(d, wavelength=CU_KA):
    """Угол 2θ (град) по межплоскостному расстоянию."""
    if d <= 0:
        return float("nan")
    x = wavelength / (2.0 * d)
    if x >= 1.0:
        return float("nan")
    return 2.0 * math.degrees(math.asin(x))


def _savgol_coeffs(window, poly):
    """Коэффициенты сглаживающего фильтра Савицкого–Голея (по центру)."""
    half = window // 2
    x = np.arange(-half, half + 1, dtype=float)
    A = np.vander(x, poly + 1, increasing=True)
    pinv = np.linalg.pinv(A)
    return pinv[0]


def smooth(y, window=7, poly=3):
    """Сглаживание Савицкого–Голея с сохранением краёв."""
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


def _prominences(y):
    """Prominence каждого локального максимума (упрощённый, но корректный)."""
    n = len(y)
    prom = np.zeros(n)
    # индексы локальных максимумов
    peaks = [i for i in range(1, n - 1) if y[i] >= y[i - 1] and y[i] > y[i + 1]]
    if y[0] > y[1]:
        peaks.insert(0, 0)
    if y[n - 1] >= y[n - 2]:
        peaks.append(n - 1)
    for i in peaks:
        # влево — до более высокой точки
        left = 0.0
        j = i - 1
        lmin = y[i]
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
        if j < 0:
            left = y[i]
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
         wavelength=CU_KA, max_peaks=120):
    """Найти пики.

    min_prominence — минимальная prominence в долях от максимума сигнала;
    min_height     — минимальная абсолютная высота в долях от максимума;
    smooth_window  — окно сглаживания (точек), 0/1 — без сглаживания.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 7:
        return []
    y = y - np.min(y)
    ymax = float(np.max(y))
    if ymax <= 0:
        return []

    if smooth_window and smooth_window > 2:
        ys = smooth(y, smooth_window)
    else:
        ys = y
    ys = np.clip(ys, 0, None)

    prom, idxs = _prominences(ys)
    peaks = []
    for i in idxs:
        if prom[i] < min_prominence * ymax:
            continue
        if ys[i] < min_height * ymax:
            continue
        # базовая линия под пиком: минимум в окрестности ±(окно поиска)
        half = max(3, int(round(0.25 / max(1e-6, (x[1] - x[0])))))  # ≈ ±0.25°
        lo, hi = max(0, i - half), min(len(y), i + half + 1)
        base = float(np.min(ys[lo:hi]))
        height = float(ys[i] - base)
        if height <= 0:
            continue
        fw = _fwhm(x, ys, i, height / 2.0)
        tt = _refine_position(x, ys, i)
        peaks.append({
            "two_theta": float(tt),
            "d": d_from_2theta(tt, wavelength),
            "height": height,
            "prominence": float(prom[i]),
            "fwhm": fw,
        })

    if not peaks:
        return []
    # убрать дубли ближе 0.05°
    peaks.sort(key=lambda p: p["two_theta"])
    merged = []
    for p in peaks:
        if merged and abs(p["two_theta"] - merged[-1]["two_theta"]) < 0.05:
            if p["height"] > merged[-1]["height"]:
                merged[-1] = p
            continue
        merged.append(p)

    imax = max(p["height"] for p in merged)
    for p in merged:
        p["intensity"] = 100.0 * p["height"] / imax
    merged.sort(key=lambda p: -p["intensity"])
    merged = merged[:max_peaks]
    merged.sort(key=lambda p: p["two_theta"])
    for k, p in enumerate(merged, 1):
        p["n"] = k
    return merged
