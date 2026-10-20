#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Разбор файлов .raw дифрактометра Thermo ARL X'TRA (ПО WinXRD 2.0).

Блок данных = N записей по 16 байт, little-endian:
    int32 интенсивность | int32 служебный | double 2theta
Перед блоком лежит запись-заголовок (0, N, стартовый_угол), N = число точек.
Библиотек не требует — только стандартный Python.
"""

import math
import re
import struct
from pathlib import Path

LAMBDA_CU_KA1 = 1.5406  # Cu K-alpha1, A


def find_block(b):
    """Возвращает (данные, шаг, служебное_поле) или None.

    Блок ищем по сигнатуре: запись-заголовок (0, N, стартовый_угол), сразу за
    ней N записей по 16 байт с ровным шагом по 2θ. Требуем, чтобы первые два
    поля заголовка были именно (0, N) — иначе в бинарный поток случайно
    попадают записи данных (там второе поле = dwell) и разбор уходит не туда.
    """
    best = None
    for off in range(0, len(b) - 16):
        _, n, start = struct.unpack_from("<iid", b, off)
        if not (100 < n < 500000 and 1.0 < start < 90.0):
            continue
        if off + 16 * (n + 1) > len(b):
            continue
        head0, head1, _ = struct.unpack_from("<iid", b, off)
        if head0 != 0 or head1 != n:
            continue
        pts = [struct.unpack_from("<iid", b, off + 16 * k) for k in range(n + 1)]
        # pts[0] — заголовок (0, N, start); сами точки начинаются с pts[1]
        ang = [p[2] for p in pts[1:]]
        if len(ang) < 3:
            continue
        diffs = [round(ang[i + 1] - ang[i], 6) for i in range(len(ang) - 1)]
        # модальный шаг: в конце блока прибор пишет служебную запись, её отсекаем
        from collections import Counter
        step, cnt = Counter(diffs).most_common(1)[0]
        if step <= 0:
            continue
        keep = 1
        for k in range(len(diffs)):
            if diffs[k] != step:
                break
            keep += 1
        ang = ang[:keep]
        pts = pts[:keep + 1]
        if len(ang) < 100:
            continue
        span = ang[-1] - ang[0]
        if span < 1.0 or span > 180.0:
            continue
        cand = ([(p[2], p[0]) for p in pts[1:]], step, pts[1][1], len(ang))
        if best is None or cand[3] > best[3]:
            best = cand
    if best is None:
        return None
    return best[0], best[1], best[2]


def read_meta(path):
    """Вытаскивает из файла служебные строки прибора (анод, ПО, детектор)."""
    b = Path(path).read_bytes()[:4096]
    meta = {}
    for key, pattern in (("anode", rb"Copper|Cobalt|Molybdenum|Chromium"),
                         ("software", rb"WinXRD[ \x2d0-9.]*"),
                         ("detector", rb"Solid State \w+|Scintillation"),
                         ("atmosphere", rb"air|vacuum|nitrogen")):
        m = re.search(pattern, b)
        if m:
            meta[key] = m.group().decode("latin-1").strip()
    return meta


def parse(path):
    """-> (данные [(2theta, интенсивность)], шаг, метаданные)."""
    b = Path(path).read_bytes()
    found = find_block(b)
    if not found:
        raise ValueError("не найден блок данных — это не файл WinXRD .raw")
    data, step, dwell = found
    meta = read_meta(path)
    meta["step"] = step
    meta["dwell"] = dwell
    meta["points"] = len(data)
    return data, step, meta


def d_spacing(two_theta):
    """Межплоскостное расстояние, A (закон Вульфа-Брэгга)."""
    return LAMBDA_CU_KA1 / (2 * math.sin(math.radians(two_theta / 2)))


def find_peaks(data, min_intensity=None, top=8):
    """Локальные максимумы; возвращает [(2theta, интенсивность)] по убыванию интенсивности."""
    imax = max(i for _, i in data)
    if min_intensity is None:
        min_intensity = max(40, 0.03 * imax)
    peaks = []
    for k in range(2, len(data) - 2):
        window = [data[k + j][1] for j in range(-2, 3)]
        if data[k][1] != max(window) or data[k][1] < min_intensity:
            continue
        if peaks and data[k][0] - peaks[-1][0] < 0.35:
            if data[k][1] > peaks[-1][1]:
                peaks[-1] = [data[k][0], data[k][1]]
        else:
            peaks.append([data[k][0], data[k][1]])
    return sorted(peaks, key=lambda p: -p[1])[:top]


def write_xy(data, src, out=None):
    out = Path(out) if out else Path(src).with_suffix(".xy")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        for ang, inten in data:
            f.write("{:8.2f}{:12d}\n".format(ang, inten))
    return out


def write_csv(data, src, out=None):
    out = Path(out) if out else Path(src).with_suffix(".csv")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("two_theta_deg,d_A,intensity\n")
        for ang, inten in data:
            f.write("{:.2f},{:.4f},{}\n".format(ang, d_spacing(ang), inten))
    return out


def convert_data(path):
    """Полная конвертация в текстовые форматы. Возвращает (xy, csv, meta)."""
    data, step, meta = parse(path)
    return write_xy(data, path), write_csv(data, path), meta


def main(argv=None):
    """CLI: перетаскивание файлов или запуск из командной строки.

    python winxrd2xy.py файл.raw [ещё.raw ...]
    """
    import argparse
    import sys

    try:
        from console import setup as _setup_console

        _setup_console()
    except Exception:  # noqa: BLE001
        for stream in (sys.stdout, sys.stderr):
            try:
                stream.reconfigure(errors="replace")
            except Exception:  # noqa: BLE001
                pass
    ap = argparse.ArgumentParser(description="WinXRD .raw -> .xy + .csv")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    bad = 0
    for f in a.files:
        p = Path(f)
        if not p.is_file():
            print("нет файла: {}".format(p))
            bad += 1
            continue
        try:
            data, step, meta = parse(p)
            xy = write_xy(data, p)
            csv = write_csv(data, p)
        except Exception as e:  # noqa: BLE001
            print("ошибка {}: {}".format(p.name, e))
            bad += 1
            continue
        print("{}: точек {} ({:.2f} - {:.2f} град, шаг {}) -> {}, {}".format(
            p.name, len(data), data[0][0], data[-1][0], step,
            Path(xy).name, Path(csv).name))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
