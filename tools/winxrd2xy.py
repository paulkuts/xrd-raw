#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Разбор файлов .raw дифрактометра Thermo ARL X'TRA (ПО WinXRD 2.0).

Блок данных = N записей по 16 байт, little-endian:
    int32 интенсивность | int32 служебный (время на точку) | double 2theta
Перед блоком лежит запись-заголовок (0, N, стартовый_угол). В реальных файлах
N включает последнюю служебную запись (её 2θ = 260° и она отбрасывается).

Длина волны берётся из единого источника `tools/radiation.py`: по умолчанию
Cu Kα средневзвешенная (1.54184 Å) — та же, в которой посчитана база эталонов.
Метаданные прибора дают анод, но НЕ дают спектральную модель (Kα1 / Kα1+Kα2),
поэтому модель здесь не угадывается.

Библиотек не требует — только стандартный Python.
"""

import re
import struct
from collections import Counter
from pathlib import Path

try:  # обычный запуск: tools/ лежит в sys.path (скрипт, меню, тесты)
    import radiation
    from radiation import BASE_WAVELENGTH, d_from_two_theta as _rad_d
except ImportError:  # запуск как пакет tools.winxrd2xy
    from . import radiation
    from .radiation import BASE_WAVELENGTH, d_from_two_theta as _rad_d

# Только для обратной совместимости внешних вызовов; источник истины — radiation.
LAMBDA_CU_KA1 = radiation.ANODES["Cu"]["ka1"]

DATA_RECORD = 16
HEADER = struct.Struct("<iid")
MIN_POINTS = 100          # меньше — это не полный скан
MAX_POINTS = 200000       # защита от мусорного заголовка (память)
MAX_FILE_BYTES = 64 * 1024 * 1024
MIN_SPAN = 0.5            # градусов 2θ
MAX_SPAN = 180.0
SERVICE_ANGLE = 170.0     # служебные записи прибора лежат далеко за рабочим углом
START_TOL = 5.0           # °: стартовый угол заголовка vs первый угол данных
MIN_ANGLE = 0.5           # °2θ: физически осмысленный минимум
_ZERO_RUN = re.compile(rb"\x00{4,}")   # заголовок всегда начинается с int32 0

# Сигнатуры чужих файлов для понятной ошибки вместо «не найден блок».
MAGIC = (
    (b"\x89PNG\r\n\x1a\n", "PNG"),
    (b"\xff\xd8\xff", "JPEG"),
    (b"GIF87a", "GIF"),
    (b"GIF89a", "GIF"),
    (b"%PDF", "PDF"),
    (b"PK\x03\x04", "ZIP/DOCX/XLSX"),
    (b"\x1f\x8b", "GZIP"),
    (b"II*\x00", "TIFF"),
    (b"MM\x00*", "TIFF"),
    (b"RIFF", "RIFF (WAV/AVI)"),
)


class _Reject(Exception):
    """Кандидат на блок данных не прошёл проверку."""

    def __init__(self, reason, detail):
        Exception.__init__(self, reason)
        self.reason = reason
        self.detail = detail


# --------------------------------------------------------------------- утилиты
def _err(name, msg):
    return ValueError("{}: {}".format(name, msg))


def _looks_like_image(b):
    for magic, label in MAGIC:
        if b.startswith(magic):
            return label
    if len(b) > 14 and b[:2] == b"BM" and b[6:10] == b"\x00\x00\x00\x00":
        return "BMP"
    return None


def _looks_like_text(b):
    """Текстовый файл (.xy/.csv/отчёт), а не бинарный .raw."""
    sample = b[:65536]
    if b"\x00" in sample:
        return False
    try:
        text = sample.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = sample.decode("cp1251")
        except UnicodeDecodeError:
            return False
    if not text:
        return False
    printable = sum(1 for ch in text if ch.isprintable() or ch in "\r\n\t")
    return printable >= 0.95 * len(text)


def _tol(step):
    return max(1e-6, 1e-3 * abs(step))


# ------------------------------------------------------------- проверка блока
def _validate(b, off, n, start):
    """Полная проверка кандидата; вернуть dict с точками или бросить _Reject."""
    end = off + DATA_RECORD * (n + 1)
    if end > len(b):
        have = max(0, (len(b) - off) // DATA_RECORD - 1)
        # заголовок обязан стоять рядом со своими данными: стартовый угол
        # и первый реально лежащий в файле угол отличаются не сильно
        if have >= 1:
            ang0 = struct.unpack_from("<d", b, off + DATA_RECORD + 8)[0]
            if not (MIN_ANGLE <= ang0 <= SERVICE_ANGLE) or abs(ang0 - start) > START_TOL:
                raise _Reject("mismatch", {"start": start, "ang0": ang0})
        raise _Reject("truncated", {"declared": n, "available": have,
                                    "start": start, "offset": off})

    recs = list(struct.iter_unpack("<iid", b[off + DATA_RECORD:end]))
    ang = [r[2] for r in recs]
    dwell = recs[0][1]
    if not (MIN_ANGLE <= ang[0] <= SERVICE_ANGLE) or abs(ang[0] - start) > START_TOL:
        raise _Reject("mismatch", {"start": start, "ang0": ang[0]})

    diffs = [ang[i + 1] - ang[i] for i in range(n - 1)]

    # дубликаты 2θ и нарушение порядка — до оценки шага
    for k, d in enumerate(diffs):
        if abs(d) <= 1e-12:
            raise _Reject("duplicates", {"angle": ang[k], "index": k, "start": start})
        if d < 0:
            raise _Reject("unordered", {"index": k, "a": ang[k], "b": ang[k + 1],
                                        "start": start})

    counts = Counter(round(d, 6) for d in diffs)
    step, _cnt = counts.most_common(1)[0]
    if step <= 0:
        raise _Reject("uneven", {"declared": n, "keep": 1, "step": step, "start": start})

    tol = _tol(step)
    keep = 1
    for d in diffs:
        if abs(d - step) <= tol:
            keep += 1
        else:
            break

    if keep < MIN_POINTS:
        raise _Reject("too_few", {"declared": n, "keep": keep, "step": step,
                                  "start": start})

    # ровный участок короче объявленного блока — шаг непостоянен
    if keep < 0.9 * n:
        raise _Reject("uneven", {"declared": n, "keep": keep, "step": step,
                                 "start": start, "broken_at": ang[keep - 1]})

    span = ang[keep - 1] - ang[0]
    if span < MIN_SPAN or span > MAX_SPAN:
        raise _Reject("span", {"declared": n, "keep": keep, "span": span,
                               "start": start})

    # Кусок внутри более длинного ровного ряда — это не начало блока.
    # Признак: блок кончился ровно на объявленном N, а следующий байт-в-байт
    # продолжает тот же шаг и имеет то же время на точку (то есть это данные,
    # а не заголовок следующего блока).
    nxt = off + DATA_RECORD * (keep + 1)
    if keep == n and nxt + DATA_RECORD <= len(b):
        i2, d2, a2 = HEADER.unpack_from(b, nxt)
        if i2 != 0 and d2 == dwell and abs(a2 - (ang[keep - 1] + step)) <= tol:
            raise _Reject("slice", {"declared": n, "keep": keep, "start": start})

    points = [(ang[k], int(recs[k][0])) for k in range(keep)]
    return {"points": points, "step": float(step), "dwell": int(dwell),
            "keep": keep, "declared": int(n), "start": float(start)}


def _scan_blocks(b):
    """Все кандидаты + причины отказа (для понятного сообщения об ошибке).

    Перебираем только позиции внутри runs из нулей: первое поле заголовка —
    int32 0, а побайтовый перебор всего файла на больших сканах заметно
    медленнее. Кандидаты внутри run проверяются все (заголовок может начинаться
    не с первого нуля серии).
    """
    valid = []
    rejects = {"truncated": [], "duplicates": [], "unordered": [],
               "uneven": [], "too_few": [], "span": [], "slice": [],
               "mismatch": []}
    if len(b) < DATA_RECORD * 3:
        return valid, rejects
    for m in _ZERO_RUN.finditer(b):
        hi = min(m.end() - 3, len(b) - DATA_RECORD + 1)
        for off in range(m.start(), max(m.start(), hi)):
            a, n, start = HEADER.unpack_from(b, off)
            if a != 0 or n < 3 or n > MAX_POINTS:
                continue
            if not (MIN_ANGLE <= start <= SERVICE_ANGLE):
                continue
            try:
                valid.append(_validate(b, off, n, start))
            except _Reject as r:
                rejects[r.reason].append(r.detail)
    return valid, rejects


def _explain(name, rejects):
    """Понятный русский текст ошибки по накопленным причинам отказа."""
    def top(reason, key):
        return max(rejects[reason], key=key)

    if rejects["truncated"]:
        d = top("truncated", lambda x: x["declared"])
        return _err(name, (
            "файл обрезан или повреждён — заголовок блока объявляет {dec} точек "
            "(старт {st:.2f}°), но в файле есть только {have} записей данных. "
            "Скопируйте файл с прибора заново.").format(
                dec=d["declared"], st=d["start"], have=d["available"]))
    if rejects["duplicates"]:
        d = top("duplicates", lambda x: x["index"])
        return _err(name, (
            "повреждён блок данных: угол 2θ = {a:.4f}° повторяется дважды "
            "(позиция {i}). Такой файл не может быть корректным сканом.").format(
                a=d["angle"], i=d["index"]))
    if rejects["unordered"]:
        d = top("unordered", lambda x: x["index"])
        return _err(name, (
            "нарушен порядок точек по 2θ: на позиции {i} угол уменьшается "
            "({a:.4f}° -> {b:.4f}°). Ожидается возрастание 2θ.").format(
                i=d["index"], a=d["a"], b=d["b"]))
    if rejects["uneven"]:
        d = top("uneven", lambda x: x["declared"] - x["keep"])
        return _err(name, (
            "неравномерный шаг по 2θ: объявлено {dec} точек с шагом {st:.4f}°, "
            "но ровный участок — только {keep} точек. Файл повреждён.").format(
                dec=d["declared"], st=d["step"], keep=d["keep"]))
    if rejects["too_few"]:
        d = top("too_few", lambda x: x["keep"])
        return _err(name, (
            "слишком мало точек: {keep} (минимум {minp}). Это не полный скан — "
            "файл обрезан или это не дифрактограмма.").format(
                keep=d["keep"], minp=MIN_POINTS))
    if rejects["span"]:
        d = top("span", lambda x: abs(x["span"]))
        return _err(name, (
            "неправдоподобный диапазон скана: {sp:.2f}° по 2θ "
            "(допустимо {lo}..{hi}°).").format(sp=d["span"], lo=MIN_SPAN, hi=MAX_SPAN))
    if rejects["slice"]:
        d = top("slice", lambda x: x["keep"])
        return _err(name, (
            "в файле нет заголовка блока данных: найден только фрагмент "
            "({keep} точек) внутри более длинного ряда. Файл повреждён.").format(
                keep=d["keep"]))
    return _err(name, (
        "не найден блок данных WinXRD: нет записи-заголовка (0, N, 2θ) с "
        "равномерным шагом по 2θ. Это не .raw файл прибора WinXRD."))


def find_block(b):
    """Возвращает (данные, шаг, время_на_точку) или None.

    Оставлено для совместимости; подробные ошибки даёт parse().
    """
    valid, _rejects = _scan_blocks(b)
    if not valid:
        return None
    best = max(valid, key=lambda c: c["keep"])
    return best["points"], best["step"], best["dwell"]


def read_meta(path):
    """Вытаскивает из файла служебные строки прибора (анод, ПО, детектор)."""
    head = Path(path).read_bytes()[:4096]
    meta = {}
    for key, pattern in (("anode", rb"Copper|Cobalt|Molybdenum|Chromium|Iron"),
                         ("software", rb"WinXRD[ \x2d0-9.]*"),
                         ("detector", rb"Solid State \w+|Scintillation"),
                         ("atmosphere", rb"air|vacuum|nitrogen")):
        m = re.search(pattern, head)
        if m:
            meta[key] = m.group().decode("latin-1").strip()
    return meta


def parse(path):
    """-> (данные [(2theta, интенсивность)], шаг, метаданные).

    Бросает ValueError с понятным русским текстом, если файл пустой, чужой
    (картинка/текст), обрезанный, с неравномерным шагом, дубликатами 2θ или
    слишком короткий. Данные не изменяются: возвращаются как записаны.
    """
    p = Path(path)
    b = p.read_bytes()
    name = p.name

    if not b:
        raise _err(name, "файл пустой (0 байт) — это не файл WinXRD .raw")
    if len(b) > MAX_FILE_BYTES:
        raise _err(name, "файл слишком большой ({} МБ) для .raw WinXRD".format(
            len(b) // (1024 * 1024)))
    image = _looks_like_image(b)
    if image:
        raise _err(name, (
            "это файл изображения ({}) — переименованная картинка, а не .raw "
            "WinXRD. Нужен исходный файл прибора.").format(image))
    if _looks_like_text(b):
        raise _err(name, (
            "это текстовый файл, а не бинарный .raw WinXRD (похоже на .xy, "
            ".csv или отчёт). Нужен исходный .raw с прибора."))
    if len(b) < DATA_RECORD * (MIN_POINTS + 1):
        raise _err(name, (
            "файл слишком мал: {} байт, для блока из {} точек нужно минимум {} "
            "байт.").format(len(b), MIN_POINTS, DATA_RECORD * (MIN_POINTS + 1)))

    valid, rejects = _scan_blocks(b)
    if not valid:
        raise _explain(name, rejects)
    best = max(valid, key=lambda c: c["keep"])

    # Обрезанный файл: есть заголовок, который объявляет заведомо больше точек,
    # чем лежит в файле. Иначе можно молча принять короткий фрагмент за скан.
    # Берём только правдоподобные случаи: большая часть объявленного блока
    # действительно есть в файле, а старт близок к найденным данным.
    first = best["points"][0][0]
    trunc = [r for r in rejects["truncated"]
             if r["declared"] > best["keep"] + MIN_POINTS
             and r["available"] >= 0.5 * r["declared"]
             and abs(r["start"] - first) <= 15.0]
    if trunc:
        raise _explain(name, {"truncated": trunc, "duplicates": [],
                              "unordered": [], "uneven": [], "too_few": [],
                              "span": [], "slice": []})

    data = best["points"]
    meta = read_meta(p)
    meta["step"] = best["step"]
    meta["dwell"] = best["dwell"]
    meta["points"] = len(data)
    meta["points_declared"] = best["declared"]
    meta["two_theta_min"] = data[0][0]
    meta["two_theta_max"] = data[-1][0]
    return data, best["step"], meta


# ------------------------------------------------------------------ пересчёты
def d_spacing(two_theta, wavelength=None):
    """Межплоскостное расстояние, Å (закон Вульфа—Брэгга).

    wavelength=None — Cu Kα средневзвешенная из tools/radiation.py (1.54184 Å),
    то есть та же длина волны, что и в peakpick.d_from_2theta и в базе эталонов.
    """
    if wavelength is None:
        wavelength = BASE_WAVELENGTH
    return _rad_d(two_theta, wavelength)


def find_peaks(data, min_intensity=None, top=8):
    """Локальные максимумы; [(2theta, интенсивность)] по убыванию интенсивности."""
    pts = list(data)
    if len(pts) < 5:
        return []
    imax = max(i for _, i in pts)
    if imax <= 0:
        return []
    if min_intensity is None:
        min_intensity = max(40, 0.03 * imax)
    peaks = []
    for k in range(2, len(pts) - 2):
        window = [pts[k + j][1] for j in range(-2, 3)]
        if pts[k][1] != max(window) or pts[k][1] < min_intensity:
            continue
        if peaks and pts[k][0] - peaks[-1][0] < 0.35:
            if pts[k][1] > peaks[-1][1]:
                peaks[-1] = [pts[k][0], pts[k][1]]
        else:
            peaks.append([pts[k][0], pts[k][1]])
    return sorted(peaks, key=lambda p: -p[1])[:top]


def write_xy(data, src, out=None):
    out = Path(out) if out else Path(src).with_suffix(".xy")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        for ang, inten in data:
            f.write("{:8.2f}{:12d}\n".format(float(ang), int(round(inten))))
    return out


def write_csv(data, src, out=None, wavelength=None):
    out = Path(out) if out else Path(src).with_suffix(".csv")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("two_theta_deg,d_A,intensity\n")
        for ang, inten in data:
            f.write("{:.2f},{:.4f},{}\n".format(
                float(ang), d_spacing(ang, wavelength), int(round(inten))))
    return out


def convert_data(path, wavelength=None):
    """Полная конвертация в текстовые форматы. Возвращает (xy, csv, meta)."""
    data, step, meta = parse(path)
    return (write_xy(data, path), write_csv(data, path, wavelength=wavelength), meta)


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
