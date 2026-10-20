#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Генератор тестового .raw в формате WinXRD (Thermo ARL X'TRA).

Собирает синтетическую дифрактограмму из эталонных фаз базы + шум и пишет
файл в том же формате, что и прибор: заголовок + блок записей по 16 байт.
Нужен для проверки конвейера, когда под рукой нет реального образца.

  ./codenv/bin/python make_test_raw.py кальцит кварц -o test_mix.raw
"""
import argparse
import gzip
import json
import math
import struct
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent


def load_db():
    with gzip.open(str(BASE / "phases_db.json.gz"), "rt", encoding="utf-8") as f:
        return json.load(f)


def find_phase(db, term):
    term = term.lower()
    for p in db["phases"]:
        name = (p.get("name") or "").lower()
        if term in name or term == p["formula"].lower():
            return p
    return None


def synth(phases, start=5.0, stop=90.0, step=0.02, noise=1.0, seed=7):
    import random

    rnd = random.Random(seed)
    n = int(round((stop - start) / step)) + 1
    y = [30.0] * n
    for ph in phases:
        for tt, d, inten in ph["peaks"]:
            if not (start < tt < stop):
                continue
            fwhm = 0.09 + 0.0008 * tt
            sigma = fwhm / 2.355
            k0 = int(round((tt - start) / step))
            for k in range(max(0, k0 - 40), min(n, k0 + 41)):
                x = start + k * step
                y[k] += inten * 4.0 * math.exp(-0.5 * ((x - tt) / sigma) ** 2)
    y = [max(0.0, v + rnd.gauss(0, noise)) for v in y]
    return start, step, [int(round(v * 10)) for v in y]


def write_raw(path, start, step, ints, dwell=150):
    head = (b"Test sample\r\nThermo ARL X'TRA\r\nWinXRD 2.0-6\r\nCopper\r\n"
            b"Solid State Thermo\r\nair\r\nOperator: synthetic\r\n")
    data = bytearray(head)
    data += struct.pack("<iid", 0, len(ints), start)
    for k, v in enumerate(ints):
        data += struct.pack("<iid", v, dwell, start + k * step)
    Path(path).write_bytes(bytes(data))
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phases", nargs="+", help="названия/формулы фаз из базы")
    ap.add_argument("-o", "--out", default="test_mix.raw")
    ap.add_argument("--noise", type=float, default=1.0)
    a = ap.parse_args()
    db = load_db()
    chosen = []
    for t in a.phases:
        p = find_phase(db, t)
        if p is None:
            print("нет фазы: {}".format(t))
            continue
        chosen.append(p)
        print("фаза: {} {} {}".format(p.get("name", "")[:40], p["formula"], p.get("sg", "")))
    if not chosen:
        sys.exit("ни одной фазы не найдено")
    start, step, ints = synth(chosen, noise=a.noise)
    write_raw(a.out, start, step, ints)
    print("записано: {} ({} точек, {}–{}°)".format(a.out, len(ints), start,
                                                   round(start + (len(ints) - 1) * step, 2)))


if __name__ == "__main__":
    main()
