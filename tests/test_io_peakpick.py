#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Тесты парсера WinXRD .raw, единого излучения и поиска пиков.

Только синтетика: реальные файлы пользователя не используются. Временные
файлы — в tempfile.mkdtemp(). Запуск:

    python -m unittest discover -s tests -v
"""

import math
import os
import shutil
import struct
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(ROOT, "tools")
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)

import numpy as np  # noqa: E402

import peakpick  # noqa: E402
import radiation  # noqa: E402
import winxrd2xy  # noqa: E402

REL = 1e-9
BASE = radiation.BASE_WAVELENGTH


# --------------------------------------------------------------------- утилиты
def gauss(x, center, fwhm, height):
    sigma = fwhm / 2.3548200450309493
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def asym(x, center, fwhm_l, fwhm_r, height):
    s_l = fwhm_l / 2.3548200450309493
    s_r = fwhm_r / 2.3548200450309493
    return np.where(x < center,
                    height * np.exp(-0.5 * ((x - center) / s_l) ** 2),
                    height * np.exp(-0.5 * ((x - center) / s_r) ** 2))


def build_raw(angles, intensities, dwell=150, preamble=True, service=True,
              header_n=None, header_start=None):
    """Собрать файл в формате WinXRD: заголовок (0, N, старт) + N записей."""
    angles = [float(a) for a in angles]
    intensities = [int(i) for i in intensities]
    step = angles[1] - angles[0] if len(angles) > 1 else 0.02
    out = bytearray()
    if preamble:
        out += b"WinXRD 2.0-6\r\nCopper\r\nSolid State Thermo\r\nair\r\n"
        out += b"1696661561\r\n" + b"\x00" * 16
    n = len(angles) if header_n is None else int(header_n)
    start = (angles[0] - step) if header_start is None else float(header_start)
    out += struct.pack("<iid", 0, n, start)
    for ang, inten in zip(angles, intensities):
        out += struct.pack("<iid", inten, dwell, ang)
    if service:
        # реальный прибор дописывает служебную запись с углом 260°
        out += struct.pack("<iid", 8, dwell, 260.0)
        out += b"\x00" * 32
    return bytes(out)


def make_spectrum(n=2500, start=10.0, step=0.02, base=100.0, noise=0.0,
                  seed=0, peaks=()):
    x = start + step * np.arange(n, dtype=float)
    y = np.full(n, float(base))
    for center, fwhm, height in peaks:
        y += gauss(x, center, fwhm, height)
    if noise:
        y = y + np.random.default_rng(seed).normal(0.0, noise, n)
    return x, y


class TmpDirCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="xrdtest_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, name, data):
        path = os.path.join(self.tmp, name)
        with open(path, "wb") as f:
            f.write(data)
        return path


# --------------------------------------------------------------- единое излучение
class TestRadiationSingleSource(unittest.TestCase):
    def test_constants_come_from_radiation(self):
        self.assertEqual(peakpick.CU_KA, radiation.BASE_WAVELENGTH)
        self.assertEqual(winxrd2xy.LAMBDA_CU_KA1, radiation.ANODES["Cu"]["ka1"])
        self.assertAlmostEqual(peakpick.CU_KA, 1.54184, places=12)
        self.assertAlmostEqual(peakpick.CU_KA1, 1.54056, places=12)

    def test_d_matches_between_modules(self):
        for tt in (4.02, 18.07, 26.66, 29.39, 32.22, 45.0, 56.0):
            a = peakpick.d_from_2theta(tt)
            b = winxrd2xy.d_spacing(tt)
            self.assertLess(abs(a - b), REL, "2θ={}: {} != {}".format(tt, a, b))
            # и то же самое при явной длине волны
            a2 = peakpick.d_from_2theta(tt, 1.54056)
            b2 = winxrd2xy.d_spacing(tt, 1.54056)
            self.assertLess(abs(a2 - b2), REL)

    def test_bragg_reference_quartz(self):
        """Кварц: d = 3.343 Å при Cu Kα (1.5418 Å) -> 2θ ≈ 26.66°."""
        two_theta = peakpick.two_theta_from_d(3.343, 1.5418)
        self.assertAlmostEqual(two_theta, 26.66, delta=0.02)
        d = winxrd2xy.d_spacing(26.66)
        self.assertAlmostEqual(d, 3.343, delta=0.005)
        # обратный пересчёт по той же длине волны — точный
        self.assertLess(abs(peakpick.two_theta_from_d(d) - 26.66), 1e-9)

    def test_ref_two_theta_recalculates_by_d(self):
        """При Kα1 позиция эталона пересчитывается ПО d, а не берётся как есть."""
        rad = radiation.Radiation("Cu", "ka1", source="user", confirmed=True)
        pos_base = radiation.two_theta_from_d(3.343, BASE)
        pos_ka1 = rad.ref_two_theta(3.343)
        self.assertAlmostEqual(pos_base, 26.6578, delta=0.01)
        self.assertAlmostEqual(pos_ka1, 26.6570, delta=0.02)
        self.assertLess(pos_ka1, pos_base)          # короче λ -> меньше угол
        self.assertLess(abs(rad.d(pos_ka1) - 3.343), 1e-9)

    def test_pick_uses_radiation_object(self):
        rad = radiation.Radiation("Cu", "ka1", source="user", confirmed=True)
        pos = rad.ref_two_theta(3.343)
        x, y = make_spectrum(n=800, start=20.0, peaks=[(pos, 0.2, 500.0)])
        peaks = peakpick.pick(x, y, radiation=rad)
        self.assertEqual(len(peaks), 1)
        self.assertAlmostEqual(peaks[0]["wavelength"], 1.54056, places=12)
        self.assertAlmostEqual(peaks[0]["d"], 3.343, delta=0.01)
        # без radiation — базовая длина волны и другой d
        plain = peakpick.pick(x, y)
        self.assertAlmostEqual(plain[0]["wavelength"], BASE, places=12)
        self.assertGreater(plain[0]["d"], peaks[0]["d"])


# ----------------------------------------------------------------------- парсер
class TestParser(TmpDirCase):
    def test_roundtrip(self):
        angles = [4.0 + 0.02 * k for k in range(300)]
        inten = [(k * 7) % 50 + 1 for k in range(300)]
        path = self.write("ok.raw", build_raw(angles, inten))
        data, step, meta = winxrd2xy.parse(path)

        self.assertEqual(len(data), 300)
        self.assertAlmostEqual(step, 0.02, places=12)
        for (ang, i), a0, i0 in zip(data, angles, inten):
            self.assertEqual(ang, a0)          # данные не искажаются
            self.assertEqual(i, i0)
            self.assertIsInstance(i, int)
        self.assertEqual(meta["points"], 300)
        self.assertEqual(meta["points_declared"], 300)
        self.assertEqual(meta["anode"], "Copper")
        self.assertEqual(meta["dwell"], 150)
        self.assertAlmostEqual(meta["two_theta_min"], 4.0, places=12)
        self.assertAlmostEqual(meta["two_theta_max"], 4.0 + 0.02 * 299, places=9)

    def test_service_record_trimmed(self):
        angles = [10.0 + 0.02 * k for k in range(200)]
        inten = [10 + k % 7 for k in range(200)]
        path = self.write("svc.raw", build_raw(angles, inten))
        data, _step, meta = winxrd2xy.parse(path)
        self.assertEqual(len(data), 200)
        self.assertAlmostEqual(data[-1][0], angles[-1], places=9)
        self.assertEqual(meta["points_declared"], 200)

    def test_parse_does_not_modify_input_file(self):
        angles = [10.0 + 0.02 * k for k in range(150)]
        blob = build_raw(angles, [k % 30 + 1 for k in range(150)])
        path = self.write("ro.raw", blob)
        winxrd2xy.parse(path)
        with open(path, "rb") as f:
            self.assertEqual(f.read(), blob)

    def test_empty_file(self):
        path = self.write("empty.raw", b"")
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("пустой", str(ctx.exception))

    def test_image_file(self):
        blob = b"\x89PNG\r\n\x1a\n" + b"\x00\x00\x00\rIHDR" + bytes(range(256)) * 20
        path = self.write("pic.raw", blob)
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("изображения", str(ctx.exception))
        self.assertIn("PNG", str(ctx.exception))

    def test_text_file(self):
        text = "# 2theta intensity\n" + "\n".join(
            "{:.2f} {}".format(10 + 0.02 * k, 100 + k) for k in range(300))
        path = self.write("scan.xy", text.encode("utf-8"))
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("текстовый", str(ctx.exception))

    def test_foreign_binary(self):
        blob = bytes(range(256)) * 16          # нет записи (0, N, угол)
        path = self.write("alien.bin", blob)
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("не найден блок данных", str(ctx.exception))

    def test_truncated_file(self):
        angles = [10.0 + 0.02 * k for k in range(600)]
        inten = [20 + k % 13 for k in range(600)]
        blob = build_raw(angles, inten)
        cut = blob[:int(len(blob) * 0.6)]
        path = self.write("cut.raw", cut)
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        msg = str(ctx.exception)
        self.assertTrue("обрезан" in msg or "поврежд" in msg, msg)
        self.assertIn("600", msg)

    def test_truncated_with_zero_counts(self):
        """Нулевые счётчики дают ложные «заголовки» внутри — всё равно обрезка."""
        angles = [10.0 + 0.02 * k for k in range(600)]
        inten = [0 if k % 17 == 0 else 30 + k % 11 for k in range(600)]
        blob = build_raw(angles, inten)
        path = self.write("cut0.raw", blob[:int(len(blob) * 0.55)])
        with self.assertRaises(ValueError):
            winxrd2xy.parse(path)

    def test_too_few_points(self):
        angles = [10.0 + 0.02 * k for k in range(60)]
        inten = [30 + k for k in range(60)]
        # длинная текстовая шапка, чтобы файл прошёл проверку размера
        blob = build_raw(angles, inten, preamble=True)
        blob = blob.replace(b"1696661561\r\n", b"1696661561\r\n" + b"A" * 2000 + b"\x00" * 8)
        path = self.write("short.raw", blob)
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("мало точек", str(ctx.exception))

    def test_uneven_step(self):
        angles = [10.0 + 0.02 * k for k in range(150)] + \
                 [13.0 + 0.05 * k for k in range(150)]
        inten = [40 + k % 9 for k in range(300)]
        path = self.write("uneven.raw", build_raw(angles, inten))
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("шаг", str(ctx.exception))

    def test_duplicate_angles(self):
        angles = [10.0 + 0.02 * k for k in range(300)]
        angles[150] = angles[149]
        inten = [40 + k % 9 for k in range(300)]
        path = self.write("dup.raw", build_raw(angles, inten))
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("повторяется", str(ctx.exception))

    def test_unordered_angles(self):
        angles = [40.0 - 0.02 * k for k in range(300)]
        inten = [40 + k % 9 for k in range(300)]
        path = self.write("rev.raw", build_raw(angles, inten))
        with self.assertRaises(ValueError) as ctx:
            winxrd2xy.parse(path)
        self.assertIn("порядок", str(ctx.exception))

    def test_missing_file(self):
        with self.assertRaises(OSError):
            winxrd2xy.parse(os.path.join(self.tmp, "nope.raw"))

    def test_find_block_compat(self):
        angles = [10.0 + 0.02 * k for k in range(200)]
        blob = build_raw(angles, [50] * 200)
        found = winxrd2xy.find_block(blob)
        self.assertIsNotNone(found)
        data, step, dwell = found
        self.assertEqual(len(data), 200)
        self.assertAlmostEqual(step, 0.02, places=12)
        self.assertEqual(dwell, 150)

    def test_write_xy_csv(self):
        angles = [10.0 + 0.02 * k for k in range(150)]
        inten = [100 + k for k in range(150)]
        path = self.write("w.raw", build_raw(angles, inten))
        data, _step, _meta = winxrd2xy.parse(path)
        xy = winxrd2xy.write_xy(data, path, os.path.join(self.tmp, "w.xy"))
        csv = winxrd2xy.write_csv(data, path, os.path.join(self.tmp, "w.csv"))
        with open(xy, encoding="utf-8") as f:
            lines = f.read().splitlines()
        self.assertEqual(len(lines), 150)
        self.assertEqual(lines[0].split(), ["10.00", "100"])
        with open(csv, encoding="utf-8") as f:
            rows = f.read().splitlines()
        self.assertEqual(rows[0], "two_theta_deg,d_A,intensity")
        self.assertEqual(len(rows), 151)
        tt, d, inten_s = rows[1].split(",")
        self.assertAlmostEqual(float(d), winxrd2xy.d_spacing(float(tt)), places=4)
        self.assertEqual(inten_s, "100")

    def test_convert_data(self):
        angles = [10.0 + 0.02 * k for k in range(150)]
        path = self.write("c.raw", build_raw(angles, [70] * 150))
        xy, csv, meta = winxrd2xy.convert_data(path)
        self.assertTrue(os.path.isfile(xy))
        self.assertTrue(os.path.isfile(csv))
        self.assertEqual(meta["points"], 150)


# ------------------------------------------------------------------- пик-пикинг
class TestPeakPick(unittest.TestCase):
    def test_noise_gives_no_peaks(self):
        for seed in range(5):
            x, y = make_spectrum(n=2500, base=100.0, noise=5.0, seed=seed)
            peaks = peakpick.pick(x, y)
            self.assertEqual(peaks, [], "seed {}: {}".format(
                seed, [round(p["two_theta"], 2) for p in peaks]))

    def test_noise_sigma_estimate(self):
        x, y = make_spectrum(n=3000, base=100.0, noise=7.0, seed=11)
        self.assertAlmostEqual(peakpick.noise_sigma(y), 7.0, delta=0.6)

    def test_adaptive_thresholds_are_the_fix(self):
        x, y = make_spectrum(n=2500, base=100.0, noise=5.0, seed=3)
        self.assertEqual(len(peakpick.pick(x, y)), 0)
        old_style = peakpick.pick(x, y, adaptive=False)
        self.assertGreater(len(old_style), 10)   # так вёл себя старый алгоритм

    def test_flat_signal(self):
        x, y = make_spectrum(n=2000, base=50.0)
        self.assertEqual(peakpick.pick(x, y), [])

    def test_sloping_background(self):
        x = 10.0 + 0.02 * np.arange(2500)
        y = np.linspace(100.0, 300.0, x.size) + \
            np.random.default_rng(3).normal(0, 3.0, x.size)
        self.assertEqual(peakpick.pick(x, y), [])

    def test_single_peak_position_d_fwhm(self):
        x, y = make_spectrum(n=2000, base=100.0, noise=5.0, seed=1,
                             peaks=[(25.0, 0.2, 500.0)])
        peaks = peakpick.pick(x, y)
        self.assertEqual(len(peaks), 1)
        p = peaks[0]
        self.assertAlmostEqual(p["two_theta"], 25.0, delta=0.02)
        self.assertAlmostEqual(p["fwhm"], 0.2, delta=0.03)
        self.assertAlmostEqual(p["intensity"], 100.0, places=6)
        self.assertAlmostEqual(p["d"], peakpick.d_from_2theta(p["two_theta"]),
                               places=12)
        self.assertGreater(p["snr"], 20.0)

    def test_wide_asymmetric_peak_fwhm(self):
        x = 10.0 + 0.02 * np.arange(2500)
        y = 100.0 + asym(x, 25.0, 0.5, 1.0, 500.0) + \
            np.random.default_rng(4).normal(0, 5.0, x.size)
        peaks = peakpick.pick(x, y)
        self.assertEqual(len(peaks), 1)
        # истинная FWHM асимметричного профиля = 0.25 + 0.50 = 0.75°
        self.assertAlmostEqual(peaks[0]["fwhm"], 0.75, delta=0.12)
        self.assertAlmostEqual(peaks[0]["two_theta"], 25.0, delta=0.05)

    def test_two_close_peaks_not_merged(self):
        x, y = make_spectrum(n=2000, base=100.0, noise=5.0, seed=5,
                             peaks=[(20.0, 0.12, 400.0), (20.35, 0.12, 300.0)])
        peaks = peakpick.pick(x, y)
        self.assertEqual(len(peaks), 2, [round(p["two_theta"], 3) for p in peaks])
        self.assertAlmostEqual(peaks[0]["two_theta"], 20.0, delta=0.03)
        self.assertAlmostEqual(peaks[1]["two_theta"], 20.35, delta=0.03)

    def test_weak_line_next_to_strong(self):
        x, y = make_spectrum(n=2000, base=100.0, noise=5.0, seed=7,
                             peaks=[(20.0, 0.15, 500.0), (22.0, 0.2, 60.0)])
        peaks = peakpick.pick(x, y)
        self.assertEqual(len(peaks), 2, [round(p["two_theta"], 2) for p in peaks])
        weak = [p for p in peaks if abs(p["two_theta"] - 22.0) < 0.1]
        self.assertEqual(len(weak), 1)
        self.assertAlmostEqual(weak[0]["intensity"], 12.0, delta=4.0)
        self.assertGreater(weak[0]["snr"], 8.0)

    def test_edge_points_are_not_peaks(self):
        x, y = make_spectrum(n=1200, base=100.0, noise=4.0, seed=9)
        y[0] = 900.0
        y[-1] = 900.0
        peaks = peakpick.pick(x, y)
        self.assertEqual(peaks, [])

    def test_peak_fields_and_order(self):
        x, y = make_spectrum(n=2500, base=100.0, noise=4.0, seed=2,
                             peaks=[(15.0, 0.2, 300.0), (25.0, 0.25, 600.0),
                                    (40.0, 0.3, 200.0)])
        peaks = peakpick.pick(x, y)
        self.assertEqual(len(peaks), 3)
        self.assertEqual([p["n"] for p in peaks], [1, 2, 3])
        self.assertEqual([p["peak_id"] for p in peaks], [0, 1, 2])
        self.assertEqual(peaks, sorted(peaks, key=lambda p: p["two_theta"]))
        for p in peaks:
            self.assertIn("wavelength", p)
            self.assertAlmostEqual(p["wavelength"], BASE, places=12)
            self.assertGreater(p["snr"], 0)
            self.assertIsNotNone(p["fwhm"])
        self.assertAlmostEqual(max(p["intensity"] for p in peaks), 100.0, places=6)

    def test_max_peaks(self):
        pk = [(12.0 + 2.0 * k, 0.2, 500.0 - 20.0 * k) for k in range(15)]
        x, y = make_spectrum(n=3000, base=100.0, noise=3.0, seed=4, peaks=pk)
        peaks = peakpick.pick(x, y, max_peaks=5)
        self.assertEqual(len(peaks), 5)
        self.assertEqual([p["peak_id"] for p in peaks], [0, 1, 2, 3, 4])

    def test_max_peaks_keeps_strongest(self):
        pk = [(12.0, 0.2, 500.0), (20.0, 0.2, 900.0), (30.0, 0.2, 300.0)]
        x, y = make_spectrum(n=2500, base=100.0, noise=3.0, seed=6, peaks=pk)
        peaks = peakpick.pick(x, y, max_peaks=2)
        self.assertEqual(len(peaks), 2)
        self.assertAlmostEqual(peaks[0]["two_theta"], 12.0, delta=0.05)
        self.assertAlmostEqual(peaks[1]["two_theta"], 20.0, delta=0.05)

    def test_old_kwargs_still_work(self):
        x, y = make_spectrum(n=2000, base=100.0, noise=4.0, seed=8,
                             peaks=[(20.0, 0.2, 500.0), (30.0, 0.2, 100.0)])
        self.assertEqual(len(peakpick.pick(x, y)), 2)
        strict = peakpick.pick(x, y, min_prominence=0.9, min_height=0.5)
        self.assertEqual(len(strict), 1)
        self.assertAlmostEqual(strict[0]["two_theta"], 20.0, delta=0.05)
        # позиционные аргументы в прежнем порядке
        pos = peakpick.pick(x, y, 0.02, 0.01, 7, BASE, 120)
        self.assertEqual(len(pos), 2)

    def test_short_input(self):
        self.assertEqual(peakpick.pick([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]), [])
        self.assertEqual(peakpick.pick([], []), [])

    def test_length_mismatch(self):
        with self.assertRaises(ValueError):
            peakpick.pick([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
                          [1.0, 2.0, 3.0])

    def test_input_arrays_not_modified(self):
        x, y = make_spectrum(n=1500, base=100.0, noise=4.0, seed=12,
                             peaks=[(25.0, 0.2, 400.0)])
        y0 = y.copy()
        x0 = x.copy()
        peakpick.pick(x, y)
        np.testing.assert_array_equal(y, y0)
        np.testing.assert_array_equal(x, x0)

    def test_smooth_does_not_mutate(self):
        y = np.array([1.0, 5.0, 3.0, 8.0, 2.0, 4.0, 6.0, 7.0, 1.0])
        y0 = y.copy()
        peakpick.smooth(y, 5)
        np.testing.assert_array_equal(y, y0)


if __name__ == "__main__":
    unittest.main()
