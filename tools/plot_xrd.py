#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Построение графика дифрактограммы из .raw (WinXRD / Thermo ARL X'TRA).

Сохраняет PNG (для экрана, тёмная тема) и SVG (векторный, для Word/диссертации).
Использует matplotlib, поэтому запускать только встроенным python этой папки.
"""

import os
import sys
from pathlib import Path

# кэш шрифтов matplotlib держим внутри папки — иначе на чужом ПК он полезет в профиль
BASE = Path(__file__).resolve().parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "tmp"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from winxrd2xy import d_spacing, find_peaks, parse  # noqa: E402

THEMES = {
    "dark": dict(bg="#12151a", panel="#161a20", curve="#4ea1ff", text="#e6edf3",
                 dim="#8b949e", grid="#2a2f37", peak="#ffcc66"),
    "light": dict(bg="#ffffff", panel="#f6f8fa", curve="#1f4fa8", text="#111418",
                  dim="#5a6472", grid="#dfe3e8", peak="#b45309"),
}


def make(path, theme="dark", outdir=None):
    data, step, meta = parse(path)
    t = THEMES.get(theme, THEMES["dark"])
    src = Path(path)
    outdir = Path(outdir) if outdir else src.parent
    name = src.stem

    x = [a for a, _ in data]
    y = [i for _, i in data]
    imax = max(y)
    top2 = max(data, key=lambda p: p[1])
    peaks = find_peaks(data)

    fig = plt.figure(figsize=(12.5, 6.8), dpi=150)
    fig.patch.set_facecolor(t["bg"])
    ax = fig.add_axes([0.075, 0.13, 0.905, 0.76])
    ax.set_facecolor(t["bg"])

    ax.plot(x, y, color=t["curve"], linewidth=1.1, zorder=3)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, imax * 1.16)
    ax.grid(color=t["grid"], linewidth=0.8, zorder=1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
    ax.tick_params(colors=t["dim"], labelsize=10)

    # подписи главных пиков: угол и межплоскостное расстояние
    for k, (ang, inten) in enumerate(peaks):
        ax.annotate("{:.2f}°\nd={:.2f} Å".format(ang, d_spacing(ang)),
                    xy=(ang, inten), xytext=(ang, inten + imax * (0.075 if k % 2 == 0 else 0.145)),
                    ha="center", va="bottom", fontsize=9, color=t["peak"], zorder=4,
                    arrowprops=dict(arrowstyle="-", color=t["peak"], linewidth=0.8,
                                    shrinkA=0, shrinkB=2))

    anode = meta.get("anode", "Cu")
    soft = meta.get("software", "WinXRD")
    title = "{}: дифрактограмма, {} Kα (λ = 1.5406 Å)".format(name, anode)
    sub = "{} | 2θ {:.2f}–{:.2f}°, шаг {:.2f}°, точек {} | максимум {} имп. при 2θ = {:.2f}°".format(
        soft, x[0], x[-1], step, len(data), imax, top2[0])
    fig.text(0.5, 0.965, title, ha="center", va="top", fontsize=15, color=t["text"])
    fig.text(0.5, 0.912, sub, ha="center", va="top", fontsize=10, color=t["dim"])
    ax.set_xlabel("2θ, градусы", color=t["text"], fontsize=12)
    ax.set_ylabel("Интенсивность, имп.", color=t["text"], fontsize=12)

    png = outdir / "{}_graph_{}.png".format(name, theme)
    svg = outdir / "{}_graph_{}.svg".format(name, theme)
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)
    return png, svg, meta


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    for arg in sys.argv[1:]:
        png, svg, meta = make(arg)
        print("график: {}".format(png))
        print("вектор: {}".format(svg))


if __name__ == "__main__":
    main()
