#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Идентификация фаз и элементов по пикам дифрактограммы .raw.

Связывает всё вместе: .raw -> пики -> база эталонов -> фазы, элементы,
график с подписями, отчёт .txt, листы Excel.
"""

import os
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "tmp"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

import peakpick  # noqa: E402
import phases  # noqa: E402
import plot_xrd  # noqa: E402
import xlsx_out  # noqa: E402
from winxrd2xy import parse  # noqa: E402


def analyze(path, tol=phases.TOL_DEFAULT, min_fom=55.0):
    """Пики + фазы + элементы для файла .raw."""
    data, step, meta = parse(path)
    x = [a for a, _ in data]
    y = [i for _, i in data]
    peaks = peakpick.pick(x, y)
    db = phases.load_db()
    res = phases.explain(peaks, db, tol=tol, min_fom=min_fom)
    return {"peaks": peaks, "result": res, "meta": meta, "data": data,
            "step": step, "db": db}


def label_of(phase):
    """Короткая подпись фазы для графика."""
    name = (phase.get("name") or "").split(" [ICSD")[0].strip()
    if name and not name.replace("-", "").replace("_", "").replace(".", "").isalnum():
        name = ""
    if name and len(name) > 22:
        name = name[:22]
    return name or phase["formula"]


def make_graph(ans, path, theme="dark"):
    """График с подписями пиков и фаз."""
    data = ans["data"]
    peaks = ans["peaks"]
    res = ans["result"]
    t = plot_xrd.THEMES.get(theme, plot_xrd.THEMES["dark"])
    src = Path(path)
    x = [a for a, _ in data]
    y = [i for _, i in data]
    imax = max(y)

    fig = plt.figure(figsize=(12.5, 7.2), dpi=150)
    fig.patch.set_facecolor(t["bg"])
    ax = fig.add_axes([0.075, 0.13, 0.905, 0.75])
    ax.set_facecolor(t["bg"])
    ax.plot(x, y, color=t["curve"], linewidth=1.0, zorder=3)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, imax * 1.30)
    ax.grid(color=t["grid"], linewidth=0.8, zorder=1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
    ax.tick_params(colors=t["dim"], labelsize=10)

    # подписываем пики, у которых есть фаза, либо самые сильные
    idx_map = {id(p): i for i, p in enumerate(peaks)}
    imax_data = max(i for _, i in data)
    shown = 0
    for p in peaks:
        i = idx_map[id(p)]
        ass = res["per_peak"].get(i)
        if not ass and p["intensity"] < 20:
            continue
        if shown > 22:
            break
        if ass:
            ph_no = sorted(ass)[0][0]
            sc = res["phases"][ph_no - 1]
            tag = label_of(sc["phase"])
        else:
            tag = "?"
        yy = p["height"] + imax_data * (0.06 + 0.045 * (shown % 3))
        ax.annotate("{:.2f}°\n{}".format(p["two_theta"], tag),
                    xy=(p["two_theta"], p["height"]), xytext=(p["two_theta"], yy),
                    ha="center", va="bottom", fontsize=8, color=t["peak"], zorder=5,
                    arrowprops=dict(arrowstyle="-", color=t["peak"], linewidth=0.7,
                                    shrinkA=0, shrinkB=2))
        shown += 1

    # легенда фаз
    legend = []
    for i, sc in enumerate(res["phases"], 1):
        legend.append("#{} {} — {}".format(i, label_of(sc["phase"]), sc["phase"]["formula"]))
    if legend:
        fig.text(0.075, 0.945, "  |  ".join(legend[:4]), ha="left", va="top",
                 fontsize=9, color=t["text"])
        if len(legend) > 4:
            fig.text(0.075, 0.917, "  |  ".join(legend[4:]), ha="left", va="top",
                     fontsize=9, color=t["text"])

    title = "{}: идентификация фаз".format(src.stem)
    sub = "2θ {:.2f}–{:.2f}°, шаг {:.2f}°, пиков {}, фаз {}".format(
        x[0], x[-1], ans["step"], len(peaks), len(res["phases"]))
    fig.text(0.5, 0.985, title, ha="center", va="top", fontsize=15, color=t["text"])
    fig.text(0.5, 0.955, sub, ha="center", va="top", fontsize=10, color=t["dim"])
    ax.set_xlabel("2θ, градусы", color=t["text"], fontsize=12)
    ax.set_ylabel("Интенсивность, имп.", color=t["text"], fontsize=12)

    png = src.parent / "{}_phases_{}.png".format(src.stem, theme)
    svg = src.parent / "{}_phases_{}.svg".format(src.stem, theme)
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)
    return png, svg


def save_report(ans, path):
    src = Path(path)
    txt = phases.format_report(ans["result"], ans["peaks"], filename=src.name)
    meta, n_db = phases.db_info(ans["db"])
    txt += "\n\nБаза эталонов: {} фаз ({}), {}".format(
        n_db, meta.get("sources", ""), meta.get("note", ""))
    out = src.parent / "{}_phases.txt".format(src.stem)
    out.write_text(txt, encoding="utf-8-sig")
    return out, txt


def run(path, theme="dark", both_themes=False, tol=phases.TOL_DEFAULT,
        min_fom=55.0, want_xlsx=True):
    """Полный проход: анализ + отчёт + график(и) + Excel."""
    ans = analyze(path, tol=tol, min_fom=min_fom)
    files = []
    themes = ("dark", "light") if both_themes else (theme,)
    for th in themes:
        png, svg = make_graph(ans, path, theme=th)
        files += [png, svg]
    rep, txt = save_report(ans, path)
    files.append(rep)
    if want_xlsx:
        xlsx = xlsx_out.make_xlsx(path, phases_result=ans["result"], peaks=ans["peaks"])
        files.append(xlsx)
    return ans, files, txt


def main():
    """CLI: без меню, для проверки и пакетной работы.

    python tools\\identify.py input\\sample.raw [--theme light] [--both] [--tol 0.3]
    """
    import argparse
    try:
        from console import setup as _setup_console

        _setup_console()
    except Exception:  # noqa: BLE001
        for stream in (sys.stdout, sys.stderr):
            try:
                stream.reconfigure(errors="replace")
            except Exception:  # noqa: BLE001
                pass
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--theme", default="dark", choices=["dark", "light"])
    ap.add_argument("--both", action="store_true")
    ap.add_argument("--tol", type=float, default=phases.TOL_DEFAULT)
    ap.add_argument("--min-fom", type=float, default=55.0)
    a = ap.parse_args()
    for f in a.files:
        p = Path(f)
        if not p.is_file():
            print("нет файла: {}".format(f))
            continue
        print("=" * 60)
        ans, files, txt = run(str(p), theme=a.theme, both_themes=a.both, tol=a.tol,
                              min_fom=a.min_fom)
        print(txt)
        print()
        for x in files:
            print("  {}".format(x))


if __name__ == "__main__":
    main()

