#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Построение графика дифрактограммы из .raw (WinXRD / Thermo ARL X'TRA).

Темы: ``dark`` (экран), ``light``, ``paper`` (светлый научный стиль для
статьи/отчёта: белый фон, тонкая тёмная кривая, светлая сетка).

Верхняя панель — **исходная** экспериментальная кривая (без сглаживания),
пики подписаны по точному углу 2θ. Если передан результат идентификации,
внизу рисуются отдельные полосы Bragg-штрихов по каждой фазе-кандидату
(коды Q, C, G, P, ... с расшифровкой в легенде под панелью, чтобы не
перекрывать кривую).

Сохраняет PNG (>= 300 dpi) и SVG (векторный). Использует matplotlib,
поэтому запускать только встроенным python этой папки.
"""

import bisect
import math
import os
import sys
from pathlib import Path

# кэш шрифтов matplotlib держим внутри папки — иначе на чужом ПК он полезет в профиль
BASE = Path(__file__).resolve().parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "tmp"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from winxrd2xy import find_peaks, parse  # noqa: E402

import result_view  # noqa: E402

export_wh = []      # для тестов/отладки: (png, svg) последнего вызова
LAST = {"labels": [], "codes": [], "n_bragg": 0, "dpi": None}

THEMES = {
    "dark": dict(
        bg="#12151a", panel="#161a20", curve="#4ea1ff", text="#e6edf3",
        dim="#8b949e", grid="#2a2f37", peak="#ffcc66", warn="#ff7b72",
        curve_lw=1.0, grid_lw=0.8, dpi=150,
        bragg=("#ffcc66", "#7ee787", "#d2a8ff", "#79c0ff", "#ffa657",
               "#a5d6ff", "#f778ba", "#c9d1d9"),
    ),
    "light": dict(
        bg="#ffffff", panel="#f6f8fa", curve="#1f4fa8", text="#111418",
        dim="#5a6472", grid="#dfe3e8", peak="#b45309", warn="#a83232",
        curve_lw=1.1, grid_lw=0.8, dpi=300,
        bragg=("#b45309", "#1f7a3f", "#6b3fa0", "#1f4fa8", "#9a3412",
               "#0f766e", "#a21caf", "#374151"),
    ),
    "paper": dict(
        bg="#ffffff", panel="#ffffff", curve="#1a1a1a", text="#111111",
        dim="#4d4d4d", grid="#e6e6e6", peak="#1a1a1a", warn="#a83232",
        curve_lw=0.8, grid_lw=0.6, dpi=300,
        bragg=("#1a1a1a", "#1f4fa8", "#a83232", "#1f7a3f", "#6b3fa0",
               "#8a5a00", "#0f766e", "#555555"),
    ),
}

DEFAULT_THEME = "dark"
MAX_LABELS = 12          # лимит только на КОЛИЧЕСТВО подписей, не на данные
LABEL_TOL = 0.55         # °, ближе этого подписи стараемся не ставить


# --------------------------------------------------------------- излучение
def _coerce_radiation(radiation, meta):
    """Radiation из объекта, dict (rad.as_dict()) или метаданных прибора."""
    try:
        import radiation as _rad
    except Exception:                                   # noqa: BLE001
        return None
    if radiation is not None and hasattr(radiation, "lam") and hasattr(radiation, "describe"):
        return radiation
    if isinstance(radiation, dict) and radiation.get("anode"):
        try:
            return _rad.Radiation.from_settings(meta or {}, anode=radiation.get("anode"),
                                                model=radiation.get("model"))
        except Exception:                               # noqa: BLE001
            pass
    try:
        return _rad.Radiation.from_settings(meta or {})
    except Exception:                                   # noqa: BLE001
        return None


def wavelength_of(rad, default=1.54184):
    return float(getattr(rad, "lam", None) or default)


def radiation_line(rad, short=True):
    """Компактная строка излучения для подписи графика.

    Явно говорит «предположение», если анод не подтверждён.
    """
    if rad is None:
        return "Cu Kα, λ = 1.5418 Å — ПРЕДПОЛОЖЕНИЕ (анод не определён)"
    label = "Kα1" if getattr(rad, "model", "ka_avg") == "ka1" else "Kα"
    src = {"user": "задано пользователем",
           "metadata": "из метаданных прибора",
           "default": "ПРЕДПОЛОЖЕНИЕ, анод не определён"}.get(rad.source, str(rad.source))
    text = "{} {}, λ = {:.4f} Å".format(rad.anode, label, rad.lam)
    if rad.confirmed:
        return text + " — " + src
    return text + " — не подтверждено: " + src


def _ref_positions(score, rad, xmin, xmax, tol=0.3):
    """Позиции рефлексов фазы в текущем излучении внутри диапазона скана."""
    phase = (score or {}).get("phase") or {}
    out = []
    for ref in phase.get("peaks") or ():
        if not isinstance(ref, (list, tuple)) or not ref:
            continue
        try:
            tt = float(ref[0])
        except (TypeError, ValueError):
            continue
        d = ref[1] if len(ref) > 1 else None
        if rad is not None and hasattr(rad, "ref_two_theta") and d:
            try:
                tt = float(rad.ref_two_theta(d, tt))
            except Exception:                           # noqa: BLE001
                pass
        if tt != tt:                                    # nan
            continue
        if xmin - tol <= tt <= xmax + tol:
            i = ref[2] if len(ref) > 2 else None
            out.append((tt, i))
    out.sort(key=lambda p: p[0])
    return out


# ------------------------------------------------------------ подписи пиков
def _curve_height(x, y, tt, window=0.06):
    """Высота кривой в окрестности угла (не сглаженная)."""
    if not x:
        return 0.0
    lo = bisect.bisect_left(x, tt - window)
    hi = bisect.bisect_right(x, tt + window)
    if lo >= hi:
        k = min(range(len(x)), key=lambda i: abs(x[i] - tt))
        return float(y[k])
    return float(max(y[lo:hi]))


def _label_priority(view, peaks, index):
    """(приоритет, -интенсивность): сильные и подтверждённые — первыми."""
    status = view.status_of(index) if view is not None else result_view.STATUS_ASSIGNED
    rank = {result_view.STATUS_ASSIGNED: 0,
            result_view.STATUS_AMBIGUOUS: 1,
            result_view.STATUS_UNEXPLAINED: 2,
            result_view.STATUS_BELOW: 3}.get(status, 4)
    if view is not None:
        score = view.phase((view.best(index) or {}).get("phase_no"))
        if score is not None and result_view.phase_status(score) == "weak":
            rank += 0.5
    return (rank, -float(peaks[index].get("intensity") or 0.0), index)


def select_labels(view, peaks, limit=MAX_LABELS, tol=LABEL_TOL):
    """Индексы пиков, которые стоит подписать (не более ``limit``)."""
    if not peaks:
        return []
    order = sorted(range(len(peaks)), key=lambda i: _label_priority(view, peaks, i))
    chosen = []
    for i in order:
        if len(chosen) >= limit:
            break
        if view is not None and view.status_of(i) == result_view.STATUS_BELOW:
            continue
        tt = float(peaks[i]["two_theta"])
        strong = float(peaks[i].get("intensity") or 0.0) >= 45.0
        if not strong and any(abs(tt - float(peaks[j]["two_theta"])) < tol for j in chosen):
            continue
        chosen.append(i)
    chosen.sort(key=lambda i: float(peaks[i]["two_theta"]))
    return chosen


def _label_color(t, view, index):
    """Цвет подписи: цвет фазы, красноватый для перекрытия, серый для «?»."""
    if view is None:
        return t["peak"]
    status = view.status_of(index)
    if status == result_view.STATUS_AMBIGUOUS:
        return t.get("warn", t["peak"])
    if status == result_view.STATUS_ASSIGNED:
        best = view.best(index) or {}
        no = best.get("phase_no")
        if no:
            return t["bragg"][(int(no) - 1) % len(t["bragg"])]
        return t["peak"]
    return t["dim"]


_RU_LABELS = None


def _ru_labels():
    """Русские названия фаз для профиля стройматериалов (id -> label)."""
    global _RU_LABELS
    if _RU_LABELS is None:
        _RU_LABELS = {}
        try:
            import building_materials as bm
            for rec in bm.PROFILE:
                _RU_LABELS[str(rec.get("id"))] = rec.get("label") or ""
        except Exception:
            pass
    return _RU_LABELS


def phase_display(ph):
    """'CaCO₃ — Кальцит' (русское название, если известно, иначе из базы)."""
    ph = ph or {}
    formula = _formula_pretty(ph.get("formula") or "")
    name = _ru_labels().get(str(ph.get("id"))) or (ph.get("name") or "")
    formula, name = formula.strip(), str(name).strip()
    # имя в базе часто = формуле (BaMg(CO3)2) — не дублируем
    raw = str(ph.get("formula") or "").strip()
    if name.replace(" ", "").lower() == raw.replace(" ", "").lower():
        name = ""
    if formula and name and name.lower() != formula.lower():
        return "{} — {}".format(formula, name)
    return formula or name


def peak_names(view, index, limit=2):
    """Формулы фаз, претендующих на пик (для подписи на графике)."""
    if view is None:
        return []
    names = []
    for entry in view.candidates(index)[:limit]:
        try:
            no = entry.get("phase_no")
        except AttributeError:      # старый кортеж (phase_no, dtheta)
            no = entry[0] if entry else None
        score = view.phase(no) if no else None
        ph = (score or {}).get("phase") or {}
        name = _formula_pretty(ph.get("formula") or "")
        if name and name not in names:
            names.append(name)
    return names


def peak_title(view, index):
    """Подпись материала для пика: 'формула — название' лучшего кандидата.

    Для перекрытия (>=2 фаз) всё равно показываем лучшего кандидата —
    альтернативы перечислены в легенде под графиком, а метка '?*'
    говорит, что пик объясняется неоднозначно.
    """
    if view is None:
        return ""
    best = view.best(index) or {}
    no = best.get("phase_no")
    score = view.phase(no) if no else None
    if not score:
        return ""
    return phase_display((score or {}).get("phase"))


def peak_tag(view, index, codes=None):
    """Метка пика: код фазы, '?*' (перекрытие) или '?' (не объяснён)."""
    if view is None:
        return ""
    status = view.status_of(index)
    if status == result_view.STATUS_AMBIGUOUS:
        return "?*"
    if status == result_view.STATUS_ASSIGNED:
        best = view.best(index) or {}
        no = best.get("phase_no")
        if codes and no and 0 < int(no) <= len(codes):
            return codes[int(no) - 1]
        return "#{}".format(no)
    return "?"


# ------------------------------------------------------------------ рисунок
def _style_axis(ax, t, paper=False):
    ax.set_facecolor(t["bg"])
    ax.grid(color=t["grid"], linewidth=t["grid_lw"], zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["text"] if paper else t["grid"])
        ax.spines[side].set_linewidth(0.8 if paper else 1.0)
    ax.tick_params(colors=t["dim"], labelsize=9.5,
                   direction="out", length=3.5, width=0.8)
    ax.tick_params(which="minor", colors=t["dim"], length=2.0, width=0.6)
    if paper:
        ax.xaxis.set_minor_locator(MultipleLocator(1.0))


def _wrap_items(items, width=150):
    """Разложить элементы легенды по строкам заданной длины."""
    lines, cur = [], ""
    for item in items:
        if cur and len(cur) + len(item) + 3 > width:
            lines.append(cur)
            cur = item
        else:
            cur = item if not cur else cur + "   |   " + item
    if cur:
        lines.append(cur)
    return lines


_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _formula_pretty(formula):
    """CaCO3 -> CaCO₃ (только для подписей на графике)."""
    return str(formula or "").translate(_SUBSCRIPT)


def _plot_note(n_peaks=None, labels_count=None):
    """Пояснение обозначений под легендой (не перекрывает кривую)."""
    return ("подписаны только пики с определённой фазой: формула — название;   "
            "перекрытие (>=2 фаз) — см. легенду\n"
            "Подписи пиков: не более 12 приоритетных (сильные и подтверждённые) "
            "из {}{}; полный список пиков, фаз и отклонений — в .txt и .xlsx".format(
                n_peaks if n_peaks is not None else "всех",
                "" if labels_count is None else " + {} фаз".format(labels_count)))


def make(path, theme=DEFAULT_THEME, outdir=None, with_phases=False,
         phases_result=None, peaks=None, radiation=None):
    """График дифрактограммы.

    path          — .raw (только чтение)
    theme         — 'dark' | 'light' | 'paper'
    outdir        — куда писать (по умолчанию рядом с исходником)
    with_phases   — рисовать фазы (Bragg-штрихи + подписи); False — «чистый» график
    phases_result — результат phases.explain() (структура по docs/CONTRACT.md)
    peaks         — список пиков (индексы должны совпадать с phases_result)
    radiation     — Radiation (или rad.as_dict()); по умолчанию из метаданных
    """
    data, step, meta = parse(path)
    theme = theme if theme in THEMES else DEFAULT_THEME
    t = THEMES[theme]
    paper = theme == "paper"
    src = Path(path)
    outdir = Path(outdir) if outdir else src.parent
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    x = [a for a, _ in data]
    y = [i for _, i in data]
    imax = max(y) if y else 1.0
    if imax <= 0:
        imax = 1.0
    rad = _coerce_radiation(radiation, meta)
    lam = wavelength_of(rad)

    # --- пики: переданные (полный список с индексами) или грубые из кривой ---
    if peaks is None:
        peaks = [{"two_theta": float(a), "intensity": None, "height": float(v)}
                 for a, v in find_peaks(data)]
        for k, p in enumerate(peaks):
            p["intensity"] = 100.0 * p["height"] / imax
            p["n"] = k + 1
    peaks = list(peaks or [])

    view = None
    codes = []
    if with_phases and phases_result is not None:
        view = result_view.PeakView(phases_result, peaks)
        scores = view.phase_scores()
        codes = result_view.phase_codes([sc.get("phase") for sc in scores])
    else:
        scores = []

    # --- Bragg-полосы -------------------------------------------------------
    bands = []
    if view is not None:
        for no, score in enumerate(scores, 1):
            refs = _ref_positions(score, rad, x[0], x[-1])
            if refs:
                bands.append((no, score, refs))
    n_bands = len(bands)

    # --- геометрия ----------------------------------------------------------
    figsize = (12.5, 8.4) if n_bands else (12.5, 6.8)
    fig = plt.figure(figsize=figsize, dpi=t["dpi"])
    fig.patch.set_facecolor(t["bg"])

    n_cols = 1 if n_bands <= 1 else (2 if n_bands == 2 else 3)
    n_rows = (n_bands + n_cols - 1) // n_cols if n_bands else 0
    note_h = 0.028 if n_bands else 0.0
    legend_h = note_h + 0.030 * n_rows

    left, right = 0.078, 0.975
    top = 0.885
    bragg_bottom, bragg_h = 0.0, 0.0
    legend_y = 0.008
    if n_bands:
        bragg_bottom = legend_y + legend_h + 0.080
        bragg_h = min(0.30, 0.048 * n_bands + 0.035)
        main_bottom = bragg_bottom + bragg_h + 0.030
    else:
        main_bottom = 0.125
    ax = fig.add_axes([left, main_bottom, right - left, max(0.30, top - main_bottom)])
    _style_axis(ax, t, paper)
    ax.plot(x, y, color=t["curve"], linewidth=t["curve_lw"], zorder=3,
            solid_joinstyle="round")

    head = imax * (1.38 if view is not None else 1.22)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, head)

    # --- подписи пиков ------------------------------------------------------
    # Подписываем ТОЛЬКО пики с определённой фазой: «формула — название».
    # Пики без фазы остаются без подписи вообще (это просто пик).
    labels = select_labels(view, peaks)
    shown = []
    for k, i in enumerate(labels):
        title = peak_title(view, i)
        if not title:
            continue
        p = peaks[i]
        tt = float(p["two_theta"])
        top_y = _curve_height(x, y, tt)
        tag = peak_tag(view, i, codes)
        prefix = tag if tag and tag not in ("?", "?*") else ""
        text = "{} {}".format(prefix, title) if prefix else title
        yy = top_y + imax * (0.035 + 0.055 * (len(shown) % 3))
        color = _label_color(t, view, i)
        ax.annotate(text, xy=(tt, top_y), xytext=(tt, yy), rotation=90,
                    ha="center", va="bottom", fontsize=8.0 if view is not None else 8.5,
                    color=color, zorder=6, clip_on=False,
                    arrowprops=dict(arrowstyle="-", color=color, linewidth=0.6,
                                    shrinkA=0, shrinkB=1))
        shown.append(tt)

    # --- Bragg-панель -------------------------------------------------------
    if n_bands:
        ax_b = fig.add_axes([left, bragg_bottom, right - left, bragg_h],
                            sharex=ax)
        ax_b.set_facecolor(t["panel"])
        ax_b.set_ylim(0.4, n_bands + 0.6)
        ax_b.set_yticks(list(range(n_bands, 0, -1)))
        ax_b.set_yticklabels([codes[no - 1] for no, _s, _r in bands], fontsize=9.5)
        for side in ("top", "right", "left"):
            ax_b.spines[side].set_visible(False)
        ax_b.spines["bottom"].set_color(t["text"] if paper else t["grid"])
        ax_b.spines["bottom"].set_linewidth(0.8)
        ax_b.tick_params(axis="y", length=0, colors=t["text"])
        ax_b.tick_params(axis="x", colors=t["dim"], labelsize=9.5,
                         direction="out", length=3.5, width=0.8)
        ax_b.grid(axis="x", color=t["grid"], linewidth=t["grid_lw"] * 0.8, zorder=0)
        ax_b.set_axisbelow(True)
        for row, (no, score, refs) in enumerate(bands):
            yy = n_bands - row
            col = t["bragg"][(no - 1) % len(t["bragg"])]
            # тонкий разделитель полос — полосы не сливаются
            ax_b.axhline(yy - 0.5, color=t["grid"], linewidth=0.5, zorder=1)
            for tt, inten in refs:
                # сила рефлекса -> толщина штриха (данные не отбрасываются)
                w = 0.7 + 1.3 * min(1.0, (float(inten) if inten else 50.0) / 100.0)
                ax_b.plot([tt, tt], [yy - 0.34, yy + 0.34], color=col,
                          linewidth=w, solid_capstyle="butt", zorder=3)
        ax_b.set_xlabel("2θ, градусы", color=t["text"], fontsize=11, labelpad=4)
        ax.tick_params(axis="x", labelbottom=False)
        ax_b.set_xlim(x[0], x[-1])
        ax_b.set_ylabel("фазы", color=t["dim"], fontsize=9)

        # --- легенда: цвет штриха = фаза, ниже панели (кривую не перекрывает) ---
        handles = [Line2D([], [], color=t["bragg"][(no - 1) % len(t["bragg"])],
                          linewidth=2.2, solid_capstyle="butt")
                   for no, _s, _r in bands]
        legend_labels = ["{} — {} ({})".format(
            codes[no - 1],
            phase_display(score.get("phase")),
            result_view.confirmation_note(score))
            for no, score, _r in bands]
        leg = fig.legend(handles, legend_labels, loc="lower left",
                         bbox_to_anchor=(left, legend_y + note_h),
                         ncol=n_cols, frameon=False, fontsize=8.4,
                         handlelength=1.5, columnspacing=2.2,
                         borderaxespad=0.0, labelspacing=0.35)
        for txt in leg.get_texts():
            txt.set_color(t["text"])
        fig.text(left, legend_y, _plot_note(len(peaks), labels_count=len(legend_labels)),
                 ha="left", va="bottom", fontsize=7.8, color=t["dim"])
    else:
        ax.set_xlabel("2θ, градусы", color=t["text"], fontsize=11)
        fig.text(left, 0.012,
                 "подписи есть только у пиков с определённой фазой; "
                 "неподписанные пики — фаза не определена (совпадения ниже порога)",
                 ha="left", va="bottom", fontsize=7.8, color=t["dim"])

    ax.set_ylabel("Интенсивность, имп.", color=t["text"], fontsize=11)

    # --- заголовок и легенда ------------------------------------------------
    if view is not None:
        title = "{}: дифрактограмма и фазы-кандидаты".format(src.stem)
    else:
        title = "{}: дифрактограмма".format(src.stem)
    sub1 = "{} | 2θ {:.2f}–{:.2f}°, шаг {:.3f}°, точек {} | максимум {} имп. при 2θ = {:.2f}°".format(
        meta.get("software", "WinXRD"), x[0], x[-1], step, len(data), int(round(imax)),
        max(data, key=lambda p: p[1])[0])
    sub2 = "Излучение: " + radiation_line(rad)
    if view is not None:
        sub2 += "   |   подписано пиков с фазой: {}".format(len(shown))
    fig.text(0.5, 0.975, title, ha="center", va="top", fontsize=14.5, color=t["text"])
    fig.text(0.5, 0.942, sub1, ha="center", va="top", fontsize=9.5, color=t["dim"])
    fig.text(0.5, 0.917, sub2, ha="center", va="top", fontsize=9.0,
             color=t["text"] if paper else t["dim"])

    # --- сохранение ---------------------------------------------------------
    kind = "phases" if (with_phases and phases_result is not None) else "graph"
    png = outdir / "{}_{}_{}.png".format(src.stem, kind, theme)
    svg = outdir / "{}_{}_{}.svg".format(src.stem, kind, theme)
    dpi = 300                     # PNG не ниже 300 dpi (требование к печати)
    fig.savefig(str(png), dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(str(svg), facecolor=fig.get_facecolor())
    plt.close(fig)

    LAST.update({"labels": shown, "codes": list(codes), "n_bragg": n_bands,
                 "dpi": dpi, "kind": kind, "theme": theme,
                 "peaks": len(peaks)})
    export_wh.append((png, svg))
    return png, svg, meta


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    theme = os.environ.get("XRD_THEME", DEFAULT_THEME)
    for arg in sys.argv[1:]:
        png, svg, meta = make(arg, theme=theme)
        print("график: {}".format(png))
        print("вектор: {}".format(svg))


if __name__ == "__main__":
    main()
