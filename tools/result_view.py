#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Единое чтение структуры результата ``phases.explain()``.

Модуль ничего не считает заново: он только приводит результат к одному виду,
в котором **все индексы — индексы полного списка пиков** (docs/CONTRACT.md,
п. 3–4). Пользуются им ``identify.py``, ``plot_xrd.py`` и ``xlsx_out.py``,
чтобы график, TXT и XLSX не расходились между собой.

Поддерживаются обе версии структуры:

* новая (контракт 1.0): ключи ``per_peak`` / ``unassigned`` / ``ambiguous`` /
  ``peak_status`` в пространстве полного списка, значения ``per_peak`` —
  словари ``{phase_no, dtheta, ref_int, overlap}``;
* старая: индексы идут в пространстве отфильтрованного списка (пики с
  ``intensity >= MIN_OBS_INT``), значения ``per_peak`` — кортежи
  ``(phase_no, dtheta)``.

Отсутствие новых ключей не ошибка: определяется «пространство» индексов,
и при необходимости индексы пересчитываются по той же маске фильтра, которую
применяет ``phases.explain``.

Зависимостей от numpy/matplotlib/openpyxl нет — только стандартная библиотека.
"""

STATUS_ASSIGNED = "assigned"
STATUS_AMBIGUOUS = "ambiguous"
STATUS_UNEXPLAINED = "unexplained"
STATUS_BELOW = "below_threshold"

STATUSES = (STATUS_ASSIGNED, STATUS_AMBIGUOUS, STATUS_UNEXPLAINED, STATUS_BELOW)

STATUS_RU = {
    STATUS_ASSIGNED: "совпал",
    STATUS_AMBIGUOUS: "перекрытие (>=2 фаз)",
    STATUS_UNEXPLAINED: "не объяснён",
    STATUS_BELOW: "под порогом поиска",
}

# ключи, которые появились вместе с контрактом 1.0 — их наличие означает,
# что индексы уже в пространстве полного списка. Только те, что однозначно
# появились вместе с новой структурой (settings/radiation могут быть добавлены
# совместимым кодом и о пространстве индексов ничего не говорят).
_CONTRACT_KEYS = ("peak_status", "ambiguous")

# явная метка пространства индексов, которую может поставить потребитель
# (identify.analyze ставит её, когда phases.explain вернул старую структуру)
SPACE_KEY = "index_space"

DEFAULT_MIN_OBS_INT = 3.0

# Короткие коды фаз: устойчивая расшифровка по минералогическому имени.
# Порядок важен — проверяется вхождение подстроки, первое совпадение выигрывает.
MINERAL_CODES = (
    ("portlandite", "P"),
    ("quartz", "Q"),
    ("calcite", "C"),
    ("gypsum", "G"),
    ("bassanite", "Bs"),
    ("anhydrite", "An"),
    ("ettringite", "E"),
    ("thaumasite", "Th"),
    ("afwillite", "Af"),
    ("tobermorite", "T"),
    ("hydrogarnet", "Hg"),
    ("hydrotalcite", "Ht"),
    ("dolomite", "D"),
    ("magnesite", "Mg"),
    ("aragonite", "Ar"),
    ("vaterite", "V"),
    ("brucite", "B"),
    ("halite", "Ha"),
    ("syngenite", "Sy"),
    ("thenardite", "Td"),
    ("mirabilite", "Mi"),
    ("hematite", "H"),
    ("magnetite", "Mt"),
    ("goethite", "Go"),
    ("anatase", "A"),
    ("rutile", "R"),
    ("mullite", "M"),
    ("wollastonite", "W"),
    ("lime", "L"),
    ("periclase", "Pe"),
    ("albite", "Ab"),
    ("microcline", "Mc"),
    ("orthoclase", "Or"),
    ("anorthite", "Ah"),
    ("kaolinite", "K"),
    ("clinochlore", "Ch"),
    ("muscovite", "Ms"),
    ("amorphous", "Am"),
)


# ------------------------------------------------------------------ мелочи
def _intensity(peak):
    try:
        return float(peak.get("intensity") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _min_obs_int(settings=None):
    """Порог интенсивности, ниже которого пики не участвуют в поиске фаз."""
    settings = settings or {}
    for key in ("min_obs_int", "min_intensity"):
        if settings.get(key) is not None:
            try:
                return float(settings[key])
            except (TypeError, ValueError):
                pass
    try:                                            # как в phases.py
        import phases as _phases

        return float(getattr(_phases, "MIN_OBS_INT", DEFAULT_MIN_OBS_INT))
    except Exception:                               # noqa: BLE001 — модуль недоступен
        return DEFAULT_MIN_OBS_INT


def normalize_entries(value):
    """Привести значение ``per_peak[i]`` к списку словарей.

    Понимает и новые словари, и старые кортежи ``(phase_no, dtheta)``.
    """
    out = []
    for item in value or ():
        if isinstance(item, dict):
            phase_no = _as_int(item.get("phase_no", item.get("phase")))
            if phase_no is None:
                continue
            try:
                dtheta = float(item.get("dtheta") or 0.0)
            except (TypeError, ValueError):
                dtheta = 0.0
            out.append({"phase_no": phase_no, "dtheta": dtheta,
                        "ref_int": item.get("ref_int"),
                        "overlap": item.get("overlap")})
        elif isinstance(item, (tuple, list)) and item:
            phase_no = _as_int(item[0])
            if phase_no is None:
                continue
            try:
                dtheta = float(item[1]) if len(item) > 1 and item[1] is not None else 0.0
            except (TypeError, ValueError):
                dtheta = 0.0
            out.append({"phase_no": phase_no, "dtheta": dtheta,
                        "ref_int": None, "overlap": None})
    out.sort(key=lambda e: (e["phase_no"], abs(e["dtheta"])))
    return out


# ------------------------------------------------------------------ фазы
def phase_status(score):
    """Статус фазы из контракта (candidate/ambiguous/weak) с разумным дефолтом."""
    status = (score or {}).get("status")
    return status if status in ("candidate", "ambiguous", "weak") else "candidate"


def phase_flags(score):
    flags = (score or {}).get("flags") or []
    return [str(f) for f in flags]


def phase_prior(score):
    prior = (score or {}).get("prior") or []
    return [str(p) for p in prior]


def phase_datasource(score):
    """'COD' | 'JARVIS-DFT' (берём из результата, иначе из записи базы)."""
    score = score or {}
    src = score.get("datasource")
    if src:
        return str(src)
    phase = score.get("phase") or {}
    if str(phase.get("src", "")).upper().startswith("JARVIS"):
        return "JARVIS-DFT"
    if phase.get("src"):
        return str(phase["src"])
    return ""


def is_dft(score):
    return "JARVIS" in phase_datasource(score).upper()


def phase_name(phase, limit=0):
    """Читаемое имя фазы без инвентарного номера; при отсутствии — формула."""
    phase = phase or {}
    name = str(phase.get("name") or "").split(" [ICSD")[0].strip()
    if name and not any(ch.isalpha() for ch in name):
        name = ""
    formula = str(phase.get("formula") or "?")
    text = name or formula
    if limit and len(text) > limit:
        text = text[: max(1, limit - 1)].rstrip() + "…"
    return text


def confirmation_note(score):
    """Короткая пометка о строгости идентификации (для подписей и легенды)."""
    notes = []
    status = phase_status(score)
    if status == "ambiguous":
        notes.append("перекрытие")
    elif status == "weak":
        notes.append("слабые линии")
    if is_dft(score):
        notes.append("DFT")
    if "FEW_DIAGNOSTIC_LINES" in phase_flags(score):
        notes.append("мало линий в диапазоне")
    notes.append("кандидат")            # без подтверждения пользователем всегда кандидат
    return ", ".join(notes)


def phase_label(score, limit=24, with_note=True):
    """Подпись фазы: имя + пометка «кандидат»/«перекрытие»."""
    score = score or {}
    label = phase_name(score.get("phase"), limit=limit)
    if with_note:
        return "{} ({})".format(label, confirmation_note(score))
    return label


def _code_from_text(text):
    words = [w for w in str(text or "").replace("-", " ").replace("_", " ").split()
             if w[:1].isalpha()]
    if len(words) >= 2:
        return "".join(w[0].upper() for w in words[:3])
    if words:
        return words[0][:2].capitalize()
    return "?"


def phase_codes(phase_dicts):
    """Устойчивые короткие коды фаз: Q, C, G, P, ... (уникальные в наборе).

    Сначала — минералогическая таблица (кварц всегда Q, кальцит всегда C),
    иначе инициалы имени/формулы. Совпадения разводятся цифрой.
    """
    codes, used = [], set()
    for phase in phase_dicts or []:
        phase = phase or {}
        name = str(phase.get("name") or "").lower()
        code = None
        for key, hint in MINERAL_CODES:
            if key in name:
                code = hint
                break
        if not code:
            code = _code_from_text(phase.get("name")) if phase.get("name") else None
        if not code or code == "?":
            formula = str(phase.get("formula") or "")
            code = "".join(ch for ch in formula if ch.isalpha())[:2].capitalize() or "?"
        base, k = code, 2
        while code in used:
            code = "{}{}".format(base, k)
            k += 1
        used.add(code)
        codes.append(code)
    return codes


# ------------------------------------------------------------------ основное
class PeakView(object):
    """Результат ``phases.explain`` в терминах полного списка пиков.

    Атрибуты:
      n          — сколько пиков было передано в explain
      legacy     — True, если результат пришёл в старой структуре
      space      — 'full' | 'filtered', пространство индексов в исходном dict
      per_peak   — {индекс полного списка: [ {phase_no, dtheta, ...}, ... ]}
      status     — {индекс: assigned|ambiguous|unexplained|below_threshold}
      assigned / ambiguous / unexplained / below — списки индексов
    """

    def __init__(self, result=None, peaks=None):
        self.result = result if isinstance(result, dict) else {}
        self.peaks = list(peaks or [])
        self.n = len(self.peaks)
        self.settings = dict(self.result.get("settings") or {})
        self.min_obs_int = _min_obs_int(self.settings)
        self._filtered = None
        self.space = self._detect_space()
        self.legacy = self.space == "filtered"
        self.per_peak = {}
        self.status = {}
        self.assigned = []
        self.ambiguous = []
        self.unexplained = []
        self.below = []
        self._build()

    # ---------------------------------------------------------- пространство
    def _filtered_indices(self):
        """Индексы пиков, которые phases.explain оставляет для поиска фаз."""
        if self._filtered is None:
            idx = [i for i, p in enumerate(self.peaks)
                   if _intensity(p) >= self.min_obs_int]
            if not idx and self.peaks:
                idx = list(range(self.n))
            self._filtered = idx
        return self._filtered

    def _legacy_shape(self):
        """Следы старой структуры: ``per_peak`` с кортежами ``(phase_no, dtheta)``.

        Контракт 1.0 (§4) требует словари ``{phase_no, dtheta, ref_int, overlap}``;
        кортежи мог выдать только старый ``phases.explain``, а он нумеровал пики
        в отфильтрованном списке (§3). Это однозначное свидетельство старого
        формата, даже если маска фильтра совпала с полным списком.
        """
        saw_tuple = False
        for value in (self.result.get("per_peak") or {}).values():
            for item in value or ():
                if isinstance(item, dict):
                    return False            # новый формат — не старый
                if isinstance(item, (tuple, list)):
                    saw_tuple = True
        return saw_tuple

    def _detect_space(self):
        declared = self.result.get(SPACE_KEY)
        if declared in ("full", "filtered"):
            return declared
        if any(key in self.result for key in _CONTRACT_KEYS):
            return "full"
        if self._legacy_shape():
            return "filtered"
        keys = list((self.result.get("per_peak") or {}).keys())
        keys += list(self.result.get("unassigned") or ())
        ints = [k for k in (_as_int(k) for k in keys) if k is not None]
        filtered = self._filtered_indices()
        # старый код не мог выдать индекс >= len(filtered)
        if ints and max(ints) >= len(filtered) and self.n > len(filtered):
            return "full"
        if len(filtered) == self.n:
            return "full"
        return "filtered"

    def to_full(self, index):
        """Индекс в исходном dict -> индекс полного списка пиков (или None)."""
        index = _as_int(index)
        if index is None:
            return None
        if self.space == "full":
            return index if 0 <= index < self.n else None
        filtered = self._filtered_indices()
        if 0 <= index < len(filtered):
            return filtered[index]
        if 0 <= index < self.n:
            return index
        return None

    # ---------------------------------------------------------------- сборка
    def _build(self):
        self.space = self._detect_space()

        # пик -> фазы
        per_peak = {}
        for key, value in (self.result.get("per_peak") or {}).items():
            full = self.to_full(key)
            if full is None:
                continue
            per_peak.setdefault(full, []).extend(normalize_entries(value))
        for full, entries in list(per_peak.items()):
            seen, unique = set(), []
            for entry in entries:
                if entry["phase_no"] in seen:
                    continue
                seen.add(entry["phase_no"])
                unique.append(entry)
            unique.sort(key=lambda e: (abs(e["dtheta"]), e["phase_no"]))
            per_peak[full] = unique
        self.per_peak = per_peak

        declared_unassigned = set()
        for key in self.result.get("unassigned") or ():
            full = self.to_full(key)
            if full is not None:
                declared_unassigned.add(full)

        observed = None
        if self.space == "filtered":
            observed = set(self._filtered_indices())

        given = {}
        raw_status = self.result.get("peak_status") or {}
        if isinstance(raw_status, dict):
            for key, value in raw_status.items():
                full = self.to_full(key)
                if full is None or value not in STATUSES:
                    continue
                given[full] = value

        declared_ambiguous = set()
        for key in self.result.get("ambiguous") or ():
            full = self.to_full(key)
            if full is not None:
                declared_ambiguous.add(full)

        for i in range(self.n):
            status = given.get(i)
            if status is None:
                entries = self.per_peak.get(i) or []
                if i in declared_ambiguous or len({e["phase_no"] for e in entries}) >= 2:
                    status = STATUS_AMBIGUOUS
                elif entries:
                    status = STATUS_ASSIGNED
                elif observed is not None and i not in observed:
                    status = STATUS_BELOW
                elif i in declared_unassigned:
                    status = STATUS_UNEXPLAINED
                elif _intensity(self.peaks[i]) < self.min_obs_int:
                    status = STATUS_BELOW
                else:
                    status = STATUS_UNEXPLAINED
            self.status[i] = status

        self.assigned = [i for i in range(self.n) if self.status[i] == STATUS_ASSIGNED]
        self.ambiguous = [i for i in range(self.n) if self.status[i] == STATUS_AMBIGUOUS]
        self.unexplained = [i for i in range(self.n) if self.status[i] == STATUS_UNEXPLAINED]
        self.below = [i for i in range(self.n) if self.status[i] == STATUS_BELOW]

    # -------------------------------------------------------------- запросы
    def status_of(self, index):
        return self.status.get(index, STATUS_UNEXPLAINED)

    def phases_of(self, index):
        """Номера фаз (1-based, как в res['phases']), претендующих на пик."""
        return [e["phase_no"] for e in self.per_peak.get(index, ())]

    def candidates(self, index):
        return list(self.per_peak.get(index, ()))

    def best(self, index):
        entries = self.per_peak.get(index) or []
        return entries[0] if entries else None

    def dtheta_of(self, index):
        best = self.best(index)
        return None if best is None else best["dtheta"]

    def is_ambiguous(self, index):
        return self.status_of(index) == STATUS_AMBIGUOUS

    def counts(self):
        return {status: sum(1 for v in self.status.values() if v == status)
                for status in STATUSES}

    def phase_scores(self):
        phases = self.result.get("phases") or []
        return [sc for sc in phases if isinstance(sc, dict)]

    def phase(self, phase_no):
        """Запись фазы по 1-based номеру (или None)."""
        phases = self.phase_scores()
        phase_no = _as_int(phase_no)
        if phase_no is not None and 1 <= phase_no <= len(phases):
            return phases[phase_no - 1]
        return None


def peak_view(result=None, peaks=None):
    """Короткий конструктор: ``PeakView(result, peaks)``."""
    return PeakView(result, peaks)


def shift_info(result):
    """(режим, величина, прирост) сдвига нуля — с дефолтами для старой структуры."""
    result = result or {}
    mode = result.get("shift_mode")
    shift = result.get("shift") or 0.0
    try:
        shift = float(shift)
    except (TypeError, ValueError):
        shift = 0.0
    if not mode:
        mode = "учтён" if abs(shift) >= 0.03 else "не применялся"
    return str(mode), shift, result.get("shift_gain")


if __name__ == "__main__":   # небольшая самопроверка структуры
    demo_peaks = [{"two_theta": 10.0, "intensity": 100.0},
                  {"two_theta": 20.0, "intensity": 1.0},
                  {"two_theta": 30.0, "intensity": 40.0}]
    demo = {"phases": [{"phase": {"formula": "SiO2", "name": "Quartz"}, "fom": 80}],
            "per_peak": {0: [{"phase_no": 1, "dtheta": 0.02}]},
            "unassigned": [2], "peak_status": {0: "assigned", 2: "unexplained"},
            "ambiguous": [], "shift_mode": "auto", "settings": {"min_obs_int": 3.0}}
    view = PeakView(demo, demo_peaks)
    print("space:", view.space, "legacy:", view.legacy, "counts:", view.counts())
    print("codes:", phase_codes([{"name": "Quartz", "formula": "SiO2"},
                                 {"name": "Calcite", "formula": "CaCO3"}]))
    print("label:", phase_label(demo["phases"][0]))
