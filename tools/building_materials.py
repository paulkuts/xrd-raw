#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Профиль типичных фаз строительных материалов (мягкий приор для phases.explain).

Что это НЕ делает:
  * это не фильтр — ни одна фаза из базы не исключается;
  * это не источник эталонов — все карточки берутся из phases_db.json.gz как есть;
  * приор — только мягкий бонус к ранжированию (+3 балла к рангу максимум) и
    пометка ``prior`` в отчёте, чтобы человек видел, что фаза «ожидаемая».

Откуда взяты идентификаторы
---------------------------
Список собран ПРОВЕРКОЙ по фактической базе (32491 фаза): для каждой карточки
сверены формула, название и d сильнейшей линии с литературным значением
(допуск ~1 %). Карточки, у которых ячейка/паттерн не соответствуют минералу,
в профиль НЕ включены — см. ``NOT_INCLUDED`` и комментарии ниже.

Что в базе отсутствует (не выдумываем)
--------------------------------------
* C4AF / браунмиллерит Ca2AlFeO5, C4AF-подобные — нет записи (единственное
  похожее — Ca2Al2FeSi3O13, это другой минерал).
* Гидрогранат C3AH6, C4AH13, моносульфоалюминат (AFm) — записей нет
  (поиск по name/formula: monosulf*, hydrogarnet, Ca4Al2... — пусто).
* Монтмориллонит, тальк, рутил, халцедон — записей с корректным паттерном нет.
* Галит NaCl: COD 9000629 даёт d(200)=2.727 Å вместо 2.821 Å (ошибка ячейки),
  поэтому в профиль не включён, хотя галит — типичная фаза высолов.
* Оксид кальция CaO: COD 9006742 «Lime» даёт d(200)=2.494 Å вместо 2.405 Å —
  не включён.
* Портландит: в базе две карточки; COD 9016940 (CaO2) битая (d сильнейшей
  линии 2.895 Å вместо 2.628 Å) — включена только COD 9000112.

Дубликаты одного минерала с разным качеством
--------------------------------------------
Приор привязан к конкретному id карточки, а не к формуле. Это осознанно:
фаза «кварц» существует в базе 165 раз, и большинство карточек — чужие
структуры с той же формулой. Совпадение по формуле давало бы бонус мусору.

Публичный API
-------------
    PROFILE                 — кортеж записей профиля
    BY_ID                   — {id карточки: запись}
    prior_weight(phase)     — 0.0 .. 1.0 (вес приора для ранжирования)
    prior_tags(phase)       — [] или ["building_materials"]
    label(phase)            — русское название или None
    group(phase)            — группа («Карбонаты», ...) или None
    verify(db)              — (найденные, отсутствующие) id по реальной базе
    summary()               — строка состава профиля для отчёта
"""

PROFILE_NAME = "building_materials"

# weight: 1.0 — эталон соответствует минералу по d и по относительным
# интенсивностям; 0.5 — прокси/полиморф/карточка с оговоркой (бонус вдвое меньше).
PROFILE = (
    # ---------------------------------------------------------- каркасные силикаты
    dict(key="quartz", id="COD 1011097", label="Кварц", formula="SiO2",
         group="Каркасные силикаты", weight=1.0,
         note="Quartz low, 26.66° (d 3.343 Å) — совпадает с литературой"),
    # ------------------------------------------------------------------ карбонаты
    dict(key="calcite", id="COD 4502443", label="Кальцит", formula="CaCO3",
         group="Карбонаты", weight=1.0,
         note="29.53° (d 3.025 Å); в базе подписана формулой, не словом Calcite"),
    dict(key="aragonite", id="COD 9015893", label="Арагонит", formula="CaCO3",
         group="Карбонаты", weight=1.0, note="26.24° (d 3.396 Å)"),
    dict(key="vaterite", id="COD 9013565", label="Ватерит", formula="CaCO3",
         group="Карбонаты", weight=0.5,
         note="полиморф CaCO3, редко; 26.64° (d 3.346 Å)"),
    dict(key="monohydrocalcite", id="COD 9010479", label="Моногидрокальцит",
         formula="CaH2CO4", group="Карбонаты", weight=0.5,
         note="20.52° (d 4.328 Å); продукт карбонизации"),
    dict(key="dolomite", id="COD 9000573", label="Доломит", formula="CaMg(CO3)2",
         group="Карбонаты", weight=1.0, note="31.03° (d 2.882 Å)"),
    dict(key="magnesite", id="COD 9007692", label="Магнезит", formula="MgCO3",
         group="Карбонаты", weight=1.0,
         note="32.69° (d 2.740 Å); карточки 9002810/9002813 дают 2.59/2.50 Å — не включены"),
    # ---------------------------------------------------------------- гидроксиды
    dict(key="portlandite", id="COD 9000112", label="Портландит", formula="Ca(HO)2",
         group="Гидроксиды", weight=1.0,
         note="18.21°/34.25° (d 4.90/2.618 Å); карточка 9016940 битая — не включена"),
    dict(key="brucite", id="COD 9003875", label="Брусит", formula="Mg(HO)2",
         group="Гидроксиды", weight=1.0,
         note="18.32°/37.39° (d 4.77/2.405 Å); 9003876/9003877 сильно смещены"),
    # ------------------------------------------------------------------ сульфаты
    dict(key="gypsum", id="COD 5000039", label="Гипс", formula="CaSO6",
         group="Сульфаты", weight=1.0,
         note="11.64°/20.74°/29.14° (d 7.59/4.28/3.06 Å); H в формуле базы опущен"),
    dict(key="bassanite", id="COD 9012208", label="Бассанит (полуводный гипс)",
         formula="Ca2S2O9", group="Сульфаты", weight=0.5,
         note="14.76°/29.76° (d 6.00/3.00 Å)"),
    dict(key="anhydrite", id="COD 9017826", label="Ангидрит", formula="CaSO4",
         group="Сульфаты", weight=1.0,
         note="25.52° (d 3.490 Å); карточка 9012207 (14.67°) не соответствует ангидриту"),
    dict(key="ettringite", id="COD 9011576", label="Эттрингит", formula="Ca6Al2S3O49",
         group="Сульфаты", weight=1.0,
         note="9.09° (d 9.726 Å); H2O в формуле базы сведена в O"),
    dict(key="thaumasite", id="COD 9010878", label="Таумазит", formula="Ca3SiCSO25",
         group="Сульфаты", weight=1.0, note="9.23° (d 9.576 Å)"),
    dict(key="thenardite", id="COD 1011184", label="Тенардит", formula="Na2SO4",
         group="Сульфаты", weight=0.5,
         note="32.34° (d 2.768 Å); типичная фаза сульфатных высолов"),
    # ----------------------------------------------------- клинкер и цементные фазы
    dict(key="alite", id="COD 1540704", label="Алит (C3S)", formula="Ca3SiO5",
         group="Клинкер", weight=1.0,
         note="32.22°/29.38° (d 2.778 Å); моноклинная форма"),
    dict(key="belite", id="COD 9012681", label="Белит (γ-C2S)", formula="Ca2SiO4",
         group="Клинкер", weight=0.5,
         note="calcio-olivine, 32.85°/29.69°; β-C2S в базе нет — прокси, полиморф γ"),
    dict(key="c3a", id="COD 9014359", label="Алюминат кальция (C3A)",
         formula="Ca3Al2O6", group="Клинкер", weight=1.0,
         note="33.24°/47.72° (d 2.696 Å)"),
    dict(key="mayenite", id="COD 1011034", label="Майенит (C12A7)",
         formula="Ca6Al7O20", group="Клинкер", weight=0.5,
         note="18.18°/33.54°; редкая фаза, обычно в высокоглинозёмистых цементах"),
    dict(key="ternesite", id="COD 9009846", label="Тернесит", formula="Ca5Si2SO12",
         group="Клинкер", weight=0.5,
         note="31.36°/27.91°; сульфосиликат клинкера (сульфоалюминатные цементы)"),
    # ------------------------------------------------------------- глинистые минералы
    dict(key="kaolinite", id="COD 1011045", label="Каолинит", formula="Al2Si2O9",
         group="Глинистые", weight=1.0,
         note="12.38° (d 7.153 Å); вторая карточка 1566358 (K-II, 7.089 Å) — 0.5"),
    dict(key="kaolinite_k2", id="COD 1566358", label="Каолинит K-II",
         formula="Al2Si2H4O9", group="Глинистые", weight=0.5,
         note="12.49° (d 7.089 Å)"),
    dict(key="illite", id="COD 9013719", label="Иллит (прокси-структура)",
         formula="KAl2(SiO3)4", group="Глинистые", weight=0.5,
         note="8.82°/24.38°; иллит — группа, карточка является её структурным прокси"),
    # -------------------------------------------------------------- серпентиновые
    dict(key="lizardite", id="COD 9004511", label="Лизардит", formula="Mg3Si2H4O9",
         group="Серпентиновые", weight=1.0,
         note="12.16°/35.97° (d 7.276 Å) — корректный серпентиновый паттерн"),
    dict(key="antigorite", id="COD 9003103", label="Антигорит",
         formula="Mg48Si34H62O147", group="Серпентиновые", weight=1.0,
         note="12.19° (d 7.261 Å)"),
    dict(key="chrysotile", id="COD 1010960", label="Хризотил", formula="Mg3Si2O9",
         group="Серпентиновые", weight=0.5,
         note="паттерн карточки anomalous: максимум 16.66° (d 5.32 Å) вместо "
              "базального 7.3 Å (12.1°, I=8); поэтому бонус половинный"),
    # ---------------------------------------------------------------- C-S-H (прокси)
    dict(key="tobermorite", id="COD 9005497", label="Тоберморит (прокси C-S-H)",
         formula="Ca2Si3O11", group="C-S-H", weight=0.5,
         note="7.86° (d 11.24 Å); C-S-H в базе нет, тоберморит — обычный прокси"),
    dict(key="xonotlite", id="COD 9008439", label="Ксонотлит", formula="Ca6Si6O19",
         group="C-S-H", weight=0.5, note="28.90°/27.45°"),
    dict(key="afwillite", id="COD 9007428", label="Афвиллит", formula="Ca3(SiO5)2",
         group="C-S-H", weight=0.5, note="28.08°/31.58°"),
    # ------------------------------------------------------------- железо, прочее
    dict(key="hematite", id="COD 9009782", label="Гематит", formula="Fe2O3",
         group="Железистые", weight=1.0, note="33.23° (d 2.696 Å)"),
    dict(key="goethite", id="COD 9003078", label="Гётит", formula="FeO2",
         group="Железистые", weight=1.0, note="21.40° (d 4.152 Å)"),
    dict(key="periclase", id="COD 9006459", label="Периклаз", formula="MgO",
         group="Железистые", weight=0.5, note="42.56° (d 2.124 Å)"),
    dict(key="anatase", id="COD 9008216", label="Анатаз", formula="TiO2",
         group="Железистые", weight=0.5,
         note="25.18° (d 3.537 Å); карточка 1010942 даёт 3.466 Å — не включена"),
)

# Записи, которые пришлось исключить, и причина — для отчёта и для проверки.
NOT_INCLUDED = (
    ("C4AF / браунмиллерит Ca2AlFeO5", "записи в базе нет"),
    ("C3AH6 (гидрогранат), C4AH13", "записей в базе нет"),
    ("Моносульфоалюминат (AFm)", "записей в базе нет (поиск monosulf*/hydrogarnet пуст)"),
    ("Галит NaCl (COD 9000629)", "d(200)=2.727 Å вместо 2.821 Å — ошибка ячейки"),
    ("Оксид кальция CaO (COD 9006742)", "d(200)=2.494 Å вместо 2.405 Å — ошибка ячейки"),
    ("Портландит COD 9016940 (CaO2)", "битая карточка: d=2.895 вместо 2.628 Å"),
    ("Монтмориллонит, тальк, рутил", "записей с корректным паттерном нет"),
)

BY_ID = {e["id"]: e for e in PROFILE}

GROUPS = []
for _e in PROFILE:
    if _e["group"] not in GROUPS:
        GROUPS.append(_e["group"])


def _entry(phase):
    if not phase:
        return None
    return BY_ID.get(phase.get("id"))


def prior_weight(phase):
    """Вес приора 0.0..1.0. 0.0 — фазы нет в профиле."""
    e = _entry(phase)
    return float(e["weight"]) if e else 0.0


def prior_tags(phase):
    """Пометка для отчёта: [] или ["building_materials"]."""
    return [PROFILE_NAME] if _entry(phase) else []


def label(phase):
    """Русское название фазы из профиля (или None)."""
    e = _entry(phase)
    return e["label"] if e else None


def group(phase):
    """Группа профиля (или None)."""
    e = _entry(phase)
    return e["group"] if e else None


def note(phase):
    """Пояснение к карточке профиля (или None)."""
    e = _entry(phase)
    return e["note"] if e else None


def verify(db):
    """Сверить профиль с реальной базой.

    Возвращает (found, missing) — два списка id. Пустой missing означает, что
    все карточки профиля действительно есть в базе.
    """
    ids = set()
    for ph in db.get("phases", ()):
        ids.add(ph.get("id"))
    found = [e["id"] for e in PROFILE if e["id"] in ids]
    missing = [e["id"] for e in PROFILE if e["id"] not in ids]
    return found, missing


def summary():
    """Короткая строка состава профиля — для отчёта."""
    counts = {}
    for e in PROFILE:
        counts[e["group"]] = counts.get(e["group"], 0) + 1
    parts = ["{}: {}".format(g, counts[g]) for g in GROUPS if g in counts]
    return "профиль {} — {} карточек ({})".format(PROFILE_NAME, len(PROFILE), "; ".join(parts))


if __name__ == "__main__":
    import phases as _phases

    print(summary())
    for e in PROFILE:
        print("  {:<11s} {:<14s} {:<22s} {:<16s} w={}".format(
            e["id"], e["formula"][:14], e["label"][:22], e["group"], e["weight"]))
    _db = _phases.load_db()
    _found, _missing = verify(_db)
    print("найдено в базе: {}/{}, отсутствует: {}".format(
        len(_found), len(PROFILE), _missing or "нет"))
