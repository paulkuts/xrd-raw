#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""XRD-конвертер: меню для файлов .raw дифрактометра Thermo ARL X'TRA (WinXRD).

Запускается через START.bat. Полностью портативно: используется встроенный
python из папки python\\, ничего доустанавливать не нужно.
"""

import os
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent      # корень портативной папки
TOOLS = Path(__file__).resolve().parent
INPUT = BASE / "input"
PYEXE = BASE / "python" / "python.exe"

sys.path.insert(0, str(TOOLS))
from winxrd2xy import convert_data, parse, write_csv, write_xy   # noqa: E402

LINE = "=" * 62
last_result = None


def head(text):
    # безопасно: строка константная, без подстановки пользовательских данных
    os.system("cls" if os.name == "nt" else "clear")
    print(LINE)
    print("  XRD-конвертер  |  .raw от Thermo ARL X'TRA (WinXRD 2.0)")
    print(LINE)
    if text:
        print(text)
        print()


def ask(prompt, default=None):
    try:
        ans = input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        return None
    return ans or default


def size_str(p):
    kb = Path(p).stat().st_size / 1024
    return "{:.0f} КБ".format(kb)


def scan_raw():
    """Ищем .raw в папке input, в корне инструмента и в текущей папке."""
    found = []
    for folder in (INPUT, BASE, Path.cwd()):
        for pattern in ("*.raw", "*.RAW"):
            for f in sorted(folder.glob(pattern)):
                try:
                    r = f.resolve()
                except OSError:
                    continue
                if r not in [x.resolve() for x in found] and r.is_file():
                    found.append(f)
    return found


def pick_file(prefill=None):
    if prefill:
        p = Path(prefill)
        if p.is_file():
            return p
        print("Файл не найден: {}".format(prefill))
    files = scan_raw()
    if not files:
        print("В папке input\\ нет файлов .raw.")
        print("Положите файл в input\\ или впишите полный путь ниже.")
        print()
        path = ask("Путь к файлу .raw (Enter - отмена): ")
        if not path:
            return None
        p = Path(path.strip('"').strip("'"))
        return p if p.is_file() else None

    print("Найденные файлы .raw:")
    for i, f in enumerate(files, 1):
        where = "input\\" if f.parent == INPUT else str(f.parent)
        print("  {}. {:<28} {}   [{}]".format(i, f.name, size_str(f), where))
    print()
    print("  0. вписать путь вручную")
    print("  Enter. отмена")
    print()
    ans = ask("Номер файла: ")
    if not ans:
        return None
    if ans == "0":
        path = ask("Путь к файлу .raw: ")
        if not path:
            return None
        p = Path(path.strip('"').strip("'"))
        return p if p.is_file() else None
    if ans.isdigit() and 1 <= int(ans) <= len(files):
        return files[int(ans) - 1]
    print("Не понял номер.")
    return None


def pick_theme():
    print()
    print("Тема графика:")
    print("  1. тёмная  - для экрана, презентаций (по умолчанию)")
    print("  2. светлая - для печати, диссертации, статей")
    ans = ask("Тема (Enter - тёмная): ")
    return "light" if ans == "2" else "dark"


def show_result(files):
    global last_result
    print()
    for f in files:
        print("  ГОТОВО: {}  ({})".format(f, size_str(f)))
    if files:
        last_result = Path(files[0]).parent


def print_summary(src):
    data, step, meta = parse(src)
    imax = max(i for _, i in data)
    top = max(data, key=lambda p: p[1])
    print("  точек: {}, 2θ {:.2f}-{:.2f} град, шаг {:.2f} град, анод {}".format(
        len(data), data[0][0], data[-1][0], step, meta.get("anode", "Cu")))
    print("  максимум: {} имп. при 2θ = {:.2f}°".format(imax, top[0]))
    return data


def do_table(src):
    """Таблица Excel (.xlsx) + CSV."""
    import xlsx_out
    data, step, meta = parse(src)
    csv = write_csv(data, src)
    xlsx = xlsx_out.make_xlsx(src)
    show_result([xlsx, csv])
    return [xlsx, csv]


def do_xy(src):
    """Данные для XRD-пакетов: .xy + .csv."""
    data, step, meta = parse(src)
    xy = write_xy(data, src)
    csv = write_csv(data, src)
    show_result([xy, csv])
    return [xy, csv]


def do_graph(src, theme):

    png, svg, meta = plot_xrd.make(src, theme=theme)
    show_result([png, svg])
    return [png, svg]


def do_identify(src, theme="dark", both_themes=False):
    """Идентификация фаз и элементов по пикам."""
    import identify
    ans, files, txt = identify.run(src, theme=theme, both_themes=both_themes,
                                   want_xlsx=True)
    print()
    print(txt)
    print()
    show_result([str(f) for f in files])
    return files


def do_everything(src, theme="dark", both_themes=False):
    """Всё сразу: .xlsx + .csv + .xy + PNG + SVG + идентификация фаз."""
    import identify

    import xlsx_out
    data, step, meta = parse(src)
    xy = write_xy(data, src)
    csv = write_csv(data, src)
    try:
        ans, idfiles, _txt = identify.run(src, theme=theme,
                                          both_themes=both_themes, want_xlsx=False)
        xlsx = xlsx_out.make_xlsx(src, phases_result=ans["result"], peaks=ans["peaks"])
        results = [xlsx, csv, xy] + [str(f) for f in idfiles]
    except Exception as e:                      # noqa: BLE001 — без фаз, но данные не теряем
        print("  (идентификация фаз не удалась: {})".format(e))
        xlsx = xlsx_out.make_xlsx(src)
        results = [xlsx, csv, xy]
    # графики уже построены идентификацией — с подписями фаз, второй («чистый») вариант не нужен
    show_result(results)
    return results


def run_action(kind, src, theme=None):
    print_summary(src)
    print()
    print("Работаю...")
    if kind == "identify":
        do_identify(src, theme or "dark")
    elif kind == "table":
        do_table(src)
    elif kind == "xy":
        do_xy(src)
    elif kind == "graph":
        do_graph(src, theme)
    elif kind == "batch":
        do_everything(src)
    else:
        do_everything(src, theme)


def open_folder(folder):
    if not folder or not Path(folder).exists():
        print("Пока нечего открывать.")
        return
    if os.name == "nt":
        subprocess.Popen(['explorer', str(Path(folder))])
    else:
        print("Папка: {}".format(folder))


def main():
    # Windows-консоль бывает в cp866/cp1251/utf-8: подстраиваемся под её кодовую страницу
    try:
        from console import setup as _setup_console

        _setup_console()
    except Exception:  # noqa: BLE001
        for stream in (sys.stdout, sys.stderr):
            try:
                stream.reconfigure(errors="replace")
            except Exception:  # noqa: BLE001
                pass
    prefill = sys.argv[1] if len(sys.argv) > 1 else None
    if prefill and not Path(prefill).is_file():
        prefill = None
    while True:
        head("")
        print("  Что сделать с файлом .raw?   (от полного к простому)")
        print()
        print("  1. ИДЕНТИФИКАЦИЯ ФАЗ И ЭЛЕМЕНТОВ: пики, фазы, график, .txt + .xlsx")
        print("  2. Пакетно: все .raw из папки input  ->  всё сразу + идентификация")
        print("  3. Всё сразу:  .xlsx + .csv + .xy + PNG + SVG (выбор темы графика)")
        print("  4. Таблица Excel:  .xlsx (лист с данными + родная диаграмма Excel) + .csv")
        print("  5. Данные для XRD-пакетов:  .xy (2θ, интенсивность) + .csv")
        print("     (HighScore, Match!, DIFFRAC.EVA, FullProf, GSAS-II, Origin)")
        print("  6. Только график:  PNG + SVG (тёмная или светлая тема)")
        print("  7. Открыть папку с последним результатом")
        print("  0. Выход")
        print()
        ans = ask("Номер действия (Enter - выход): ")
        if ans is None or ans == "0":
            break

        if ans == "7":
            head("")
            open_folder(last_result)
            ask("Enter - назад в меню: ")
            continue

        try:
            if ans == "2":
                files, seen = [], set()
                # на Windows glob не различает регистр: *.raw и *.RAW дают дубли — гасим
                for f in list(INPUT.glob("*.raw")) + list(INPUT.glob("*.RAW")):
                    key = f.resolve()
                    if key not in seen:
                        seen.add(key)
                        files.append(f)
                files.sort()
                if not files:
                    head("В папке input\\ нет файлов .raw.")
                    ask("Enter - назад в меню: ")
                    continue
                head("Обрабатываю {} файл(ов) из input\\ (обе темы графиков)...".format(len(files)))
                for f in files:
                    print("\n--> {}".format(f.name))
                    run_action("batch", f)
                ask("\nEnter - назад в меню: ")
                continue

            if ans not in ("1", "3", "4", "5", "6"):
                continue

            head("")
            src = pick_file(prefill)
            if not src:
                continue

            print()
            print("Файл: {}".format(src))
            theme = pick_theme() if ans in ("1", "3", "6") else None
            print()
            run_action({"1": "identify", "3": "all", "4": "table",
                        "5": "xy", "6": "graph"}[ans], src, theme)

            print()
            if ask("Открыть папку с результатом? (1 - да, Enter - нет): ") == "1":
                open_folder(last_result)
            ask("Enter - назад в меню: ")

        except ValueError as e:
            head("Ошибка: {}".format(e))
            ask("Enter - назад в меню: ")
        except Exception as e:                      # noqa: BLE001 — показываем юзеру причину
            head("Не получилось: {}: {}".format(type(e).__name__, e))
            ask("Enter - назад в меню: ")

    print("\nПока.")


if __name__ == "__main__":
    main()
