# -*- coding: utf-8 -*-
"""Согласование вывода Python с кодовой страницей консоли Windows.

Проблема: cmd после `chcp 65001` ждёт UTF-8, а встроенный Python пишет в
ANSI (cp1251) — в окне получаются кракозябры. Функция спрашивает у Windows
текущую кодовую страницу вывода и переключает stdout/stderr на неё.
На Linux/macOS ничего не меняет.
"""
import sys

_MAP = {65001: "utf-8", 866: "cp866", 1251: "cp1251", 855: "cp855", 437: "cp437"}


def console_encoding():
    """Кодовая страница консоли вывода ('utf-8', 'cp866', ...) или None."""
    try:
        import ctypes

        cp = int(ctypes.windll.kernel32.GetConsoleOutputCP())
    except Exception:  # noqa: BLE001  (не Windows или нет консоли)
        return None
    if cp <= 0:
        return None
    return _MAP.get(cp, "cp{}".format(cp))


def setup():
    """Вызывать в самом начале CLI-скрипта, до любого вывода."""
    enc = console_encoding()
    for stream in (sys.stdout, sys.stderr):
        try:
            if enc:
                stream.reconfigure(encoding=enc, errors="replace")
            else:
                stream.reconfigure(errors="replace")
        except Exception:  # noqa: BLE001
            pass
    return enc
