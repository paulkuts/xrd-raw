# -*- coding: utf-8 -*-
"""Поиск фаз в дампе по набору элементов и диапазону сильнейшей линии."""
import json
import sys


def top(r):
    return max(r["peaks"], key=lambda x: x[1])


def main():
    need = set(sys.argv[1].split(",")) if len(sys.argv) > 1 else set()
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
    rs = [json.loads(l) for l in open("cod_phases.jsonl", encoding="utf-8")]
    hits = []
    for r in rs:
        if need and not need.issubset(set(r["elements"])):
            continue
        d1 = top(r)[0]
        if lo <= d1 <= hi:
            hits.append((d1, r))
    hits.sort(key=lambda x: -x[0])
    print("совпадений:", len(hits))
    for d1, r in hits[:25]:
        print("  {:<24} {:<14} {:<30} d1={:6.3f}".format(
            r["formula"][:24], r["id"], (r.get("name") or "")[:30], d1))


if __name__ == "__main__":
    main()
