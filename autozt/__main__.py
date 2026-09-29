# -*- coding: utf-8 -*-
"""python -m autozt 入口。"""
import os
import sys

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass

# `python -m autozt agent ...` 与 bin/autozt、pip 入口保持同一路由实现。
from autozt.cli import route_subcommand as _route

_rc = _route(sys.argv[1:])
if _rc is not None:
    raise SystemExit(_rc)

from autozt import main

if __name__ == "__main__":
    try:
        main()
        sys.stdout.flush()
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)
