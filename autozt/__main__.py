# -*- coding: utf-8 -*-
"""python -m autozt 入口。"""
import os
import sys

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass

# `python -m autozt agent ...` 与 bin/autozt 保持同一入口；只识别裸命令词，
# 避免材料名恰好叫 ``agent`` 或 ``mcp`` 时误路由。
_VALUE_FLAGS = {'-c', '--config', '-tt', '-p', '-j', '-job', '-status', '--status',
                '-x', '--exclude', '--host', '-u', '--user'}


def _bare_index(name):
    skip = False
    for index, token in enumerate(sys.argv[1:]):
        if skip:
            skip = False
            continue
        if token in _VALUE_FLAGS:
            skip = True
            continue
        if token == name:
            return index + 1
    return None


_i_mcp = _bare_index('mcp')
_i_agent = _bare_index('agent')
if _i_mcp is not None and (_i_agent is None or _i_mcp < _i_agent):
    from autozt import mcp as _mcp_intercept
    raise SystemExit(_mcp_intercept.main(sys.argv[_i_mcp + 1:]))
if _i_agent is not None:
    from autozt import agent_cli as _agent_intercept
    raise SystemExit(_agent_intercept.main(sys.argv[1:]))

from autozt import main

if __name__ == "__main__":
    try:
        main()
        sys.stdout.flush()
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)
