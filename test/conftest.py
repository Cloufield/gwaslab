"""Ensure in-repo ``src/`` is preferred over any installed ``gwaslab`` package."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from timing_log import get_timing_log, is_xdist_controller, is_xdist_worker, outcome_to_status

ROOT = Path(__file__).resolve().parent.parent
SRC = str(ROOT / "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


def pytest_sessionstart(session: pytest.Session) -> None:
    config = session.config
    if is_xdist_worker(config):
        get_timing_log().start_run(
            "pytest",
            mode="xdist_worker",
            worker_id=config.workerinput["workerid"],
        )
    elif is_xdist_controller(config):
        get_timing_log().start_run(
            "pytest",
            mode="xdist_controller",
            worker_count=config.getoption("numprocesses"),
        )
    else:
        get_timing_log().start_run("pytest")


def pytest_runtest_logreport(report: pytest.TestReport) -> None:
    if report.when == "call":
        should_log = True
    elif report.when == "setup" and report.outcome != "passed":
        should_log = True
    else:
        should_log = False

    if not should_log:
        return

    duration = report.duration if report.duration is not None else 0.0
    status = outcome_to_status(report.outcome)
    get_timing_log().record(duration, status, report.nodeid)


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    get_timing_log().finish_run()
