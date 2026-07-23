"""Shared per-test timing log for unittest and pytest runners."""

from __future__ import annotations

import re
import time
import unittest
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Literal, Optional

LOG_PATH = Path(__file__).resolve().parent / "test_timings.log"
_PART_GLOB = "test_timings.log.*.part"
_TIMING_LINE_RE = re.compile(r"^([\d.]+)s\s+(\S+)\s+(.+)$")

_PYTEST_OUTCOME_TO_STATUS = {
    "passed": "PASS",
    "failed": "FAIL",
    "skipped": "SKIP",
    "error": "ERROR",
    "xfailed": "XFAIL",
    "xpassed": "XPASS",
}

RunMode = Literal["serial", "xdist_worker", "xdist_controller"]


def part_path(worker_id: str) -> Path:
    return LOG_PATH.parent / f"{LOG_PATH.name}.{worker_id}.part"


def parse_timing_line(line: str) -> tuple[float, str, str] | None:
    match = _TIMING_LINE_RE.match(line.strip())
    if match is None:
        return None
    return float(match.group(1)), match.group(2), match.group(3)


def format_timing_line(duration_seconds: float, status: str, test_id: str) -> str:
    return f"{duration_seconds:.3f}s  {status:<6}  {test_id}"


def is_xdist_worker(config) -> bool:
    return hasattr(config, "workerinput")


def is_xdist_controller(config) -> bool:
    if is_xdist_worker(config):
        return False
    if config.pluginmanager.hasplugin("xdist") is None:
        return False
    return config.getoption("dist", default="no") != "no"


def _cleanup_part_files() -> None:
    for part in LOG_PATH.parent.glob(_PART_GLOB):
        part.unlink(missing_ok=True)


class TimingLog:
    """Timing log for serial pytest, xdist workers, or xdist controller merge."""

    def __init__(self) -> None:
        self._fp: Optional[IO[str]] = None
        self._mode: RunMode = "serial"
        self._header_line: Optional[str] = None
        self.count = 0
        self.total_seconds = 0.0

    def start_run(
        self,
        runner_name: str,
        *,
        mode: RunMode = "serial",
        worker_id: Optional[str] = None,
        worker_count: Optional[int] = None,
    ) -> None:
        self._mode = mode
        self.count = 0
        self.total_seconds = 0.0
        started = datetime.now(timezone.utc).isoformat()

        if mode == "xdist_controller":
            _cleanup_part_files()
            workers = worker_count if worker_count is not None else "?"
            self._header_line = f"# run started {started} runner={runner_name} workers={workers}"
            self._fp = None
            return

        if mode == "xdist_worker":
            if worker_id is None:
                raise ValueError("worker_id is required for xdist_worker mode")
            self._fp = part_path(worker_id).open("w", encoding="utf-8")
            return

        self._header_line = f"# run started {started} runner={runner_name}"
        self._fp = LOG_PATH.open("w", encoding="utf-8")
        self._fp.write(f"{self._header_line}\n")
        self._fp.flush()

    def record(self, duration_seconds: float, status: str, test_id: str) -> None:
        if self._mode == "xdist_controller":
            return
        if self._fp is None:
            return
        self._fp.write(format_timing_line(duration_seconds, status, test_id) + "\n")
        self._fp.flush()
        self.count += 1
        self.total_seconds += duration_seconds

    def finish_run(self) -> None:
        if self._mode == "xdist_worker":
            if self._fp is not None:
                self._fp.close()
                self._fp = None
            return

        if self._mode == "xdist_controller":
            self._merge_xdist_parts()
            return

        if self._fp is None:
            return
        self._fp.write(f"# total: {self.count} tests, {self.total_seconds:.3f}s\n")
        self._fp.close()
        self._fp = None

    def _merge_xdist_parts(self) -> None:
        entries: list[tuple[float, str, str]] = []
        part_files = sorted(LOG_PATH.parent.glob(_PART_GLOB))

        for part in part_files:
            for raw_line in part.read_text(encoding="utf-8").splitlines():
                parsed = parse_timing_line(raw_line)
                if parsed is not None:
                    entries.append(parsed)

        entries.sort(key=lambda item: item[2])
        count = len(entries)
        total_seconds = sum(entry[0] for entry in entries)
        header = self._header_line or "# run started unknown runner=pytest"

        with LOG_PATH.open("w", encoding="utf-8") as fp:
            fp.write(f"{header}\n")
            for duration, status, test_id in entries:
                fp.write(format_timing_line(duration, status, test_id) + "\n")
            fp.write(f"# total: {count} tests, {total_seconds:.3f}s\n")

        _cleanup_part_files()
        self.count = count
        self.total_seconds = total_seconds


_timing_log = TimingLog()


def get_timing_log() -> TimingLog:
    return _timing_log


def outcome_to_status(outcome: str) -> str:
    return _PYTEST_OUTCOME_TO_STATUS.get(outcome, outcome.upper())


class TimedTextTestResult(unittest.TextTestResult):
    """TextTestResult that records per-test duration to the shared timing log."""

    def startTest(self, test: unittest.TestCase) -> None:  # noqa: N802
        self._test_started_at = time.perf_counter()
        super().startTest(test)

    def _elapsed(self) -> float:
        return time.perf_counter() - self._test_started_at

    def _log(self, test: unittest.TestCase, status: str) -> None:
        get_timing_log().record(self._elapsed(), status, test.id())

    def addSuccess(self, test: unittest.TestCase) -> None:  # noqa: N802
        super().addSuccess(test)
        self._log(test, "PASS")

    def addFailure(self, test: unittest.TestCase, err) -> None:  # noqa: N802
        super().addFailure(test, err)
        self._log(test, "FAIL")

    def addError(self, test: unittest.TestCase, err) -> None:  # noqa: N802
        super().addError(test, err)
        self._log(test, "ERROR")

    def addSkip(self, test: unittest.TestCase, reason: str) -> None:  # noqa: N802
        super().addSkip(test, reason)
        self._log(test, "SKIP")

    def addExpectedFailure(self, test: unittest.TestCase, err) -> None:  # noqa: N802
        super().addExpectedFailure(test, err)
        self._log(test, "XFAIL")

    def addUnexpectedSuccess(self, test: unittest.TestCase) -> None:  # noqa: N802
        super().addUnexpectedSuccess(test)
        self._log(test, "XPASS")


class TimedTextTestRunner(unittest.TextTestRunner):
    """TextTestRunner that uses TimedTextTestResult."""

    def __init__(self, **kwargs) -> None:
        kwargs.setdefault("resultclass", TimedTextTestResult)
        super().__init__(**kwargs)
