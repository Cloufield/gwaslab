"""Simple test runner for the repository.

Options:
- --start-dir: directory containing tests (default: test/)
- --pattern: discovery glob for test files (default: test_*.py)
- --verbosity: output verbosity for the runner (int)
- --failfast: stop on first failure/error
- --buffer: capture test output (stdout/stderr) during run
- --k: substring filter for test ids, similar to pytest -k
- --shuffle: run tests in random order; use with --seed for reproducibility
- --seed: integer seed to make --shuffle deterministic
"""

import os
import sys
import unittest
import argparse
import random


def _iter_tests(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            for t in _iter_tests(item):
                yield t
        else:
            yield item


def _build_suite(start_dir, pattern, k=None, shuffle=False, seed=None):
    loader = unittest.defaultTestLoader
    discovered = loader.discover(start_dir=start_dir, pattern=pattern)
    if k is None and not shuffle:
        return discovered
    cases = list(_iter_tests(discovered))
    if k is not None:
        cases = [t for t in cases if k in t.id()]
    if shuffle:
        if seed is not None:
            random.seed(seed)
        random.shuffle(cases)
    return unittest.TestSuite(cases)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the gwaslab test suite",
        epilog="Example: python test/run_all_tests.py --k Sumstats --shuffle --seed 42",
    )
    parser.add_argument(
        "--start-dir",
        default=os.path.dirname(__file__),
        help="Directory to start discovery from (defaults to the test folder)",
    )
    parser.add_argument(
        "--pattern",
        default="test_*.py",
        help="Glob pattern for discovering tests (default: test_*.py)",
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=2,
        help="Verbosity level for the test runner (0-3)",
    )
    parser.add_argument(
        "--failfast",
        action="store_true",
        help="Stop on first failure or error",
    )
    parser.add_argument(
        "--buffer",
        action="store_true",
        help="Capture stdout/stderr during tests",
    )
    parser.add_argument(
        "--k",
        help="Substring filter applied to test ids (similar to pytest -k)",
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Shuffle test execution order",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Seed for shuffling to make order reproducible",
    )
    args = parser.parse_args(argv)

    root = os.path.abspath(os.path.join(args.start_dir, ".."))
    src = os.path.join(root, "src")
    if src not in sys.path:
        sys.path.insert(0, src)

    suite = _build_suite(args.start_dir, args.pattern, k=args.k, shuffle=args.shuffle, seed=args.seed)
    result = unittest.TextTestRunner(verbosity=args.verbosity, failfast=args.failfast, buffer=args.buffer).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
