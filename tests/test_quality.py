from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_ruff_check_passes() -> None:
    assert_command_succeeds([sys.executable, "-m", "ruff", "check", "."], cwd=REPO_ROOT)


def test_mypy_check_passes() -> None:
    assert_command_succeeds([sys.executable, "-m", "mypy"], cwd=REPO_ROOT)


def assert_command_succeeds(command: list[str], *, cwd: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, (
        f"{' '.join(command)} failed with exit code {completed.returncode}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
