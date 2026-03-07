"""Validate static frontend assets (syntax checks)."""
import subprocess
import shutil

import pytest

from pathlib import Path

JS_DIR = Path(__file__).resolve().parent.parent / "app" / "static" / "js"
JS_FILES = sorted(JS_DIR.glob("*.js"))


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
@pytest.mark.parametrize("js_file", JS_FILES, ids=lambda p: p.name)
def test_js_syntax(js_file):
    """Each JS file must pass Node.js syntax check (node --check)."""
    result = subprocess.run(
        ["node", "--check", str(js_file)],
        capture_output=True, text=True, timeout=10,
    )
    assert result.returncode == 0, (
        f"{js_file.name} has syntax errors:\n{result.stderr}"
    )
