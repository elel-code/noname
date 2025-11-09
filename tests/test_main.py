import io
import sys
from pathlib import Path

# 确保可从仓库根目录导入 main.py
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import main


def test_main_hello_mode_prints_expected():
    buf = io.StringIO()
    old = sys.stdout
    try:
        sys.stdout = buf
        rc = main.main(["--mode", "hello"])
    finally:
        sys.stdout = old
    out = buf.getvalue().strip()
    assert rc == 0
    assert out == "Hello from noname!"
