# 项目概览
- 目的：一个极简的 Python 项目模板，主入口 `main.py` 打印 "Hello from noname!"，适合用于拓展功能或快速验证 Python 环境。
- 技术栈：Python 3.14（`pyproject.toml` 指定 `requires-python >=3.14`），采用 `pyproject` 原生包配置，依赖为空。项目通过 `.venv` 管理虚拟环境。
- 结构：根目录包含 `pyproject.toml`、`main.py`、`.gitignore`、`.python-version` 和空的 `README.md`；源代码即 `main.py`，目前没有子目录或多个模块。
- 运行：可在项目根激活 `.venv`（如 `.\.venv\Scripts\Activate.ps1`）后执行 `python main.py`，或直接 `python main.py`（确保 `python >=3.14`）。