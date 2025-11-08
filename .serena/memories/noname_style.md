# 编码与风格规范
- 缩进：使用 4 个空格，不混用制表符；关键块（函数、条件）前后留空行以增加可读性。
- 命名：模块/文件小写，函数和变量使用蛇形 (`snake_case`)，常量使用全大写。
- 类型与注释：函数可选写简单 docstring；鼓励保持 `main.py` 逻辑清晰、语句简洁；暂无特定 linter，建议遵循 PEP 8。
- 工具：未配置格式化器，可在 `.venv` 中使用 `python -m pip install black` 后运行 `black main.py` 进行格式化；如果增加依赖，更新 `pyproject.toml`。