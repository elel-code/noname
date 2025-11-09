from __future__ import annotations

"""CLI 入口模块

- hello 模式：最小可运行输出（便于只读/无依赖环境自检）
- demo 模式：委派到 app.run_demo 执行单/多目标流程
"""

import argparse
import sys
from typing import TYPE_CHECKING

from app import run_demo

if TYPE_CHECKING:  # 仅类型检查与提示，避免运行时硬依赖
    from algorithms import Ingredient


def main(argv: list[str] | None = None) -> int:
    """CLI 入口：
    - 默认 hello 模式：打印最小可运行输出，满足快速自检。
    - demo 模式：运行算法演示（需要依赖与可写目录）。
    """

    parser = argparse.ArgumentParser(description="noname CLI")
    parser.add_argument("--mode", choices=("hello", "demo"), default="hello", help="运行模式")
    parser.add_argument("--ingredients", type=str, default=None, help="成分数据文件（CSV/JSON）")
    parser.add_argument(
        "--format", dest="ingredients_format", choices=("auto", "csv", "json"), default="auto", help="成分文件格式"
    )
    # 输出/可视化/算法控件
    parser.add_argument("--output-dir", type=str, default=None, help="可视化输出目录（默认取配置 viz.output_dir）")
    parser.add_argument("--no-fig", action="store_true", help="禁用可视化输出")
    parser.add_argument("--save-metrics", type=str, default=None, help="将 NSGA 候选指标保存为 CSV")
    parser.add_argument("--precision", type=int, default=3, help="打印指标的小数位数（默认 3）")
    parser.add_argument("--dpi", type=int, default=None, help="可视化输出 DPI（默认取配置 viz.dpi 或 200）")
    parser.add_argument("--colormap", type=str, default=None, help="NSGA 前沿图的 colormap（默认取配置 viz.colormap 或 viridis）")
    parser.add_argument("--font", type=str, default=None, help="首选中文字体名（可选）")
    parser.add_argument("--font-dir", type=str, default=None, help="字体文件目录（可选，作为回退）")
    parser.add_argument("--max-candidates", type=int, default=None, help="NSGA 候选解数量上限（优先于配置 nsga2.max_candidates）")
    # 可视化可调
    parser.add_argument("--figs", type=str, default=None, help="图表：pareto,bar,pie,hist,density,radar（逗号分隔）")
    parser.add_argument("--img-formats", type=str, default=None, help="输出格式：png,svg,pdf（逗号分隔）")
    parser.add_argument("--no-annotate", action="store_true", help="禁用帕累托图高亮点标注")
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("--responsive", action="store_true", help="启用图表自适应尺寸（默认取配置）")
    grp.add_argument("--no-responsive", action="store_true", help="禁用图表自适应尺寸")
    parser.add_argument("--title-prefix", type=str, default=None, help="为图表标题添加前缀文本")
    # 文本输出增强
    parser.add_argument("--summary-top-k", type=int, default=3, help="打印 Top-K 汇总（ACI 与最低毒性），默认 3，0 表示关闭")
    parser.add_argument("--export-json", type=str, default=None, help="导出汇总 JSON 路径（可选）")

    args = parser.parse_args(argv)

    if args.mode == "hello":
        print("Hello from noname!")
        return 0

    # demo 模式
    provided: list[Ingredient] | None = None
    if args.ingredients:
        try:
            from algorithms import load_ingredients

            fmt = None if args.ingredients_format == "auto" else args.ingredients_format
            provided = load_ingredients(args.ingredients, fmt=fmt)
        except Exception as e:  # noqa: BLE001 - 面向用户的直观提示
            print(f"加载成分文件失败：{e}", file=sys.stderr)
            return 3
    return run_demo(provided, args=args)


if __name__ == "__main__":  # pragma: no cover - 直接作为脚本运行
    raise SystemExit(main())
