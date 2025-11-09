from __future__ import annotations

import math
import warnings
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from algorithms import CandidateSolution, Ingredient

HighlightPoint = tuple[str, float, float]
DEFAULT_OUTPUT_DIR = Path("artifacts")
PREFERRED_FONTS: tuple[str, ...] = (
    "Microsoft YaHei",
    "SimHei",
    "Source Han Sans CN",
    "Noto Sans CJK SC",
    "PingFang SC",
)
FALLBACK_FONT_FILES: tuple[str, ...] = (
    "SourceHanSansCN-Regular.otf",
    "SourceHanSansSC-Regular.otf",
    "NotoSansSC-Regular.otf",
)

__all__ = [
    "configure_matplotlib_fonts",
    "render_nsga_visualizations",
    "ensure_output_dir",
]


def ensure_output_dir(output_dir: Path | None = None) -> Path:
    """确保图表输出目录存在。"""

    target = output_dir or DEFAULT_OUTPUT_DIR
    target.mkdir(parents=True, exist_ok=True)
    return target


def configure_matplotlib_fonts(
    preferred_fonts: Sequence[str] | None = None,
    fallback_dir: Path | None = None,
) -> str | None:
    """配置 Matplotlib 中文字体，返回最终生效的字体名称。"""

    preferred_fonts = tuple(preferred_fonts or PREFERRED_FONTS)
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    chosen_font = next((font for font in preferred_fonts if font in available_fonts), None)
    if chosen_font is None:
        font_dir = fallback_dir or (Path(__file__).resolve().parent / "resources" / "fonts")
        for filename in FALLBACK_FONT_FILES:
            font_path = font_dir / filename
            if not font_path.exists():
                continue
            font_manager.fontManager.addfont(str(font_path))
            chosen_font = font_manager.FontProperties(fname=str(font_path)).get_name()
            break
    if chosen_font:
        plt.rcParams["font.family"] = chosen_font
    else:
        warnings.warn(
            (
                "未找到可用的中文字体，图表中文字可能会出现乱码。"
                "请安装微软雅黑/黑体或将字体文件放到 resources/fonts/ 目录。"
            ),
            UserWarning,
            stacklevel=2,
        )
    plt.rcParams["axes.unicode_minus"] = False
    return chosen_font


def _generate_label_offsets(count: int, base_radius: float = 18) -> list[tuple[float, float]]:
    """根据数量生成均匀分布的标签偏移量，避免注释互相遮挡。"""

    if count <= 0:
        return []
    offsets: list[tuple[float, float]] = []
    for idx in range(count):
        angle = 2 * math.pi * idx / count
        dx = base_radius * math.cos(angle)
        dy = base_radius * math.sin(angle)
        offsets.append((dx, dy))
    return offsets


def _annotate_highlight_points(ax, highlights: list[HighlightPoint]) -> None:
    """在图表中绘制带箭头注释的高亮点。"""

    if not highlights:
        return
    offsets = _generate_label_offsets(len(highlights))
    for (label, toxicity, aci), (dx, dy) in zip(highlights, offsets):
        ha = "left" if dx >= 0 else "right"
        va = "bottom" if dy >= 0 else "top"
        ax.annotate(
            label,
            xy=(toxicity, aci),
            xycoords="data",
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=9,
            ha=ha,
            va=va,
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#FF6B6B", lw=0.8, alpha=0.85),
            arrowprops=dict(arrowstyle="->", color="#FF6B6B", lw=0.8, shrinkA=2, shrinkB=2),
        )


def plot_nsga_front(
    metrics: np.ndarray,
    highlights: Iterable[HighlightPoint],
    output_path: Path,
    *,
    dpi: int = 200,
    colormap: str = "viridis",
) -> None:
    """绘制 NSGA-II 帕累托前沿图。"""

    metrics = np.asarray(metrics)
    if not metrics.size:
        return
    aci_scores = [-row[0] for row in metrics]
    toxicities = [row[1] for row in metrics]
    fig, ax = plt.subplots(figsize=(6, 4))
    scatter = ax.scatter(
        toxicities,
        aci_scores,
        c=aci_scores,
        cmap=colormap,
        edgecolors="black",
        alpha=0.8,
    )
    ax.set_xlabel("肝毒性（越低越好）")
    ax.set_ylabel("ACI（越高越好）")
    ax.set_title("NSGA-II 帕累托前沿")
    fig.colorbar(scatter, ax=ax, label="ACI")
    highlight_points = list(highlights or [])
    if highlight_points:
        ax.scatter(
            [toxicity for _, toxicity, _ in highlight_points],
            [aci for _, _, aci in highlight_points],
            color="#FF6B6B",
            marker="*",
            s=140,
            edgecolors="white",
            linewidths=1.2,
            zorder=5,
        )
        _annotate_highlight_points(ax, highlight_points)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_candidate_mix_bar(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    output_path: Path,
    *,
    dpi: int = 200,
) -> None:
    """绘制最佳候选的条形占比图。"""

    normalized = candidate.with_normalized()
    names: list[str] = []
    values: list[float] = []
    for idx in normalized.iter_selected_indices():
        names.append(ingredients[idx].name)
        values.append(normalized.proportions[idx] * 100)
    if not names:
        return
    fig, ax = plt.subplots(figsize=(5, 3 + 0.2 * len(names)))
    bars = ax.barh(names, values, color="#6495ED")
    ax.set_xlabel("占比 (%)")
    ax.set_title("最佳候选配比")
    for bar, value in zip(bars, values):
        ax.text(value + 1, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%", va="center")
    ax.set_xlim(0, max(100, max(values) + 10))
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_candidate_mix_pie(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    output_path: Path,
    *,
    dpi: int = 200,
) -> None:
    """绘制最佳候选的饼图占比。"""

    normalized = candidate.with_normalized()
    labels: list[str] = []
    values: list[float] = []
    for idx in normalized.iter_selected_indices():
        labels.append(ingredients[idx].name)
        values.append(normalized.proportions[idx] * 100)
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.pie(values, labels=labels, autopct="%1.1f%%", startangle=90, pctdistance=0.8)
    ax.set_title("最佳候选配比（饼图）")
    ax.axis("equal")
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def render_nsga_visualizations(
    best_candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    metrics: np.ndarray,
    highlight_points: Iterable[HighlightPoint] | None = None,
    output_dir: Path | None = None,
    *,
    dpi: int = 200,
    colormap: str = "viridis",
) -> dict[str, Path]:
    """统一生成 NSGA 相关图表并返回路径。"""

    target_dir = ensure_output_dir(output_dir)
    pareto_path = target_dir / "nsga_pareto.png"
    bar_path = target_dir / "best_mix_bar.png"
    pie_path = target_dir / "best_mix_pie.png"
    plot_nsga_front(metrics, list(highlight_points or []), pareto_path, dpi=dpi, colormap=colormap)
    plot_candidate_mix_bar(best_candidate, ingredients, bar_path, dpi=dpi)
    plot_candidate_mix_pie(best_candidate, ingredients, pie_path, dpi=dpi)
    return {
        "pareto": pareto_path,
        "mix_bar": bar_path,
        "mix_pie": pie_path,
    }
