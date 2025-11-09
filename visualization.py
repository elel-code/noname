from __future__ import annotations

"""可视化模块

- 负责单目标（GA/PSO）与多目标（NSGA-II）的图表生成；
- 封装统一的输出目录创建、字体配置与多格式保存；
- 设计成按需调用，便于在无图形环境下禁用（--no-fig）。
"""

import math
import warnings
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from algorithms import CandidateSolution, Ingredient, evaluate_metrics

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
    "render_single_objective_visualizations",
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
    output_path: Path | None,
    *,
    dpi: int = 200,
    colormap: str = "viridis",
) -> plt.Figure | None:
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
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def plot_candidate_mix_bar(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    output_path: Path | None,
    *,
    dpi: int = 200,
    responsive: bool = True,
) -> plt.Figure | None:
    """绘制最佳候选的条形占比图。"""

    normalized = candidate.with_normalized()
    names: list[str] = []
    values: list[float] = []
    for idx in normalized.iter_selected_indices():
        names.append(ingredients[idx].name)
        values.append(normalized.proportions[idx] * 100)
    if not names:
        return
    # 响应式高度：随项数变化；固定模式下采用较小增长系数
    growth = 0.25 if responsive else 0.12
    height = 3 + growth * len(names)
    fig, ax = plt.subplots(figsize=(5.6, height))
    bars = ax.barh(names, values, color="#6495ED")
    ax.set_xlabel("占比 (%)")
    ax.set_title("最佳候选配比")
    for bar, value in zip(bars, values):
        ax.text(value + 1, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%", va="center")
    ax.set_xlim(0, max(100, max(values) + 10))
    fig.tight_layout()
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def plot_candidate_mix_pie(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    output_path: Path | None,
    *,
    dpi: int = 200,
) -> plt.Figure | None:
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
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def plot_nsga_density(
    metrics: np.ndarray,
    output_path: Path | None,
    *,
    dpi: int = 200,
    colormap: str = "viridis",
) -> plt.Figure | None:
    """绘制 NSGA-II 解的二维密度（hexbin）。"""

    metrics = np.asarray(metrics)
    if not metrics.size:
        return None
    aci_scores = [-row[0] for row in metrics]
    toxicities = [row[1] for row in metrics]
    fig, ax = plt.subplots(figsize=(6, 4))
    hb = ax.hexbin(toxicities, aci_scores, gridsize=25, cmap=colormap, mincnt=1)
    ax.set_xlabel("肝毒性（越低越好）")
    ax.set_ylabel("ACI（越高越好）")
    ax.set_title("NSGA-II 解密度")
    fig.colorbar(hb, ax=ax, label="计数")
    fig.tight_layout()
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def plot_metrics_hist(
    metrics: np.ndarray,
    output_path: Path | None,
    *,
    dpi: int = 200,
) -> plt.Figure | None:
    """绘制 ACI 与肝毒性分布直方图。"""

    metrics = np.asarray(metrics)
    if not metrics.size:
        return None
    aci_scores = np.array([-row[0] for row in metrics])
    toxicities = np.array([row[1] for row in metrics])
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    axes[0].hist(aci_scores, bins=20, color="#4C9AFF", alpha=0.85)
    axes[0].set_title("ACI 分布")
    axes[0].set_xlabel("ACI")
    axes[0].set_ylabel("频数")
    axes[1].hist(toxicities, bins=20, color="#FFB020", alpha=0.85)
    axes[1].set_title("肝毒性分布")
    axes[1].set_xlabel("肝毒性")
    axes[1].set_ylabel("频数")
    fig.tight_layout()
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def plot_candidate_radar(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    output_path: Path | None,
    *,
    dpi: int = 200,
) -> plt.Figure | None:
    """绘制候选解综合指标雷达图（归一化展示）。"""

    # 计算原始指标
    aci_raw, toxicity, penalty, _ = evaluate_metrics(candidate, ingredients)
    # 归一化到 0..1（粗略、可视对比用）
    aci_norm = (aci_raw + 1.0) / 2.0  # 假设 ACI 原始在约 [-1,1]
    tox_norm = 1.0 - min(max(toxicity, 0.0), 1.0)
    pen_norm = 1.0 - min(max(penalty, 0.0), 1.0)
    labels = ["ACI", "(反)毒性", "(反)惩罚"]
    values = [aci_norm, tox_norm, pen_norm]
    values += values[:1]
    angles = np.linspace(0, 2 * np.pi, len(labels) + 1)
    fig = plt.figure(figsize=(4.8, 4.8))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(angles, values, color="#3FB27F", linewidth=2)
    ax.fill(angles, values, color="#3FB27F", alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_yticklabels([])
    ax.set_title("候选解雷达图")
    fig.tight_layout()
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def _savefig_multi(fig: plt.Figure, base: Path, formats: Sequence[str], *, dpi: int) -> list[Path]:
    """以多格式保存图像，返回保存路径列表（不关闭 fig）。"""
    saved: list[Path] = []
    for fmt in formats:
        path = base.with_suffix("." + fmt.lower())
        fig.savefig(path, dpi=dpi)
        saved.append(path)
    plt.close(fig)
    return saved


def _apply_title_prefix(fig: plt.Figure, prefix: str | None) -> None:
    """为图中所有子图标题添加前缀。"""
    if not prefix:
        return
    for ax in fig.axes:
        title = ax.get_title()
        if title:
            ax.set_title(f"{prefix}{title}")


def plot_convergence(
    series: Sequence[float] | None,
    label: str,
    output_path: Path | None,
    *,
    dpi: int = 200,
) -> plt.Figure | None:
    """绘制单目标优化收敛曲线（若数据不足则跳过）。"""
    if not series or len(series) < 2:
        return None
    xs = np.arange(1, len(series) + 1)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(xs, series, color="#2E7D32", linewidth=1.8)
    ax.set_xlabel("迭代/世代")
    ax.set_ylabel("目标值（ACI - w*毒性）")
    ax.set_title(f"{label} 收敛曲线")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if output_path is None:
        return fig
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return None


def render_nsga_visualizations(
    best_candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    metrics: np.ndarray,
    highlight_points: Iterable[HighlightPoint] | None = None,
    output_dir: Path | None = None,
    *,
    dpi: int = 200,
    colormap: str = "viridis",
    figs: Sequence[str] | None = None,
    formats: Sequence[str] | None = None,
    annotate: bool = True,
    responsive: bool = True,
    title_prefix: str | None = None,
) -> dict[str, Path]:
    """统一生成 NSGA 相关图表并返回路径。"""
    target_dir = ensure_output_dir(output_dir)
    all_supported: tuple[str, ...] = ("pareto", "bar", "pie", "hist", "density", "radar")
    figs_to_draw = tuple([f for f in (figs or ("pareto", "bar", "pie")) if f in all_supported])
    results: dict[str, Path] = {}
    save_formats = list(formats or ["png"])

    # 1) 帕累托前沿
    if "pareto" in figs_to_draw:
        effective_highlights = list(highlight_points or []) if annotate else []
        fig = plot_nsga_front(metrics, effective_highlights, None, dpi=dpi, colormap=colormap)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "nsga_pareto", save_formats, dpi=dpi)
            if saved:
                results["pareto"] = saved[0]

    # 2) 候选配比（条形）
    if "bar" in figs_to_draw:
        fig = plot_candidate_mix_bar(best_candidate, ingredients, None, dpi=dpi, responsive=responsive)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "best_mix_bar", save_formats, dpi=dpi)
            if saved:
                results["mix_bar"] = saved[0]

    # 3) 候选配比（饼图）
    if "pie" in figs_to_draw:
        fig = plot_candidate_mix_pie(best_candidate, ingredients, None, dpi=dpi)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "best_mix_pie", save_formats, dpi=dpi)
            if saved:
                results["mix_pie"] = saved[0]

    # 4) 直方图
    if "hist" in figs_to_draw:
        fig = plot_metrics_hist(metrics, None, dpi=dpi)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "metrics_hist", save_formats, dpi=dpi)
            if saved:
                results["hist"] = saved[0]

    # 5) 密度图
    if "density" in figs_to_draw:
        fig = plot_nsga_density(metrics, None, dpi=dpi, colormap=colormap)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "nsga_density", save_formats, dpi=dpi)
            if saved:
                results["density"] = saved[0]

    # 6) 雷达图
    if "radar" in figs_to_draw:
        fig = plot_candidate_radar(best_candidate, ingredients, None, dpi=dpi)
        if fig is not None:
            _apply_title_prefix(fig, title_prefix)
            saved = _savefig_multi(fig, target_dir / "candidate_radar", save_formats, dpi=dpi)
            if saved:
                results["radar"] = saved[0]

    return results


def render_single_objective_visualizations(
    best_ga: CandidateSolution,
    best_pso: CandidateSolution,
    ingredients: Sequence[Ingredient],
    ga_convergence: Sequence[float] | None,
    pso_convergence: Sequence[float] | None,
    output_dir: Path | None = None,
    *,
    dpi: int = 200,
    formats: Sequence[str] | None = None,
    title_prefix: str | None = None,
) -> dict[str, Path]:
    """为单目标 GA/PSO 生成独立可视化：配比图 + 收敛曲线。"""
    target_dir = ensure_output_dir(output_dir)
    save_formats = list(formats or ["png"])
    results: dict[str, Path] = {}

    # GA 配比
    fig = plot_candidate_mix_bar(best_ga, ingredients, None, dpi=dpi, responsive=True)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "ga_best_mix_bar", save_formats, dpi=dpi)
        if saved:
            results["ga_mix_bar"] = saved[0]
    fig = plot_candidate_mix_pie(best_ga, ingredients, None, dpi=dpi)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "ga_best_mix_pie", save_formats, dpi=dpi)
        if saved:
            results["ga_mix_pie"] = saved[0]

    # PSO 配比
    fig = plot_candidate_mix_bar(best_pso, ingredients, None, dpi=dpi, responsive=True)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "pso_best_mix_bar", save_formats, dpi=dpi)
        if saved:
            results["pso_mix_bar"] = saved[0]
    fig = plot_candidate_mix_pie(best_pso, ingredients, None, dpi=dpi)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "pso_best_mix_pie", save_formats, dpi=dpi)
        if saved:
            results["pso_mix_pie"] = saved[0]

    # 收敛曲线
    fig = plot_convergence(list(ga_convergence or []), "GA", None, dpi=dpi)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "ga_convergence", save_formats, dpi=dpi)
        if saved:
            results["ga_convergence"] = saved[0]
    fig = plot_convergence(list(pso_convergence or []), "PSO", None, dpi=dpi)
    if fig is not None:
        _apply_title_prefix(fig, title_prefix)
        saved = _savefig_multi(fig, target_dir / "pso_convergence", save_formats, dpi=dpi)
        if saved:
            results["pso_convergence"] = saved[0]

    return results
