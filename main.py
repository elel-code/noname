from pathlib import Path
import math
import warnings

import matplotlib.pyplot as plt
from matplotlib import font_manager

from algorithms import (
    CandidateSolution,
    Ingredient,
    PymooNSGAII,
    PymooSingleObjectiveGA,
    PySwarmsPSO,
    evaluate_metrics,
    select_candidate_by_weighted_sum,
)


def sample_ingredients() -> list[Ingredient]:
    """返回示例成分列表，用于算法演示。"""
    return [
        Ingredient(
            "黄芩素",
            smiles="OC1=CC=C2C(=C1)C(=O)C=C(O)C(=O)O2",
            hepatotoxicity_score=0.3,
            ob=35.0,
            dl=0.18,
            mw=270.2,
            alogp=2.1,
            h_don=3,
            h_acc=5,
            caco2=0.82,
            bbb=-0.3,
            fasa=0.33,
            hl=6.5,
        ),
        Ingredient(
            "川芎嗪",
            smiles="CN1C=NC(=CN=C1N(C)C)N(C)C",
            hepatotoxicity_score=0.2,
            ob=32.0,
            dl=0.09,
            mw=136.2,
            alogp=0.3,
            h_don=1,
            h_acc=2,
            caco2=1.1,
            bbb=0.2,
            fasa=0.28,
            hl=3.5,
        ),
        Ingredient(
            "黄连碱",
            smiles="COC1=CC=CC=C1O",
            ob=26.0,
            dl=0.12,
            mw=336.4,
            alogp=1.5,
            h_don=0,
            h_acc=4,
            caco2=0.4,
            bbb=-0.7,
            fasa=0.48,
            hl=8.0,
        ),
        Ingredient(
            "丹参酮",
            smiles="CC1=C2C(=O)C=C(C2=CC3=CC=CC=C31)C",
            hepatotoxicity_score=0.15,
            ob=45.0,
            dl=0.26,
            mw=294.3,
            alogp=3.1,
            h_don=0,
            h_acc=4,
            caco2=1.2,
            bbb=0.5,
            fasa=0.22,
            hl=10.5,
        ),
        Ingredient(
            "甘草酸",
            smiles="CC1=CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C(=O)O",
            hepatotoxicity_score=0.4,
            ob=28.0,
            dl=0.11,
            mw=822.0,
            alogp=1.2,
            h_don=5,
            h_acc=16,
            caco2=0.2,
            bbb=-1.2,
            fasa=0.62,
            hl=4.0,
        ),
        Ingredient(
            "白芍苷",
            smiles="COC1=CC=C2C(=O)C3=C(C=C(C(=C3O2)O)O)O1",
            hepatotoxicity_score=0.18,
            ob=38.0,
            dl=0.14,
            mw=480.4,
            alogp=-1.1,
            h_don=6,
            h_acc=11,
            caco2=0.35,
            bbb=-0.9,
            fasa=0.58,
            hl=5.5,
        ),
    ]


def describe_candidate(label: str, candidate: CandidateSolution, ingredients: list[Ingredient]) -> None:
    """打印候选结果的指标与具体成分选取。

    Args:
        label: 输出标签，如算法名称。
        candidate: 要展示的组合。
        ingredients: 成分清单。
    """
    aci_penalized, toxicity, penalty, aci_raw = evaluate_metrics(candidate, ingredients)
    normalized = candidate.with_normalized()
    mix = []
    for idx in normalized.iter_selected_indices():
        percentage = normalized.proportions[idx] * 100
        mix.append(f"{ingredients[idx].name}({percentage:.1f}%)")
    mix_repr = "[" + ", ".join(mix) + "]" if mix else "[]"
    # 对外展示时使用 penalized ACI，避免占比过度集中。
    print(
        f"{label}：ACI_pen {aci_penalized:.3f} (raw {aci_raw:.3f})，"
        f"肝毒性 {toxicity:.3f}，惩罚 {penalty:.2f}，配伍 {mix_repr}"
    )


def ensure_output_dir() -> Path:
    output_dir = Path("artifacts")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def configure_matplotlib_fonts() -> None:
    """确保 Matplotlib 在保存图表时可正常渲染中文。"""

    preferred_fonts = (
        "Microsoft YaHei",
        "SimHei",
        "Source Han Sans CN",
        "Noto Sans CJK SC",
        "PingFang SC",
    )
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    chosen_font = next((font for font in preferred_fonts if font in available_fonts), None)
    if chosen_font is None:
        font_dir = Path(__file__).resolve().parent / "resources" / "fonts"
        fallback_files = (
            "SourceHanSansCN-Regular.otf",
            "SourceHanSansSC-Regular.otf",
            "NotoSansSC-Regular.otf",
        )
        for filename in fallback_files:
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
                "未找到可用的中文字体，图表中文字可能仍会出现乱码。"
                "请安装微软雅黑/黑体或将字体放到 resources/fonts/ 下再运行脚本。"
            ),
            UserWarning,
            stacklevel=2,
        )
    plt.rcParams["axes.unicode_minus"] = False


def _generate_label_offsets(count: int, base_radius: float = 18) -> list[tuple[float, float]]:
    """根据数量生成均匀分布的标签偏移量，避免中文注释重叠。"""

    if count <= 0:
        return []
    offsets: list[tuple[float, float]] = []
    for idx in range(count):
        angle = 2 * math.pi * idx / count
        dx = base_radius * math.cos(angle)
        dy = base_radius * math.sin(angle)
        offsets.append((dx, dy))
    return offsets


def annotate_highlight_points(ax, highlights: list[tuple[str, float, float]]) -> None:
    """使用带箭头的注释呈现高亮点，缓解文字拥挤。"""

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


def plot_nsga_front(metrics, highlights: list[tuple[str, float, float]], output_path: Path) -> None:
    if not metrics.size:
        return
    aci_scores = [-row[0] for row in metrics]
    toxicities = [row[1] for row in metrics]
    fig, ax = plt.subplots(figsize=(6, 4))
    scatter = ax.scatter(
        toxicities,
        aci_scores,
        c=aci_scores,
        cmap="viridis",
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
        annotate_highlight_points(ax, highlight_points)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_candidate_mix_bar(candidate: CandidateSolution, ingredients: list[Ingredient], output_path: Path) -> None:
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
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_candidate_mix_pie(candidate: CandidateSolution, ingredients: list[Ingredient], output_path: Path) -> None:
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
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_single_objective_algorithms(ingredients: list[Ingredient]) -> tuple[list[CandidateSolution], float]:
    """执行 pymoo GA 与 pyswarms PSO，并返回两个最优解用于对比。

    Args:
        ingredients: 成分列表。

    Returns:
        [GA 最优, PSO 最优]。
    """
    ga = PymooSingleObjectiveGA(ingredients)
    pso = PySwarmsPSO(ingredients)
    # 两种算法分别返回单目标最优解，用以展示不同启发式在相同指标下的表现差异。
    best_ga = ga.run()
    best_pso = pso.run()
    return [best_ga, best_pso], ga.toxicity_weight


def run_nsga(ingredients: list[Ingredient], toxicity_weight: float, highlight_points: list[tuple[str, float, float]] | None = None) -> None:
    """运行 NSGA-II 并打印非支配集候选的指标，并按权重挑选一个候选。

    Args:
        ingredients: 成分列表。

    Returns:
        None
    """
    nsga = PymooNSGAII(ingredients)
    solutions, metrics = nsga.run()
    print("\nNSGA-II 非支配集候选：")
    for idx, (candidate, metric) in enumerate(zip(solutions, metrics), 1):
        aci = -metric[0]
        toxicity = metric[1]
        penalty = evaluate_metrics(candidate, ingredients)[2]
        # 输出时保留罚分，便于肉眼观察可行/不可行组合差异。
        print(f"  {idx}. ACI {aci:.3f}，肝毒性 {toxicity:.3f}，惩罚 {penalty:.2f}")
    best_candidate, best_metrics = select_candidate_by_weighted_sum(solutions, ingredients, toxicity_weight)
    print(
        f"\n按 toxicity_weight={toxicity_weight:.2f} 的加权结果："
        f"ACI_pen {best_metrics['aci_penalized']:.3f} (raw {best_metrics['aci_raw']:.3f})，"
        f"肝毒性 {best_metrics['toxicity']:.3f}，惩罚 {best_metrics['penalty']:.2f}"
    )
    describe_candidate("NSGA-权重选择", best_candidate, ingredients)
    output_dir = ensure_output_dir()
    pareto_path = output_dir / "nsga_pareto.png"
    mix_path = output_dir / "best_mix_bar.png"
    highlight_points = list(highlight_points or [])
    highlight_points.append(("NSGA-权重", best_metrics["toxicity"], best_metrics["aci_raw"]))
    plot_nsga_front(metrics, highlight_points, pareto_path)
    plot_candidate_mix_bar(best_candidate, ingredients, mix_path)
    pie_path = output_dir / "best_mix_pie.png"
    plot_candidate_mix_pie(best_candidate, ingredients, pie_path)
    print(f"\n可视化已保存：{pareto_path}、{mix_path}、{pie_path}")


def main() -> None:
    """主执行逻辑：分别演示单目标与 NSGA-II 多目标输出。

    Returns:
        None
    """
    ingredients = sample_ingredients()
    single_objective_results, toxicity_weight = run_single_objective_algorithms(ingredients)
    print("单目标算法结果：")
    for label, candidate in zip(("遗传算法", "粒子群算法"), single_objective_results):
        describe_candidate(label, candidate, ingredients)
    highlight_points: list[tuple[str, float, float]] = []
    for label, candidate in zip(("GA 最优", "PSO 最优"), single_objective_results):
        aci, toxicity, _, _ = evaluate_metrics(candidate, ingredients)
        highlight_points.append((label, toxicity, aci))
    # 单目标权重会作为多目标加权依据，形成闭环体验。
    run_nsga(ingredients, toxicity_weight, highlight_points)


configure_matplotlib_fonts()


if __name__ == "__main__":
    main()
