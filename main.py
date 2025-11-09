from __future__ import annotations

import argparse
import sys

# 说明：
# - 将重依赖（pymoo/pyswarms/matplotlib 等）的导入延迟到函数内部，
#   使得项目在未安装这些依赖的环境中仍可默认运行（hello 模式）。
# - 若使用 demo 模式，则按需导入并运行原有优化与可视化流程。


def sample_ingredients() -> list["Ingredient"]:
    """返回示例成分列表，用于算法演示。"""

    # 延迟导入，避免默认模式下的依赖报错
    from algorithms import Ingredient

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


def describe_candidate(label: str, candidate: "CandidateSolution", ingredients: list["Ingredient"]) -> None:
    """打印候选结果的指标与具体成分选取。

    Args:
        label: 输出标签，如算法名称。
        candidate: 要展示的组合。
        ingredients: 成分清单。
    """
    # 延迟导入以避免默认模式触发依赖
    from algorithms import evaluate_metrics

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

def run_single_objective_algorithms(ingredients: list["Ingredient"]) -> tuple[list["CandidateSolution"], float]:
    """执行 pymoo GA 与 pyswarms PSO，并返回两个最优解用于对比。

    Args:
        ingredients: 成分列表。

    Returns:
        [GA 最优, PSO 最优]。
    """
    # 延迟导入
    from algorithms import PymooSingleObjectiveGA, PySwarmsPSO

    ga = PymooSingleObjectiveGA(ingredients)
    pso = PySwarmsPSO(ingredients)
    # 两种算法分别返回单目标最优解，用以展示不同启发式在相同指标下的表现差异。
    best_ga = ga.run()
    best_pso = pso.run()
    return [best_ga, best_pso], ga.toxicity_weight
def run_nsga(
    ingredients: list["Ingredient"],
    toxicity_weight: float,
    highlight_points: list[tuple[str, float, float]] | None = None,
) -> None:
    """运行 NSGA-II 并打印非支配集候选的指标，并按权重挑选一个候选。

    Args:
        ingredients: 成分列表。

    Returns:
        None
    """
    # 延迟导入
    from algorithms import PymooNSGAII, evaluate_metrics, select_candidate_by_weighted_sum
    from visualization import render_nsga_visualizations

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
    highlight_points = list(highlight_points or [])
    highlight_points.append(("NSGA-权重", best_metrics["toxicity"], best_metrics["aci_raw"]))
    artifacts = render_nsga_visualizations(best_candidate, ingredients, metrics, highlight_points)
    print(
        f"\n可视化已保存：{artifacts['pareto']}、{artifacts['mix_bar']}、{artifacts['mix_pie']}"
    )

def run_demo() -> int:
    """运行原演示（包含算法与可视化）。"""

    try:
        from visualization import configure_matplotlib_fonts
    except Exception as e:  # noqa: BLE001 - 直接面向用户的报错提示
        print(
            "缺少可视化或相关依赖，无法运行 demo。\n"
            "请在可写环境下安装依赖后重试（例如：pymoo、pyswarms、matplotlib 等）",
            file=sys.stderr,
        )
        return 2

    configure_matplotlib_fonts()

    ingredients = sample_ingredients()
    single_objective_results, toxicity_weight = run_single_objective_algorithms(ingredients)
    print("单目标算法结果：")
    for label, candidate in zip(("遗传算法", "粒子群算法"), single_objective_results):
        describe_candidate(label, candidate, ingredients)

    highlight_points: list[tuple[str, float, float]] = []
    # 延迟导入评估函数，避免默认模式触发依赖
    from algorithms import evaluate_metrics

    for label, candidate in zip(("GA 最优", "PSO 最优"), single_objective_results):
        _, toxicity, _, aci_raw = evaluate_metrics(candidate, ingredients)
        highlight_points.append((label, toxicity, aci_raw))
    # 单目标权重会作为多目标加权依据，形成闭环体验。
    run_nsga(ingredients, toxicity_weight, highlight_points)
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI 入口：
    - 默认 hello 模式：打印最小可运行输出，满足快速自检。
    - demo 模式：运行原有算法演示（需要额外依赖与可写目录）。
    """

    parser = argparse.ArgumentParser(description="noname CLI")
    parser.add_argument(
        "--mode",
        choices=("hello", "demo"),
        default="hello",
        help="运行模式：hello（默认）或 demo（需要依赖）",
    )
    args = parser.parse_args(argv)

    if args.mode == "hello":
        print("Hello from noname!")
        return 0

    # demo 模式
    return run_demo()


if __name__ == "__main__":  # pragma: no cover - 直接作为脚本运行
    raise SystemExit(main())
