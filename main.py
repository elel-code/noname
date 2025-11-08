from algorithms import (
    CandidateSolution,
    Ingredient,
    PymooNSGAII,
    PymooSingleObjectiveGA,
    PySwarmsPSO,
    evaluate_metrics,
)


def sample_ingredients() -> list[Ingredient]:
    """返回示例成分列表，用于算法演示。"""
    return [
        Ingredient("黄芩素", "黄芩", hepatotoxicity_score=0.3, synergy_baseline=0.9, ob=35.0, dl=0.18),
        Ingredient("川芎嗪", "川芎", hepatotoxicity_score=0.2, synergy_baseline=0.85, ob=32.0, dl=0.09),
        Ingredient("黄连碱", "黄连", hepatotoxicity_score=0.5, synergy_baseline=0.95, ob=26.0, dl=0.12),
        Ingredient("丹参酮", "丹参", hepatotoxicity_score=0.15, synergy_baseline=0.7, ob=45.0, dl=0.26),
        Ingredient("甘草酸", "甘草", hepatotoxicity_score=0.4, synergy_baseline=0.75, ob=28.0, dl=0.11),
        Ingredient("白芍苷", "白芍", hepatotoxicity_score=0.18, synergy_baseline=0.68, ob=38.0, dl=0.14),
    ]


def describe_candidate(label: str, candidate: CandidateSolution, ingredients: list[Ingredient]) -> None:
    """打印候选结果的三项指标与具体成分选取。

    Args:
        label: 输出标签，如算法名称。
        candidate: 要展示的组合。
        ingredients: 成分清单。
    """
    synergy, toxicity, penalty = evaluate_metrics(candidate, ingredients)
    selected = [ing.name for select, ing in zip(candidate.selects, ingredients) if select]
    print(f"{label}：协同 {synergy:.3f}，肝毒性 {toxicity:.3f}，惩罚 {penalty:.2f}，配伍 {selected}")


def run_single_objective_algorithms(ingredients: list[Ingredient]) -> list[CandidateSolution]:
    """执行 pymoo GA 与 pyswarms PSO，并返回两个最优解用于对比。

    Args:
        ingredients: 成分列表。

    Returns:
        [GA 最优, PSO 最优]。
    """
    ga = PymooSingleObjectiveGA(ingredients)
    pso = PySwarmsPSO(ingredients)
    best_ga = ga.run()
    best_pso = pso.run()
    return [best_ga, best_pso]


def run_nsga(ingredients: list[Ingredient]) -> None:
    """运行 NSGA-II 并打印非支配集候选的指标。

    Args:
        ingredients: 成分列表。

    Returns:
        None
    """
    nsga = PymooNSGAII(ingredients)
    solutions, metrics = nsga.run()
    print("\nNSGA-II 非支配集候选：")
    for idx, (candidate, metric) in enumerate(zip(solutions, metrics), 1):
        synergy = -metric[0]
        toxicity = metric[1]
        penalty = evaluate_metrics(candidate, ingredients)[2]
        print(f"  {idx}. 协同 {synergy:.3f}，肝毒性 {toxicity:.3f}，惩罚 {penalty:.2f}")


def main() -> None:
    """主执行逻辑：分别演示单目标与 NSGA-II 多目标输出。

    Returns:
        None
    """
    ingredients = sample_ingredients()
    single_objective_results = run_single_objective_algorithms(ingredients)
    print("单目标算法结果：")
    for label, candidate in zip(("遗传算法", "粒子群算法"), single_objective_results):
        describe_candidate(label, candidate, ingredients)
    run_nsga(ingredients)


if __name__ == "__main__":
    main()
