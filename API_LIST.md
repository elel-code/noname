# API 列表

## algorithms.py
- `Ingredient`：冻结数据类，存储成分名、母药名、SwissADME 肝毒性分数、协同基线、OB/DL。
- `CandidateSolution`：混合染色体，使用位串 `select_bits` + `array('f')` 存配比，提供 `from_components`、`iter_selects`、`iter_selected_indices`、`normalized_proportions`、`with_normalized`。
- `ensure_selection(selects, ingredients)`：保证选择向量至少包含一个 `True`；若全 False，则启用协同基线最高的成分以保持可复现。
- `vector_to_candidate(vector, ingredients)`：将 `[0,1]` 向量截取为选择+配比（只归一一次），返回 `CandidateSolution`。
- `compute_loew_synergy(candidate, ingredients)`：在归一化组合上计算 Loewe 协同加权得分。
- `compute_hepatotoxicity_score(candidate, ingredients)`：汇总加权的 hepatotoxicity_score。
- `diversity_penalty(candidate, ingredients)`：若某味药占比 ≥80%，返回 0.1；否则 0。
- `evaluate_metrics(candidate, ingredients)`：一次性返回协同、肝毒性、惩罚三要素（内部进行单次归一）。
- `single_objective_score(candidate, ingredients, toxicity_weight=1.0)`：协同减去毒性的单目标适应度。
- `multi_objective_score(candidate, ingredients)`：返回 `(协同, 肝毒性)`，供 NSGA-II 排序。
- `CombinationProblem`：`pymoo` 的 `ElementwiseProblem` 实现，支持 single/multi 目标。
- `PymooSingleObjectiveGA`：调用 `pymoo.algorithms.soo.nonconvex.ga.GA` 处理协同-毒性最小化。
- `PymooNSGAII`：调用 `pymoo.algorithms.moo.nsga2.NSGA2` 生成 Pareto 解。
- `PySwarmsPSO`：使用 `pyswarms.single.global_best.GlobalBestPSO` 求单目标最优。

## main.py
- `sample_ingredients()`：返回示例成分列表（黄芩素、川芎嗪、黄连碱等）。
- `describe_candidate(label, candidate, ingredients)`：打印协同、毒性、惩罚与选中成分（通过 `iter_selected_indices` 获取）。
- `run_single_objective_algorithms(ingredients)`：执行 `PymooSingleObjectiveGA` 与 `PySwarmsPSO` 并返回两个最优解。
- `run_nsga(ingredients)`：运行 `PymooNSGAII` 并打印非支配集候选的协同/毒性/惩罚。
- `main()`：入口，读取示例成分、输出单目标与 NSGA-II 结果。
