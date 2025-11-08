# API 列表

> 配置：`config/algorithms.json` 软编码 GA/NSGA-II/PSO 的 generations、population、迭代次数及 PSO options，`algorithms.py` 会在 `_load_algorithm_config()` 中读取（缺省使用内置 `_DEFAULT_ALGO_CONFIG`）。

## algorithms.py
- `Ingredient`：冻结数据类，存储成分名、母药名、SwissADME 肝毒性分数、协同基线以及 ADME 特征（MW、ALogP、HDon/HAcc、OB、Caco-2、BBB、DL、FASA-、HL）。
- `CandidateSolution`：混合染色体，使用位串 `select_bits` + `array('f')` 存配比，提供 `from_components`、`iter_selects`、`iter_selected_indices`、`normalized_proportions`、`with_normalized`。
- `ensure_selection(selects, ingredients)`：保证选择向量至少包含一个 `True`；若全 False，则启用协同基线最高的成分以保持可复现。
- `vector_to_candidate(vector, ingredients)`：将 `[0,1]` 向量截取为选择+配比（只归一一次），返回 `CandidateSolution`。
- `compute_desirability_score(ingredient)`：将单体 ADME 特征映射到 0-1 并取几何平均。
- `compute_complementarity_score(candidate, ingredients)`：基于特征方差衡量组合互补性。
- `compute_aci_score(candidate, ingredients)`：融合单体可用性与互补项，生成 ACI。
- `compute_hepatotoxicity_score(candidate, ingredients)`：汇总加权的 hepatotoxicity_score。
- `diversity_penalty(candidate, ingredients)`：若某味药占比 ≥80%，返回 0.1；否则 0。
- `evaluate_metrics(candidate, ingredients)`：一次性返回 ACI、肝毒性、惩罚三要素（内部进行单次归一）。
- `single_objective_score(candidate, ingredients, toxicity_weight=1.0)`：ACI 减去毒性的单目标适应度。
- `multi_objective_score(candidate, ingredients)`：返回 `(ACI, 肝毒性)`，供 NSGA-II 排序。
- `_load_algorithm_config()`：读取/缓存 `config/algorithms.json`（缺失退回默认配置）。
- `CombinationProblem`：`pymoo` 的 `ElementwiseProblem` 实现，支持 single/multi 目标。
- `PymooSingleObjectiveGA`：调用 `pymoo.algorithms.soo.nonconvex.ga.GA` 处理 ACI-毒性优化，默认参数来自配置。
- `PymooNSGAII`：调用 `pymoo.algorithms.moo.nsga2.NSGA2` 生成 Pareto 解（最大化 ACI、最小化毒性），默认 generations/pop_size 由配置驱动。
- `PySwarmsPSO`：使用 `pyswarms.single.global_best.GlobalBestPSO` 求单目标最优，迭代数/粒子数/`options` 由配置文件注入。

## main.py
- `sample_ingredients()`：返回示例成分列表（黄芩素、川芎嗪、黄连碱等）。
- `describe_candidate(label, candidate, ingredients)`：打印 ACI、毒性、惩罚与选中成分（通过 `iter_selected_indices` 获取）。
- `run_single_objective_algorithms(ingredients)`：执行 `PymooSingleObjectiveGA` 与 `PySwarmsPSO` 并返回两个最优解。
- `run_nsga(ingredients)`：运行 `PymooNSGAII` 并打印非支配集候选的 ACI/毒性/惩罚。
- `main()`：入口，读取示例成分、输出单目标与 NSGA-II 结果。
