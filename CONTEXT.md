# 上下文概述

1. 当前仓库是一个简化的 Python 模板，主入口 `main.py` 调用基于 `pymoo` + `pyswarms` 的优化插件，展示中药成分复方的单目标与多目标选配流程。
2. 核心算法集中在 `algorithms.py`：`Ingredient`/`CandidateSolution` 表征混合染色体，后者使用位串（bytearray）存储选择掩码 + `array('f')` 存配比，以降低成千上万成分时的内存占用；`vector_to_candidate` 将 `[0,1]` 向量映射为选中状态与配比，评估流程计算基于 ADME 的 ACI（含 coverage 惩罚与互补项配置）、肝毒性与基于配比的多样性惩罚。
3. 单目标通过 `PymooSingleObjectiveGA`（GA）和 `PySwarmsPSO`（PSO）运行，返回 ACI-毒性的综合得分；多目标通过 `PymooNSGAII`（NSGA-II）输出 Pareto 解并记录 ACI/毒性指标，同时可用 `select_candidate_by_weighted_sum` 在帕累托前沿中按权重挑选解（惩罚项已融入适应度）。
4. 所有算法依赖 `numpy>=1.26`、`pymoo>=0.6.1.5,<0.7`、`pyswarms>=1.3.0`，已在 `.venv` 成功安装并通过 `python main.py` 验证；`pyproject.toml` 记录了这些依赖。
5. `ensure_selection` 采用确定性策略：若初始染色体全部为 False，就启用协同基线最高的成分保证可复现。
6. 算法参数（GA/NSGA-II/PSO 的代数、种群、迭代次数、超参数）现集中在 `config/algorithms.json` 中进行软编码；若文件缺失则回退至 `_DEFAULT_ALGO_CONFIG`。
