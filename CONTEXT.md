# 上下文概述

1. 当前仓库是一个简化的 Python 模板，主入口 `main.py` 调用基于 `pymoo` + `pyswarms` 的优化插件，展示中药成分复方的单目标与多目标选配流程。
2. 核心算法集中在 `algorithms.py`：`Ingredient`/`CandidateSolution` 表征混合染色体，`vector_to_candidate` 将 `[0,1]` 向量映射为选中状态与配比；`CandidateSolution` 现使用位串（bytearray）存储选择掩码、`array('f')` 存储配比以降低 1 万级成分时的内存压力，并在入口归一一次即可。
3. 单目标通过 `PymooSingleObjectiveGA`（GA）和 `PySwarmsPSO`（PSO）运行，返回协同-毒性的综合得分；多目标通过 `PymooNSGAII`（NSGA-II）输出 Pareto 解并记录协同/毒性指标（多样性惩罚按配比统计）。
4. 所有算法依赖 `numpy>=1.26`、`pymoo>=0.6.1.5,<0.7`、`pyswarms>=1.3.0`，已在 `.venv` 成功安装并通过 `python main.py` 验证；`pyproject.toml` 记录了这些依赖。
5. `ensure_selection` 采用确定性策略：若初始染色体全部为 False，就启用协同基线最高的成分，保证大规模数据场景下的可复现行为。
