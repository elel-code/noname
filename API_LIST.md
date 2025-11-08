# API 文档

本文档提供了 `algorithms.py` 和 `main.py` 中主要类和函数的详细说明。

---

## `algorithms.py` - 核心算法与评估

### 数据结构

-   **`Ingredient`**: 一个 `dataclass`，用于存储单个成分的所有相关信息。
    -   **字段**: `name`, `host_drug`, `hepatotoxicity_score`, `synergy_baseline`, 以及多个 ADME 特征（`MW`, `AlogP`, 等）。
    -   **作用**: 作为整个系统的基础数据单元。

-   **`CandidateSolution`**: 表示一个候选的成分组合。
    -   **结构**: 内部包含一个布尔列表 `select_bits` (决定哪些成分被选中) 和一个浮点数数组 `proportions` (存储各成分的配比)。
    -   **主要方法**:
        -   `from_components()`: 从指定的成分和配比创建实例。
        -   `iter_selects()`: 迭代选中的成分及其配比。
        -   `normalized_proportions()`: 返回归一化后的配比。

### ACI 与肝毒性评估函数

这些函数接收一个 `CandidateSolution` 和成分列表，计算特定的性能指标。

-   `compute_desirability_score(ingredient)`: 计算单个成分的“可用性”分数。
-   `compute_complementarity_score(candidate, ingredients)`: 计算组合的“互补性”分数。
-   `compute_aci_score(candidate, ...)`: 结合可用性和互补性，计算最终的 ACI 分数。
-   `compute_hepatotoxicity_score(candidate, ...)`: 计算组合的加权肝毒性分数。

### 优化目标函数

这些函数定义了遗传算法和粒子群优化的目标。

-   **`single_objective_score(candidate, ...)`**: 单目标优化函数。
    -   **公式**: `ACI - (toxicity_weight * 肝毒性)`
    -   **用途**: 用于 GA 和 PSO。

-   **`multi_objective_score(candidate, ...)`**: 多目标优化函数。
    -   **返回**: `(-ACI, 肝毒性)`
    -   **注意**: 返回 `-ACI` 是因为 `pymoo` 默认最小化所有目标。
    -   **用途**: 用于 NSGA-II。

### `pymoo` 与 `pyswarms` 封装

-   **`CombinationProblem`**: `pymoo` 的 `Problem` 子类，将上述评估函数与优化框架连接起来。

-   **`PymooSingleObjectiveGA`**: 使用 `pymoo` 的遗传算法寻找单目标最优解。
-   **`PymooNSGAII`**: 使用 `pymoo` 的 NSGA-II 算法寻找多目标帕累托前沿。
-   **`PySwarmsPSO`**: 使用 `pyswarms` 的粒子群算法寻找单目标最优解。

---

## `main.py` - 程序入口与流程控制

### 核心执行函数

-   **`run_single_objective_algorithms(ingredients)`**:
    -   **流程**:
        1.  初始化并运行 `PymooSingleObjectiveGA`。
        2.  初始化并运行 `PySwarmsPSO`。
        3.  打印两种算法找到的最佳候选解。

-   **`run_nsga(ingredients, toxicity_weight)`**:
    -   **流程**:
        1.  初始化并运行 `PymooNSGAII`，得到帕累托前沿上的多个候选解。
        2.  打印所有非支配解。
        3.  使用 `select_candidate_by_weighted_sum` 从帕累托前沿中选择一个“折衷”最优解。

### 辅助函数

-   `sample_ingredients()`: 提供一组硬编码的示例成分数据，用于快速测试。
-   `describe_candidate(label, candidate, ingredients)`: 一个格式化打印函数，用于清晰地展示候选解的各项指标和成分配比。

### 程序主入口

-   **`main()`**:
    1.  调用 `sample_ingredients()` 加载数据。
    2.  执行 `run_single_objective_algorithms`。
    3.  执行 `run_nsga`。
