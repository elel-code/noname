from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA as PymooGAAlg
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pyswarms.single.global_best import GlobalBestPSO


@dataclass(frozen=True)
class Ingredient:
    """描述单个中药成分的静态属性。

    Attributes:
        name: 成分名称。
        herb: 所属中药名。
        hepatotoxicity_score: SwissADME 肝毒性打分（越高越差）。
        synergy_baseline: 协同贡献基线（用于加权模型）。
        ob: 口服生物利用度（可选）。
        dl: 类药性指标（可选）。
    """

    name: str
    herb: str
    hepatotoxicity_score: float
    synergy_baseline: float
    ob: float = 0.0
    dl: float = 0.0


@dataclass
class CandidateSolution:
    """混合染色体：布尔选择 + 实数剂量。

    Attributes:
        selects: 每个成分是否被选中（True/False）。
        proportions: 与 selects 同步的混合剂量向量（非归一化）。
    """

    selects: List[bool]
    proportions: List[float]

    def normalized_proportions(self) -> List[float]:
        total = 0.0
        for sel, prop in zip(self.selects, self.proportions):
            if sel:
                total += prop
        if total <= 0:
            selected = [sel for sel in self.selects if sel]
            if not selected:
                return [0.0] * len(self.proportions)
            default = 1.0 / len(selected)
            return [default if sel else 0.0 for sel in self.selects]
        return [prop / total if sel else 0.0 for sel, prop in zip(self.selects, self.proportions)]

    def with_normalized(self) -> "CandidateSolution":
        return CandidateSolution(list(self.selects), self.normalized_proportions())


def ensure_selection(selects: List[bool]) -> List[bool]:
    """保证至少选中一个成分。"""
    if any(selects):
        return selects
    copy = list(selects)
    copy[np.random.randint(len(copy))] = True
    return copy


def vector_to_candidate(vector: Sequence[float], ingredients: Sequence[Ingredient]) -> CandidateSolution:
    """将 pymoo/pyswarms 的 [0,1] 向量映射为混合染色体。

    Args:
        vector: 兼具“是否选中”和“配比”信息的连续向量。
        ingredients: 组成向量所对应的成分列表，长度应为 vector 长度的一半。

    Returns:
        归一化后的 CandidateSolution。
    """
    half = len(vector) // 2
    selects = [float(value) > 0.5 for value in vector[:half]]
    selects = ensure_selection(selects)
    proportions = [float(value) for value in vector[half:]]
    proportions = [max(0.01, min(1.0, value)) for value in proportions]
    return CandidateSolution(selects, proportions).with_normalized()


def compute_loew_synergy(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """Loewe 指标的近似：配比加权协同基线。

    Args:
        candidate: 当前组合。
        ingredients: 成分数据集。

    Returns:
        加权协同得分。
    """
    normalized = candidate.with_normalized()
    synergy = 0.0
    for select, proportion, ingredient in zip(normalized.selects, normalized.proportions, ingredients):
        if select:
            synergy += ingredient.synergy_baseline * proportion
    return synergy


def compute_hepatotoxicity_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """计算组合的加权肝毒性得分。

    Args:
        candidate: 当前组合。
        ingredients: 成分数据集。

    Returns:
        配比加权后的 hepatotoxicity_score。
    """
    normalized = candidate.with_normalized()
    score = 0.0
    for select, proportion, ingredient in zip(normalized.selects, normalized.proportions, ingredients):
        if select:
            score += ingredient.hepatotoxicity_score * proportion
    return score


def diversity_penalty(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """若同一母药占比 ≥80%，返回 0.1 惩罚，否则 0。

    Args:
        candidate: 当前组合。
        ingredients: 成分信息（用于获取 herb 名称）。

    Returns:
        惩罚值（0.1 或 0）。
    """
    normalized = candidate.with_normalized()
    herb_proportions: dict[str, float] = {}
    for select, proportion, ingredient in zip(normalized.selects, normalized.proportions, ingredients):
        if not select:
            continue
        herb_proportions[ingredient.herb] = herb_proportions.get(ingredient.herb, 0.0) + proportion
    if not herb_proportions:
        return 0.0
    majority_ratio = max(herb_proportions.values())
    return 0.1 if majority_ratio >= 0.8 else 0.0


def evaluate_metrics(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float, float]:
    """返回 (协同, 肝毒性, 惩罚)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (协同得分, 肝毒性, 惩罚)
    """
    synergy = compute_loew_synergy(candidate, ingredients)
    toxicity = compute_hepatotoxicity_score(candidate, ingredients)
    penalty = diversity_penalty(candidate, ingredients)
    synergy_with_penalty = synergy * (1 - penalty)
    return synergy_with_penalty, toxicity, penalty


def single_objective_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient], toxicity_weight: float = 1.0) -> float:
    """单目标评价：协同 - toxicity_weight * 肝毒性。

    Args:
        candidate: 当前候选。
        ingredients: 成分列表。
        toxicity_weight: 毒性在适应度中的相对权重。

    Returns:
        复方综合得分（越高越优）。
    """
    synergy, toxicity, _ = evaluate_metrics(candidate, ingredients)
    return synergy - toxicity_weight * toxicity


def multi_objective_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float]:
    """返回 (协同, 肝毒性)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (协同得分, 肝毒性)。
    """
    synergy, toxicity, _ = evaluate_metrics(candidate, ingredients)
    return synergy, toxicity


class CombinationProblem(ElementwiseProblem):
    """pymoo 兼容问题定义，用于描述复方的选择+配比变量与协同/毒性目标。

    Args:
        ingredients: 候选成分列表。
        mode: \"single\" 表示单目标（协同 - 毒性），\"multi\" 表示多目标。
        toxicity_weight: 仅在单目标模式下生效，表示肝毒性惩罚权重。
    """

    def __init__(self, ingredients: Sequence[Ingredient], mode: str = "single", toxicity_weight: float = 1.0):
        n_var = len(ingredients) * 2
        n_obj = 1 if mode == "single" else 2
        super().__init__(n_var=n_var, n_obj=n_obj, xl=0.0, xu=1.0)
        self.ingredients = ingredients
        self.mode = mode
        self.toxicity_weight = toxicity_weight

    def _evaluate(self, x, out, *args, **kwargs):
        """在 pymoo 运行过程中评估单个向量。

        Args:
            x: 当前向量（shape: n_var）。
            out: 输出 dict，pymoo 要求填入 \"F\"。
        """
        candidate = vector_to_candidate(x, self.ingredients)
        synergy, toxicity = multi_objective_score(candidate, self.ingredients)
        if self.mode == "single":
            out["F"] = np.array([- (synergy - self.toxicity_weight * toxicity)])
        else:
            out["F"] = np.array([-synergy, toxicity])


class PymooSingleObjectiveGA:
    """使用 pymoo GA 处理单目标优化（协同 - 毒性）。"""

    def __init__(self, ingredients: Sequence[Ingredient], toxicity_weight: float = 1.0, generations: int = 25, population_size: int = 32):
        """配置 GA 所需组件。

        Args:
            ingredients: 成分列表。
            toxicity_weight: 毒性惩罚权重。
            generations: 演化代数。
            population_size: 种群规模。
        """
        self.problem = CombinationProblem(ingredients, mode="single", toxicity_weight=toxicity_weight)
        self.algorithm = PymooGAAlg(pop_size=population_size, eliminate_duplicates=True)
        self.generations = generations

    def run(self) -> CandidateSolution:
        """执行单目标 GA 并返回最优 candidate。"""
        result = minimize(self.problem, self.algorithm, termination=("n_gen", self.generations), verbose=False)
        return vector_to_candidate(result.X, self.problem.ingredients)


class PymooNSGAII:
    """基于 pymoo NSGA-II 的多目标优化器。"""

    def __init__(self, ingredients: Sequence[Ingredient], generations: int = 30, population_size: int = 40):
        """初始化 NSGA-II 算法。

        Args:
            ingredients: 成分列表。
            generations: 演化代数。
            population_size: 种群规模。
        """
        self.problem = CombinationProblem(ingredients, mode="multi")
        self.algorithm = NSGA2(pop_size=population_size)
        self.generations = generations

    def run(self) -> Tuple[List[CandidateSolution], np.ndarray]:
        """执行 NSGA-II 并返回 candidate 与目标值矩阵。

        Returns:
            solutions: 输入向量转化的 CandidateSolution 列表。
            np.ndarray: 每个 solution 的目标值（协同取负 + 肝毒性）。
        """
        result = minimize(self.problem, self.algorithm, termination=("n_gen", self.generations), verbose=False)
        solutions = [vector_to_candidate(x, self.problem.ingredients) for x in result.X]
        return solutions, result.F


class PySwarmsPSO:
    """借助 pyswarms 实现的单目标粒子群优化。"""

    def __init__(self, ingredients: Sequence[Ingredient], iterations: int = 40, swarm_size: int = 24):
        """配置 PSO 超参数。

        Args:
            ingredients: 成分列表。
            iterations: 迭代次数。
            swarm_size: 粒子数量。
        """
        self.ingredients = ingredients
        self.iterations = iterations
        self.swarm_size = swarm_size
        self.dimensions = len(ingredients) * 2
        self.options = {"c1": 1.2, "c2": 1.2, "w": 0.6}

    def _fitness(self, swarm: np.ndarray) -> np.ndarray:
        """pyswarms 回调：计算每个粒子的目标值（越小越优）。"""
        scores = []
        for particle in swarm:
            candidate = vector_to_candidate(particle, self.ingredients)
            score = single_objective_score(candidate, self.ingredients)
            scores.append(-score)
        return np.array(scores)

    def run(self) -> CandidateSolution:
        """执行 PSO 并返回最优 candidate。"""
        optimizer = GlobalBestPSO(
            n_particles=self.swarm_size,
            dimensions=self.dimensions,
            options=self.options,
            bounds=(np.zeros(self.dimensions), np.ones(self.dimensions)),
        )
        _, best_pos = optimizer.optimize(self._fitness, iters=self.iterations, verbose=False)
        return vector_to_candidate(best_pos, self.ingredients)
