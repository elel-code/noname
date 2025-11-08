from __future__ import annotations

import math
from array import array
import copy
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA as PymooGAAlg
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pyswarms.single.global_best import GlobalBestPSO

CONFIG_PATH = Path(__file__).resolve().parent / "config" / "algorithms.json"
_DEFAULT_ALGO_CONFIG = {
    "ga": {"generations": 25, "population_size": 32, "toxicity_weight": 1.0},
    "nsga2": {"generations": 30, "population_size": 40},
    "pso": {
        "iterations": 40,
        "swarm_size": 24,
        "options": {"c1": 1.2, "c2": 1.2, "w": 0.6},
    },
    "aci": {
        "lambda_weight": 0.3,
        "complement_variance_max": 0.25,
        "coverage_floor": 0.0,
    },
}
_ALGO_CONFIG_CACHE: dict | None = None


def _deep_update(base: dict, overrides: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_update(result[key], value)
        else:
            result[key] = value
    return result


def _load_algorithm_config() -> dict:
    global _ALGO_CONFIG_CACHE
    if _ALGO_CONFIG_CACHE is not None:
        return _ALGO_CONFIG_CACHE
    config = copy.deepcopy(_DEFAULT_ALGO_CONFIG)
    try:
        with CONFIG_PATH.open("r", encoding="utf-8") as fp:
            user_config = json.load(fp)
            config = _deep_update(config, user_config)
    except FileNotFoundError:
        pass
    _ALGO_CONFIG_CACHE = config
    return config


def _get_aci_settings() -> dict:
    """便捷访问 ACI 相关配置。"""

    return _load_algorithm_config().get("aci", _DEFAULT_ALGO_CONFIG["aci"])


@dataclass(frozen=True)
class Ingredient:
    """描述单个中药成分的静态属性。

    Attributes:
        name: 成分名称。
        herb: 所属中药名。
        hepatotoxicity_score: 外部模型/实验给出的肝毒性打分（越高越差）。
        synergy_baseline: 协同贡献基线（用于加权模型）。
        ob: Oral Bioavailability，若缺失可置 None。
        dl: Drug Likeness，若缺失可置 None。
        mw: 分子量 (MW)。
        alogp: ALogP（脂溶性）。
        h_don: 氢键供体数量。
        h_acc: 氢键受体数量。
        caco2: Caco-2 渗透性（可用 logPapp）。
        bbb: 血脑屏障穿透（logBB）。
        fasa: FASA- 极性表面积比。
        hl: 半衰期（小时）。
    """

    name: str
    herb: str
    hepatotoxicity_score: float
    synergy_baseline: float
    ob: float | None = None
    dl: float | None = None
    mw: float | None = None
    alogp: float | None = None
    h_don: float | None = None
    h_acc: float | None = None
    caco2: float | None = None
    bbb: float | None = None
    fasa: float | None = None
    hl: float | None = None


@dataclass(frozen=True)
class FeatureWindow:
    """用于定义 ADME 特征的“宜居区”。

    note 会描述该特征的单位或取值口径，方便事后定位。
    """

    lower: float
    upper: float
    softness: float
    note: str = ""


ADME_FEATURE_ATTRS: dict[str, str] = {
    "mw": "mw",
    "alogp": "alogp",
    "hdon": "h_don",
    "hacc": "h_acc",
    "ob": "ob",
    "caco2": "caco2",
    "bbb": "bbb",
    "dl": "dl",
    "fasa": "fasa",
    "hl": "hl",
}

FEATURE_WINDOWS: dict[str, FeatureWindow] = {
    "mw": FeatureWindow(180.0, 500.0, 120.0, note="Dalton"),
    "alogp": FeatureWindow(-0.5, 5.0, 1.5, note="ALogP"),
    "hdon": FeatureWindow(0.0, 5.0, 2.0, note="count"),
    "hacc": FeatureWindow(0.0, 10.0, 3.0, note="count"),
    "ob": FeatureWindow(30.0, 70.0, 15.0, note="OB%"),
    "caco2": FeatureWindow(0.4, 1.5, 0.3, note="log10(Papp, cm/s)"),
    "bbb": FeatureWindow(-1.0, 1.0, 0.5, note="logBB or mapped分类"),
    "dl": FeatureWindow(0.18, 0.5, 0.1, note="DL"),
    "fasa": FeatureWindow(0.2, 0.6, 0.15, note="FASA-"),
    "hl": FeatureWindow(3.0, 12.0, 3.0, note="hours"),
}

ACI_LAMBDA = 0.3
_SAFE_EPS = 1e-6
_TOTAL_FEATURES = len(ADME_FEATURE_ATTRS)


BITS_PER_BYTE = 8


def _pack_bools(values: Sequence[bool]) -> bytearray:
    data = bytearray(math.ceil(len(values) / BITS_PER_BYTE))
    for idx, value in enumerate(values):
        if value:
            data[idx // BITS_PER_BYTE] |= 1 << (idx % BITS_PER_BYTE)
    return data


def _bit_is_set(bits: bytearray, idx: int) -> bool:
    return bool(bits[idx // BITS_PER_BYTE] & (1 << (idx % BITS_PER_BYTE)))


@dataclass
class CandidateSolution:
    """混合染色体：布尔选择 + 实数配比（使用位存储减少内存占用）。"""

    select_bits: bytearray
    proportions: array
    length: int

    @classmethod
    def from_components(cls, selects: Sequence[bool], proportions: Sequence[float]) -> "CandidateSolution":
        return cls(_pack_bools(selects), array("f", proportions), len(selects))

    def iter_selects(self) -> Iterable[bool]:
        for idx in range(self.length):
            yield _bit_is_set(self.select_bits, idx)

    def iter_selected_indices(self) -> Iterable[int]:
        for idx, sel in enumerate(self.iter_selects()):
            if sel:
                yield idx

    def normalized_proportions(self) -> array:
        selects = list(self.iter_selects())
        total = sum(prop for sel, prop in zip(selects, self.proportions) if sel)
        normalized = array("f", [0.0] * self.length)
        if total <= 0:
            selected_indices = [idx for idx, sel in enumerate(selects) if sel]
            if not selected_indices:
                return normalized
            default = 1.0 / len(selected_indices)
            for idx in selected_indices:
                normalized[idx] = default
            return normalized
        for idx, (sel, prop) in enumerate(zip(selects, self.proportions)):
            if sel:
                normalized[idx] = prop / total
        return normalized

    def with_normalized(self) -> "CandidateSolution":
        return CandidateSolution(bytearray(self.select_bits), self.normalized_proportions(), self.length)


def ensure_selection(selects: List[bool], ingredients: Sequence[Ingredient]) -> List[bool]:
    """保证至少选中一个成分；若全为 False，则启用协同基线最高的成分以保持可复现。"""
    if any(selects):
        return selects
    copy = list(selects)
    best_index = max(range(len(ingredients)), key=lambda idx: ingredients[idx].synergy_baseline)
    copy[best_index] = True
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
    selects = ensure_selection(selects, ingredients)
    proportions = [float(value) for value in vector[half:]]
    proportions = [max(0.01, min(1.0, value)) for value in proportions]
    return CandidateSolution.from_components(selects, proportions)


def _window_score(value: float, window: FeatureWindow) -> float:
    """根据 FeatureWindow 计算单一特征的得分。"""

    if window.lower <= value <= window.upper:
        return 1.0
    delta = window.lower - value if value < window.lower else value - window.upper
    return math.exp(-((delta / max(window.softness, _SAFE_EPS)) ** 2))


def _prepare_feature_value(feature: str, value: float | None) -> float | None:
    """根据特征类型调整原始值（如 BBB 的 0/1 映射）。"""

    if value is None:
        return None
    if feature == "bbb":
        if isinstance(value, bool):
            return 0.5 if value else -0.5
        if isinstance(value, (int, float)):
            if math.isclose(value, 0.0, rel_tol=1e-3, abs_tol=1e-3):
                return -0.5
            if math.isclose(value, 1.0, rel_tol=1e-3, abs_tol=1e-3):
                return 0.5
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in {"y", "yes", "pass", "true", "+"}:
                return 0.5
            if lowered in {"n", "no", "fail", "false", "-"}:
                return -0.5
    return float(value)


def _normalized_feature_value(feature: str, value: float | None) -> float | None:
    """将原始特征值映射到 0-1 区间；缺失值返回 None。"""

    prepared = _prepare_feature_value(feature, value)
    if prepared is None:
        return None
    window = FEATURE_WINDOWS.get(feature)
    if window is None:
        return None
    return _window_score(prepared, window)


def compute_desirability_score(ingredient: Ingredient) -> float:
    """计算单体可用性分（几何平均）。"""

    scores: list[float] = []
    for feature, attr in ADME_FEATURE_ATTRS.items():
        value = getattr(ingredient, attr)
        normalized = _normalized_feature_value(feature, value)
        if normalized is not None:
            scores.append(max(normalized, _SAFE_EPS))
    if not scores:
        return 0.0
    product = 1.0
    for score in scores:
        product *= score
    geo_mean = product ** (1 / len(scores))
    coverage = len(scores) / max(_TOTAL_FEATURES, 1)
    coverage_floor = max(0.0, min(1.0, _get_aci_settings().get("coverage_floor", 0.0)))
    return geo_mean * max(coverage, coverage_floor)


def compute_complementarity_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """根据特征多样性计算组合互补性。"""

    feature_diversities: list[float] = []
    variance_cap = max(
        _SAFE_EPS,
        float(_get_aci_settings().get("complement_variance_max", 0.25)),
    )
    for feature, attr in ADME_FEATURE_ATTRS.items():
        values: list[float] = []
        weights: list[float] = []
        for select, proportion, ingredient in zip(candidate.iter_selects(), candidate.proportions, ingredients):
            if not select:
                continue
            normalized = _normalized_feature_value(feature, getattr(ingredient, attr))
            if normalized is None:
                continue
            values.append(normalized)
            weights.append(proportion)
        if len(values) <= 1:
            continue
        weights_arr = np.array(weights, dtype=float)
        total = weights_arr.sum()
        if total <= 0:
            continue
        weights_arr /= total
        values_arr = np.array(values, dtype=float)
        mean = float(np.dot(weights_arr, values_arr))
        variance = float(np.dot(weights_arr, (values_arr - mean) ** 2))
        feature_diversities.append(min(1.0, variance / variance_cap))
    if not feature_diversities:
        return 0.0
    return float(sum(feature_diversities) / len(feature_diversities))


def compute_aci_score(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    lambda_weight: float | None = None,
) -> float:
    """计算基于 ADME 特征的 ACI 指标。"""

    desirabilities: list[tuple[float, float]] = []
    for select, proportion, ingredient in zip(candidate.iter_selects(), candidate.proportions, ingredients):
        if not select:
            continue
        desirability = compute_desirability_score(ingredient)
        desirabilities.append((proportion, desirability))
    if not desirabilities:
        return 0.0
    total_prop = sum(weight for weight, _ in desirabilities)
    if total_prop <= 0:
        return 0.0
    weighted_mean = sum(weight * desirability for weight, desirability in desirabilities) / total_prop
    complementarity = compute_complementarity_score(candidate, ingredients)
    if lambda_weight is None:
        lambda_weight = _get_aci_settings().get("lambda_weight", ACI_LAMBDA)
    lambda_weight = max(0.0, min(1.0, lambda_weight))
    return (1 - lambda_weight) * weighted_mean + lambda_weight * complementarity


def compute_hepatotoxicity_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """计算组合的加权肝毒性得分。

    Args:
        candidate: 已归一化组合。
        ingredients: 成分数据集。

    Returns:
        配比加权后的 hepatotoxicity_score。
    """
    score = 0.0
    for select, proportion, ingredient in zip(candidate.iter_selects(), candidate.proportions, ingredients):
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
    herb_proportions: dict[str, float] = {}
    for select, proportion, ingredient in zip(candidate.iter_selects(), candidate.proportions, ingredients):
        if not select:
            continue
        herb_proportions[ingredient.herb] = herb_proportions.get(ingredient.herb, 0.0) + proportion
    if not herb_proportions:
        return 0.0
    majority_ratio = max(herb_proportions.values())
    return 0.1 if majority_ratio >= 0.8 else 0.0


def evaluate_metrics(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float, float, float]:
    """返回 (ACI(含惩罚), 肝毒性, 惩罚, 纯 ACI)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (ACI, 肝毒性, 惩罚, ACI_raw)
    """
    normalized_candidate = candidate.with_normalized()
    aci_raw = compute_aci_score(normalized_candidate, ingredients)
    toxicity = compute_hepatotoxicity_score(normalized_candidate, ingredients)
    penalty = diversity_penalty(normalized_candidate, ingredients)
    aci_with_penalty = aci_raw * (1 - penalty)
    return aci_with_penalty, toxicity, penalty, aci_raw


def single_objective_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient], toxicity_weight: float = 1.0) -> float:
    """单目标评价：ACI - toxicity_weight * 肝毒性。

    Args:
        candidate: 当前候选。
        ingredients: 成分列表。
        toxicity_weight: 毒性在适应度中的相对权重。

    Returns:
        复方综合得分（越高越优）。
    """
    aci, toxicity, _, _ = evaluate_metrics(candidate, ingredients)
    return aci - toxicity_weight * toxicity


def multi_objective_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float]:
    """返回 (ACI, 肝毒性)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (ACI, 肝毒性)。
    """
    aci, toxicity, _, _ = evaluate_metrics(candidate, ingredients)
    return aci, toxicity


class CombinationProblem(ElementwiseProblem):
    """pymoo 兼容问题定义，用于描述复方的选择+配比变量与 ACI/毒性目标。

    Args:
        ingredients: 候选成分列表。
        mode: \"single\" 表示单目标（ACI - 毒性），\"multi\" 表示多目标。
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
        aci, toxicity = multi_objective_score(candidate, self.ingredients)
        if self.mode == "single":
            out["F"] = np.array([-(aci - self.toxicity_weight * toxicity)])
        else:
            out["F"] = np.array([-aci, toxicity])


class PymooSingleObjectiveGA:
    """使用 pymoo GA 处理单目标优化（ACI - 毒性）。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        toxicity_weight: float | None = None,
        generations: int | None = None,
        population_size: int | None = None,
    ):
        """配置 GA 所需组件；默认参数来自 config/algorithms.json。"""
        cfg = _load_algorithm_config()["ga"]
        self.toxicity_weight = toxicity_weight if toxicity_weight is not None else cfg["toxicity_weight"]
        self.generations = generations if generations is not None else cfg["generations"]
        pop = population_size if population_size is not None else cfg["population_size"]
        self.problem = CombinationProblem(ingredients, mode="single", toxicity_weight=self.toxicity_weight)
        self.algorithm = PymooGAAlg(pop_size=pop, eliminate_duplicates=True)

    def run(self) -> CandidateSolution:
        """执行单目标 GA 并返回最优 candidate。"""
        result = minimize(self.problem, self.algorithm, termination=("n_gen", self.generations), verbose=False)
        return vector_to_candidate(result.X, self.problem.ingredients)


class PymooNSGAII:
    """基于 pymoo NSGA-II 的多目标优化器。"""

    def __init__(self, ingredients: Sequence[Ingredient], generations: int | None = None, population_size: int | None = None):
        """初始化 NSGA-II；默认参数读取配置文件。"""
        cfg = _load_algorithm_config()["nsga2"]
        self.generations = generations if generations is not None else cfg["generations"]
        pop = population_size if population_size is not None else cfg["population_size"]
        self.problem = CombinationProblem(ingredients, mode="multi")
        self.algorithm = NSGA2(pop_size=pop)

    def run(self) -> Tuple[List[CandidateSolution], np.ndarray]:
        """执行 NSGA-II 并返回 candidate 与目标值矩阵。

        Returns:
            solutions: 输入向量转化的 CandidateSolution 列表。
            np.ndarray: 每个 solution 的目标值（ACI 取负 + 肝毒性）。
        """
        result = minimize(self.problem, self.algorithm, termination=("n_gen", self.generations), verbose=False)
        solutions = [vector_to_candidate(x, self.problem.ingredients) for x in result.X]
        return solutions, result.F


class PySwarmsPSO:
    """借助 pyswarms 实现的单目标粒子群优化。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        iterations: int | None = None,
        swarm_size: int | None = None,
        options: dict | None = None,
    ):
        """配置 PSO 超参数，默认来自配置文件。"""
        cfg = _load_algorithm_config()["pso"]
        self.ingredients = ingredients
        self.iterations = iterations if iterations is not None else cfg["iterations"]
        self.swarm_size = swarm_size if swarm_size is not None else cfg["swarm_size"]
        self.options = options if options is not None else cfg["options"]
        self.dimensions = len(ingredients) * 2

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


def select_candidate_by_weighted_sum(
    pareto_candidates: Sequence[CandidateSolution],
    ingredients: Sequence[Ingredient],
    toxicity_weight: float,
) -> Tuple[CandidateSolution, dict]:
    """从 Pareto 集合中按权重（ACI-毒性）挑选最优解。

    Returns:
        (最佳 candidate, 指标字典)。
    """

    if not pareto_candidates:
        raise ValueError("Pareto 集合为空，无法做权重挑选。")
    best = None
    best_score = -math.inf
    for candidate in pareto_candidates:
        score = single_objective_score(candidate, ingredients, toxicity_weight)
        if score > best_score:
            best = candidate
            best_score = score
    assert best is not None
    aci_pen, tox, penalty, aci_raw = evaluate_metrics(best, ingredients)
    return best, {
        "score": best_score,
        "aci_penalized": aci_pen,
        "aci_raw": aci_raw,
        "toxicity": tox,
        "penalty": penalty,
    }
