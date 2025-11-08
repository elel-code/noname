from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass
from typing import Callable, List, Sequence, Tuple


@dataclass(frozen=True)
class Ingredient:
    """存储单个中药成分的静态属性，供算法引用。"""
    name: str
    herb: str
    hepatotoxicity_score: float
    synergy_baseline: float
    ob: float = 0.0
    dl: float = 0.0


@dataclass
class CandidateSolution:
    selects: List[bool]
    proportions: List[float]

    """表示成分选择（二进制）与配比（实数）的混合染色体。"""

    def normalized_proportions(self) -> List[float]:
        """归一化选中成分的配比，避免总量为零。"""
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
        """返回配比已归一化的解拷贝。"""
        return CandidateSolution(list(self.selects), self.normalized_proportions())


def ensure_selection(selects: List[bool]) -> List[bool]:
    """确保至少有一个成分被选中，避免空组合。

    Args:
        selects: 原始选择布尔矩阵。

    Returns:
        至少启用一个成分的选择向量。
    """
    if any(selects):
        return selects
    copy = list(selects)
    copy[random.randrange(len(copy))] = True
    return copy


def make_random_candidate(num_ingredients: int) -> CandidateSolution:
    """生成随机的选择+配比候选并归一化。

    Args:
        num_ingredients: 可供选择的成分总数。

    Returns:
        包含随机布尔选择与实数配比的候选组合。
    """
    selects = [random.random() < 0.5 for _ in range(num_ingredients)]
    selects = ensure_selection(selects)
    proportions = [random.random() for _ in range(num_ingredients)]
    solution = CandidateSolution(selects, proportions)
    return solution.with_normalized()


def compute_loew_synergy(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """作为 Loewe Additivity 的近似：按配比加权叠加各成分协同基线贡献。

    Args:
        candidate: 当前解的选择与配比。
        ingredients: 所有可选成分的属性信息。

    Returns:
        加权后的协同得分。
    """
    normalized = candidate.with_normalized()
    synergy = 0.0
    for select, proportion, ingredient in zip(normalized.selects, normalized.proportions, ingredients):
        if select:
            synergy += ingredient.synergy_baseline * proportion
    return synergy


def compute_hepatotoxicity_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """计算组合的加权肝毒性得分，来自 SwissADME 输出。

    Args:
        candidate: 当前解的混合编码。
        ingredients: 成分的 hepatotoxicity_score。

    Returns:
        组合肝毒性加权值。
    """
    normalized = candidate.with_normalized()
    score = 0.0
    for select, proportion, ingredient in zip(normalized.selects, normalized.proportions, ingredients):
        if select:
            score += ingredient.hepatotoxicity_score * proportion
    return score


def diversity_penalty(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> float:
    """若某一味药占比 ≥80%，则附加罚分，鼓励配伍多样性。

    Args:
        candidate: 当前组合的选择情况。
        ingredients: 成分清单含所属药名。

    Returns:
        0.1 或 0.0 表示是否触发惩罚。
    """
    normalized = candidate.with_normalized()
    selected_herbs = [ingredient.herb for select, ingredient in zip(normalized.selects, ingredients) if select]
    if not selected_herbs:
        return 0.0
    counts = Counter(selected_herbs)
    majority_ratio = max(counts.values()) / len(selected_herbs)
    return 0.1 if majority_ratio >= 0.8 else 0.0


def evaluate_metrics(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float, float]:
    """返回协同、毒性与惩罚三项，用于适应度整合。

    Args:
        candidate: 混合编码组合。
        ingredients: 可选成分列表。

    Returns:
        (协同得分, 肝毒性, 惩罚)。
    """
    synergy = compute_loew_synergy(candidate, ingredients)
    toxicity = compute_hepatotoxicity_score(candidate, ingredients)
    penalty = diversity_penalty(candidate, ingredients)
    synergy_with_penalty = synergy * (1 - penalty)
    return synergy_with_penalty, toxicity, penalty


def single_objective_score(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient], toxicity_weight: float = 1.0
) -> float:
    """单目标适应度：协同减去毒性，toxicity_weight 可调。

    Args:
        candidate: 混合编码解。
        ingredients: 可选成分列表。
        toxicity_weight: 肝毒性在适应度中的权重。

    Returns:
        综合得分，越大越优。
    """
    synergy, toxicity, _ = evaluate_metrics(candidate, ingredients)
    return synergy - toxicity_weight * toxicity


def multi_objective_score(candidate: CandidateSolution, ingredients: Sequence[Ingredient]) -> Tuple[float, float]:
    """返回 (协同, 肝毒性) 供 NSGA-II 进行非支配排序。

    Args:
        candidate: 当前组合解。
        ingredients: 成分列表。

    Returns:
        (协同得分, 肝毒性)。
    """
    synergy, toxicity, _ = evaluate_metrics(candidate, ingredients)
    return synergy, toxicity


class GeneticAlgorithm:
    """实现二进制选择+实数配比的遗传算法框架。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        fitness_function: Callable[[CandidateSolution], float],
        mutation_rate: float = 0.1,
    ):
        """初始化遗传算法。

        Args:
            ingredients: 待筛选的成分列表。
            fitness_function: 评估 candidate 的打分函数。
            mutation_rate: 每个基因突变的概率。
        """
        self.ingredients = ingredients
        self.fitness_function = fitness_function
        self.mutation_rate = mutation_rate

    def run(self, generations: int = 25, population_size: int = 32) -> CandidateSolution:
        """按迭代次数演化群体，返回适应度最高的候选。

        Args:
            generations: 演化代数。
            population_size: 单代群体大小。

        Returns:
            最优候选（按给定 fitness 函数）。
        """
        population = [make_random_candidate(len(self.ingredients)) for _ in range(population_size)]
        for _ in range(generations):
            population = sorted(population, key=self.fitness_function, reverse=True)
            survivors = population[: population_size // 2]
            children = []
            while len(children) < population_size - len(survivors):
                parent_a, parent_b = random.sample(survivors, 2)
                child = self.crossover(parent_a, parent_b)
                child = self.mutate(child)
                children.append(child)
            population = survivors + children
        best = max(population, key=self.fitness_function)
        return best.with_normalized()

    def crossover(self, parent_a: CandidateSolution, parent_b: CandidateSolution) -> CandidateSolution:
        """将两个父代通过单点交叉生成子代。

        Args:
            parent_a: 父代 A。
            parent_b: 父代 B。

        Returns:
            交叉后的子代 candidate。
        """
        split = random.randrange(1, len(parent_a.selects))
        selects = parent_a.selects[:split] + parent_b.selects[split:]
        proportions = parent_b.proportions[:split] + parent_a.proportions[split:]
        return CandidateSolution(ensure_selection(selects), proportions).with_normalized()

    def mutate(self, candidate: CandidateSolution) -> CandidateSolution:
        """对选择位随机翻转，并微调配比以引入多样性。

        Args:
            candidate: 待变异个体。

        Returns:
            变异后的个体。
        """
        selects = candidate.selects[:]
        proportions = candidate.proportions[:]
        for i in range(len(selects)):
            if random.random() < self.mutation_rate:
                selects[i] = not selects[i]
        selects = ensure_selection(selects)
        for i in range(len(proportions)):
            if random.random() < self.mutation_rate:
                delta = random.uniform(-0.25, 0.25)
                proportions[i] = max(0.01, min(1.0, proportions[i] + delta))
        return CandidateSolution(selects, proportions).with_normalized()


class ParticleSwarmOptimization:
    """粒子群算法，采用同一混合向量表示选择和配比。"""

    def __init__(self, ingredients: Sequence[Ingredient], fitness_function: Callable[[CandidateSolution], float]):
        """初始化粒子群。

        Args:
            ingredients: 成分列表。
            fitness_function: 评价函数。
        """
        self.ingredients = ingredients
        self.dimension = len(ingredients) * 2
        self.fitness_function = fitness_function

    def run(self, iterations: int = 40, swarm_size: int = 24) -> CandidateSolution:
        """执行若干代迭代，并返回最佳粒子对应的候选解。

        Args:
            iterations: 总迭代次数。
            swarm_size: 粒子数量。

        Returns:
            最佳粒子的 candidate。
        """
        positions = [self._random_vector() for _ in range(swarm_size)]
        velocities = [[random.uniform(-0.2, 0.2) for _ in range(self.dimension)] for _ in range(swarm_size)]
        personal_best_positions = list(positions)
        personal_best_scores = [self._evaluate_vector(vec) for vec in positions]
        global_best_position = personal_best_positions[int(max(range(swarm_size), key=lambda i: personal_best_scores[i]))]
        global_best_score = max(personal_best_scores)
        for _ in range(iterations):
            for idx in range(swarm_size):
                velocities[idx] = self._update_velocity(velocities[idx], positions[idx], personal_best_positions[idx], global_best_position)
                positions[idx] = self._apply_velocity(positions[idx], velocities[idx])
                score = self._evaluate_vector(positions[idx])
                if score > personal_best_scores[idx]:
                    personal_best_scores[idx] = score
                    personal_best_positions[idx] = list(positions[idx])
                if score > global_best_score:
                    global_best_score = score
                    global_best_position = list(positions[idx])
        return self._vector_to_candidate(global_best_position).with_normalized()

    def _evaluate_vector(self, vector: List[float]) -> float:
        """将向量转换为 candidate 并打分。

        Args:
            vector: 当前粒子向量。

        Returns:
            适应度得分。
        """
        candidate = self._vector_to_candidate(vector)
        return self.fitness_function(candidate)

    def _random_vector(self) -> List[float]:
        """生成 [0, 1] 区间的初始粒子。

        Returns:
            随机初始化的向量。
        """
        return [random.random() for _ in range(self.dimension)]

    def _vector_to_candidate(self, vector: List[float]) -> CandidateSolution:
        """将粒子向量转为混合选择+配比候选。

        Args:
            vector: 粒子状态向量。

        Returns:
            CandidateSolution 实例。
        """
        half = len(vector) // 2
        selects = [value > 0.5 for value in vector[:half]]
        selects = ensure_selection(selects)
        proportions = [max(0.01, min(1.0, value)) for value in vector[half:]]
        return CandidateSolution(selects, proportions).with_normalized()

    def _apply_velocity(self, position: List[float], velocity: List[float]) -> List[float]:
        """将速度应用到当前向量并限制到合法区间。

        Args:
            position: 当前向量。
            velocity: 速度矢量。

        Returns:
            更新后的向量（0-1 限制）。
        """
        return [min(1.0, max(0.0, p + v)) for p, v in zip(position, velocity)]

    def _update_velocity(self, velocity: List[float], position: List[float], personal_best: List[float], global_best: List[float]) -> List[float]:
        """按照惯性 + 认知 + 社会成分更新速度向量。

        Args:
            velocity: 当前速度。
            position: 当前向量。
            personal_best: 个人最优状态。
            global_best: 群体最优状态。

        Returns:
            更新后的速度。
        """
        inertia = 0.6
        cognitive = 1.2
        social = 1.2
        new_velocity = []
        for v, x, p, g in zip(velocity, position, personal_best, global_best):
            r1, r2 = random.random(), random.random()
            cognitive_term = cognitive * r1 * (p - x)
            social_term = social * r2 * (g - x)
            new_velocity.append(inertia * v + cognitive_term + social_term)
        return new_velocity


def dominates(score_a: Tuple[float, float], score_b: Tuple[float, float]) -> bool:
    synergy_a, tox_a = score_a
    synergy_b, tox_b = score_b
    return (synergy_a >= synergy_b and tox_a <= tox_b) and (synergy_a > synergy_b or tox_a < tox_b)


def fast_non_dominated_sort(scores: List[Tuple[float, float]]) -> List[List[int]]:
    size = len(scores)
    """快速非支配排序，返回按层级列表中存放索引。"""
    dominates_list: List[List[int]] = [[] for _ in range(size)]
    dominated_count = [0] * size
    fronts: List[List[int]] = [[]]
    for p in range(size):
        for q in range(size):
            if dominates(scores[p], scores[q]):
                dominates_list[p].append(q)
            elif dominates(scores[q], scores[p]):
                dominated_count[p] += 1
        if dominated_count[p] == 0:
            fronts[0].append(p)
    current = 0
    while fronts[current]:
        next_front: List[int] = []
        for index in fronts[current]:
            for dominated in dominates_list[index]:
                dominated_count[dominated] -= 1
                if dominated_count[dominated] == 0:
                    next_front.append(dominated)
        current += 1
        fronts.append(next_front)
    return [front for front in fronts if front]


def crowding_distance(scores: List[Tuple[float, float]], front: List[int]) -> dict:
    """计算拥挤度，用于优先保留 Pareto 前沿的边界和多样性。"""
    distance = {idx: 0.0 for idx in front}
    if not front:
        return distance
    for objective_index in range(2):
        reverse = objective_index == 0
        sorted_front = sorted(front, key=lambda idx: scores[idx][objective_index], reverse=reverse)
        min_value = scores[sorted_front[-1]][objective_index]
        max_value = scores[sorted_front[0]][objective_index]
        distance[sorted_front[0]] = distance[sorted_front[-1]] = float("inf")
        if max_value == min_value:
            continue
        for i in range(1, len(sorted_front) - 1):
            prev_score = scores[sorted_front[i - 1]][objective_index]
            next_score = scores[sorted_front[i + 1]][objective_index]
            distance[sorted_front[i]] += (prev_score - next_score) / (max_value - min_value)
    return distance


class NSGAII:
    """NSGA-II 多目标优化器，基于混合染色体并利用遗传突变。"""

    def __init__(self, ingredients: Sequence[Ingredient], mutation_rate: float = 0.1):
        self.ingredients = ingredients
        self.mutator = GeneticAlgorithm(
            ingredients, lambda sol: single_objective_score(sol, ingredients), mutation_rate
        )

    def run(self, generations: int = 30, population_size: int = 40) -> List[CandidateSolution]:
        """执行 NSGA-II 主循环并返回最终候选集。

        Args:
            generations: 演化代数。
            population_size: 每代候选数量。

        Returns:
            末代候选列表。
        """
        population = [make_random_candidate(len(self.ingredients)) for _ in range(population_size)]
        for _ in range(generations):
            offspring = [self.mutator.mutate(random.choice(population)) for _ in range(population_size)]
            population = self._select_next_generation(population + offspring, population_size)
        return population

    def _select_next_generation(self, combined: List[CandidateSolution], size: int) -> List[CandidateSolution]:
        """按层级与拥挤度筛选下一代。

        Args:
            combined: 当前代 + 子代构成的整体列表。
            size: 下一代规模。

        Returns:
            按 NSGA-II 策略筛选后的候选列表。
        """
        evaluations = [(candidate, multi_objective_score(candidate, self.ingredients)) for candidate in combined]
        scores = [score for _, score in evaluations]
        fronts = fast_non_dominated_sort(scores)
        next_population: List[CandidateSolution] = []
        for front in fronts:
            if len(next_population) + len(front) > size:
                distances = crowding_distance(scores, front)
                sorted_front = sorted(front, key=lambda idx: distances[idx], reverse=True)
                needed = size - len(next_population)
                next_population.extend([evaluations[idx][0] for idx in sorted_front[:needed]])
                break
            next_population.extend([evaluations[idx][0] for idx in front])
        return next_population

    def pareto_front(self, population: List[CandidateSolution]) -> List[Tuple[CandidateSolution, Tuple[float, float]]]:
        """提取非支配解集合，便于输出协同/毒性 Trade-off。

        Args:
            population: 待分类的候选群体。

        Returns:
            非支配解列表（candidate + 对应 score）。
        """
        evaluated = [(candidate, multi_objective_score(candidate, self.ingredients)) for candidate in population]
        non_dominated: List[Tuple[CandidateSolution, Tuple[float, float]]] = []
        for candidate, score in evaluated:
            if not any(dominates(other_score, score) for _, other_score in evaluated if other_score != score):
                non_dominated.append((candidate, score))
        return non_dominated
