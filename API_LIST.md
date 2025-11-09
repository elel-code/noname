# API 文档（更新）

本文档描述 `algorithms.py` 与 `main.py` 的主要类型与函数，并补充新增的“成分读取接口”。

---

## algorithms.py

### 数据结构

- `Ingredient`（dataclass）
  - 字段：
    - 标识与来源：`name: str`、`smiles: str | None`
    - 安全性：`hepatotoxicity_score: float | None`
    - ADME：`ob, dl, mw, alogp, h_don, h_acc, caco2, bbb, fasa, hl`（均为 `float | None`）
  - 行为：
    - 若未提供 `hepatotoxicity_score` 且存在 `smiles`，在 `__post_init__` 中自动估算（优先 DILIPred，不可用则使用启发式）。

- `CandidateSolution`
  - 结构：`select_bits: bytearray`（位存储选择）、`proportions: array('f')`（配比）、`length: int`
  - 方法：`from_components(selects, proportions)`、`iter_selects()`、`iter_selected_indices()`、`normalized_proportions()`、`with_normalized()`

### 数据读取接口（新增）

- `load_ingredients_csv(path: str | Path, encoding='utf-8') -> list[Ingredient]`
  - 支持常见同义列名（不区分大小写），至少需要 `name`，且 `smiles` 与 `hepatotoxicity_score` 至少其一存在。
  - `bbb` 既可为 0/1 或 true/false，也可为 logBB 数值。

- `load_ingredients_json(path: str | Path, encoding='utf-8') -> list[Ingredient]`
  - 期望 JSON 根为数组；键名支持与 CSV 相同的同义映射。

- `load_ingredients(path: str | Path, fmt: str | None = None, encoding='utf-8') -> list[Ingredient]`
  - 通用入口；`fmt` 为 `csv/json/None(auto)`；当 `fmt=None` 时按扩展名推断。

示例：
```python
from algorithms import load_ingredients
ings = load_ingredients("resources/templates/ingredients.csv", fmt="csv")
```

### 评估与指标

- `compute_desirability_score(ingredient) -> float`
- `compute_complementarity_score(candidate, ingredients) -> float`
- `compute_aci_score(candidate, ingredients, lambda_weight=None) -> float`
- `compute_hepatotoxicity_score(candidate, ingredients) -> float`
- `evaluate_metrics(candidate, ingredients) -> (aci_with_penalty, toxicity, penalty, aci_raw)`
- `diversity_penalty(candidate, ingredients) -> float`

### 目标函数

- `single_objective_score(candidate, ingredients, toxicity_weight=1.0) -> float`
- `multi_objective_score(candidate, ingredients) -> (aci, toxicity)`

### 算法封装

- `CombinationProblem`（pymoo 兼容问题定义）
- `PymooSingleObjectiveGA(ingredients, ...)` → `.run() -> CandidateSolution`
- `PymooNSGAII(ingredients, ...)` → `.run() -> list[CandidateSolution], np.ndarray`
- `PySwarmsPSO(ingredients, ...)` → `.run() -> CandidateSolution`

注意：numpy/pymoo/pyswarms/dilipred 均惰性导入，仅在使用到相应功能时才需要环境已安装依赖。

---

## main.py

### CLI

- `--mode {hello,demo}`：默认 `hello`（无需依赖），`demo` 需安装优化与可视化依赖
- `--ingredients PATH`：成分文件路径（CSV/JSON）
- `--format {auto,csv,json}`：文件格式，默认 `auto`

### 关键函数

- `run_single_objective_algorithms(ingredients) -> (list[CandidateSolution], float)`
- `run_nsga(ingredients, toxicity_weight, highlight_points=None) -> None`
- `run_demo(ingredients: list[Ingredient] | None = None) -> int`
- `main(argv: list[str] | None = None) -> int`

示例：
```bash
python main.py --mode demo --ingredients resources/templates/ingredients.json --format json
```

---

## 附录：字段同义映射（节选）

- `name`: name/ingredient/compound/drug/名称/成分
- `smiles`: smiles/smile
- `hepatotoxicity_score`: hepatotoxicity/hepatotoxicity_score/dili/toxicity/toxicity_score
- `bbb`: bbb/logbb/bloodbrainbarrier（允许 0/1/true/false 或数值）
- 其他：`ob, dl, mw, alogp, h_don(hbd), h_acc(hba), caco2(logpapp/papp), fasa, hl(half_life/t_half)`
