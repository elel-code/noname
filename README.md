# 中药成分组合优化模板（含最小可运行入口）
## 开发与质量检查

本项目已内置 Ruff 规则，防止误用 `is/is not` 与字面量进行比较（F632）。

### 环境准备

1. 创建虚拟环境（Windows PowerShell）：
   - `py -3.11 -m venv .venv`
   - `.\.venv\Scripts\Activate.ps1`
2. 安装依赖（使用 uv 推荐）：
   - 仅主依赖：`uv sync`
   - 含开发依赖：`uv sync --group dev`

开发依赖（`ruff`、`pytest`）已在 `pyproject.toml` 的 `[dependency-groups.dev]` 中声明。

### 运行程序

- 基本验证（不依赖重库）：
  - `python main.py` → 输出 `Hello from noname!`
- 演示模式（需安装额外依赖并具备可写目录）：
  - 使用内置示例：`python main.py --mode demo`
  - 指定外部成分文件：
    - CSV：`python main.py --mode demo --ingredients resources/templates/ingredients.csv --format csv`
    - JSON：`python main.py --mode demo --ingredients resources/templates/ingredients.json --format json`
    - 自动按扩展名推断：`--format auto`（默认）

成分字段（列名不区分大小写，支持常见同义名）：
- 必需：`name`；且 `smiles` 与 `hepatotoxicity_score` 至少其一存在
- 可选：`ob, dl, mw, alogp, h_don, h_acc, caco2, bbb, fasa, hl`
  - `bbb` 可用 0/1、true/false 或 logBB 数值表示

## 接口速览（Python）

从 CSV/JSON 读取药物成分，并调用算法：

```python
from algorithms import load_ingredients, PymooNSGAII, select_candidate_by_weighted_sum

ingredients = load_ingredients("resources/templates/ingredients.csv", fmt="csv")

nsga = PymooNSGAII(ingredients)
solutions, metrics = nsga.run()
best, info = select_candidate_by_weighted_sum(solutions, ingredients, toxicity_weight=1.0)
print(info)
```

在脚本中直接运行演示：

```python
import main
main.main(["--mode", "demo", "--ingredients", "resources/templates/ingredients.json", "--format", "json"])
```

说明：依赖较重的库（numpy/pymoo/pyswarms/matplotlib/dilipred）均为惰性导入，仅在对应功能被调用时才需要安装。

### Lint 与自动修复

启用 Ruff 的 F632 规则（`is/is not` 与字面量比较）：

```
uv run ruff check . --select F632
uv run ruff check . --select F632 --fix
```

### 语法告警即错误

将 `SyntaxWarning` 视为错误进行编译检查：

```
Get-ChildItem -Recurse -Filter *.py | ForEach-Object { python -W error::SyntaxWarning -m py_compile $_.FullName }
```

### 运行测试

```
uv run pytest -q
```

测试包含：
- `tests/test_main.py`：校验 `--mode hello` 输出
- `tests/test_no_is_literal.py`：AST 扫描，确保仓库内无 `is/is not` 与字面量比较

该项目默认提供“最小可运行”的入口，便于在无依赖/只读环境中快速自检；
同时保留原有基于 `pymoo` 与 `pyswarms` 的优化演示作为可选模式。

## 快速开始（最小运行）

```bash
python main.py
# 输出：
# Hello from noname!
```

## 高级用法（演示模式，需要依赖且可写输出目录）

```bash
python main.py --mode demo
```

说明：演示模式会按需导入 `pymoo`、`pyswarms`、`matplotlib` 等依赖，并在 `artifacts/` 目录输出可视化图片。
若缺少依赖或在只读环境中运行，将给出友好提示而不会影响最小入口的可运行性。

## 项目结构

```
.
├── config/
│   └── algorithms.json  # 算法超参数与 ACI 配置
├── main.py              # 主入口：运行优化算法并输出结果
├── algorithms.py        # 核心逻辑：ACI 计算、优化问题定义与算法封装
├── ACI.md               # ADME 互补指数 (ACI) 的详细定义
├── API_LIST.md          # 主要函数与类的 API 文档
└── README.md            # 项目概览（本文）
```

## 依赖安装（可选）

如需运行演示模式，请在可写环境中安装相应依赖后重试（版本可根据实际需求调整）：

```bash
python -m pip install --upgrade pip
python -m pip install numpy pymoo pyswarms scikit-learn matplotlib dilipred
```

## ADME 互补指数 (ACI)

ACI 是一个综合指标，用于评估中药成分组合的整体质量，由两部分构成：

1. **单体可用性 (Desirability)**: 基于各成分的理化和 ADME 特征（如分子量、亲脂性、血脑屏障通透性等）计算得出的分数，反映了单个成分的“药用潜力”。
2. **组合互补性 (Complementarity)**: 衡量组合内各成分在 ADME 特征空间中的多样性。一个高度互补的组合意味着其成分在特征上各有所长，而非高度相似。

ACI 的计算公式如下：
`ACI = (1 - λ) * 平均可用性 + λ * 互补性`

其中 `λ` 是一个权重因子，用于平衡可用性与互补性的重要性（默认为 0.3）。

## 配置（config/algorithms.json）

所有算法相关的参数均在 `config/algorithms.json` 中定义；若文件缺失，则使用内置默认值。加载逻辑支持“深合并”，仅需覆盖需要变更的键即可（参见 `algorithms._deep_update`）。

示例（与当前默认一致）：

```
{
  "ga": {
    "generations": 25,
    "population_size": 32,
    "toxicity_weight": 1.0
  },
  "nsga2": {
    "generations": 30,
    "population_size": 40
  },
  "pso": {
    "iterations": 40,
    "swarm_size": 24,
    "options": { "c1": 1.2, "c2": 1.2, "w": 0.6 }
  },
  "aci": {
    "lambda_weight": 0.3,
    "complement_variance_max": 0.25,
    "coverage_floor": 0.0
  }
}
```

参数说明（按模块）：

- `ga`（单目标遗传算法，Pymoo GA）
  - `generations`（int）：进化代数；越大搜索越充分、耗时越长。建议 20–200。
  - `population_size`（int）：种群规模；影响多样性与计算量。建议 16–128。
  - `toxicity_weight`（float）：单目标适应度 `ACI - toxicity_weight * 肝毒性` 的惩罚权重。范围 [0, +∞)。

- `nsga2`（多目标，Pymoo NSGA-II）
  - `generations`（int）：进化代数。建议 20–200。
  - `population_size`（int）：种群规模。建议 32–200。

- `pso`（单目标，PySwarms PSO）
  - `iterations`（int）：迭代轮数。建议 20–200。
  - `swarm_size`（int）：粒子数量。建议 16–128。
  - `options`（object）：PSO 超参数
    - `c1`（float）：个体学习因子（cognitive）；偏大易早收敛，偏小探索更多。常见 0.5–2.5。
    - `c2`（float）：群体学习因子（social）；与 `c1` 一起决定收敛/发散趋势。常见 0.5–2.5。
    - `w`（float）：惯性权重；较大更关注全局搜索，较小更精细局部搜索。常见 0.3–0.9。

- `aci`（指标权重与归一化）
  - `lambda_weight`（float）：ACI = (1-λ)*平均可用性 + λ*互补性 的权重 λ，范围 [0,1]；默认 0.3。
  - `complement_variance_max`（float）：互补性方差归一化上限；越小越容易获得较高互补性分数。默认 0.25。
  - `coverage_floor`（float）：可用性分的特征覆盖度下限（[0,1]），用于惩罚特征缺失。默认 0.0。

使用建议：
- 小数据/快速预览：降低 `generations/iterations` 与 `population_size/swam_size`。
- 更重探索：增大 `population_size/swam_size` 与 `generations/iterations`；NSGA-II 更适合探索折衷解集。
- 偏向安全性：增大 `ga.toxicity_weight`；或在最终选择时提高权重。
- 偏向多样性：增大 `aci.lambda_weight`；或提高 `pso.options.w` 增强全局搜索。

只覆盖想变更的键即可，例如仅调 GA 的种群规模：

```
{
  "ga": { "population_size": 64 }
}
```

## 示例输出

运行 `main.py` 将依次执行单目标优化（GA 和 PSO）与多目标优化（NSGA-II），并输出类似如下的结果：

```
--- 单目标优化 (GA) ---
[GA] ACI (修正后): 0.85, 肝毒性: 0.25, 惩罚: 0.0
- 黄芩素 (Baicalein): 30%
- 黄连碱 (Berberine): 70%

--- 单目标优化 (PSO) ---
[PSO] ACI (修正后): 0.88, 肝毒性: 0.22, 惩罚: 0.0
- 黄芩素 (Baicalein): 40%
- 黄连碱 (Berberine): 60%

--- 多目标优化 (NSGA-II) 非支配集 ---
候选 1: ACI=0.92, 肝毒性=0.35
候选 2: ACI=0.85, 肝毒性=0.20
...

--- NSGA-II 帕累托前沿加权选择 ---
[NSGA-II 择优] ACI (修正后): 0.85, 肝毒性: 0.20, 惩罚: 0.0
- 黄芩素 (Baicalein): 35%
- 黄连碱 (Berberine): 65%
```
