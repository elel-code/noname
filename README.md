# 中药成分组合优化模板

该项目演示了如何基于 `pymoo` 和 `pyswarms` 在中药成分集合上进行双目标优化，旨在最大化 **ADME 互补指数 (ACI)** 并最小化 **肝毒性**。

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

## 如何运行

1. **安装依赖**:
   ```bash
   pip install -r requirements.txt
   ```
   *注意：`requirements.txt` 文件当前缺失，需要根据 `pyproject.toml` 或 `uv.lock` 手动生成。*

2. **执行优化**:
   ```bash
   python main.py
   ```

## ADME 互补指数 (ACI)

ACI 是一个综合指标，用于评估中药成分组合的整体质量，由两部分构成：

1. **单体可用性 (Desirability)**: 基于各成分的理化和 ADME 特征（如分子量、亲脂性、血脑屏障通透性等）计算得出的分数，反映了单个成分的“药用潜力”。
2. **组合互补性 (Complementarity)**: 衡量组合内各成分在 ADME 特征空间中的多样性。一个高度互补的组合意味着其成分在特征上各有所长，而非高度相似。

ACI 的计算公式如下：
`ACI = (1 - λ) * 平均可用性 + λ * 互补性`

其中 `λ` 是一个权重因子，用于平衡可用性与互补性的重要性（默认为 0.3）。

## 配置

所有算法相关的参数均在 `config/algorithms.json` 中定义。

### ACI 相关参数

- `lambda_weight`: ACI 公式中的 `λ` 权重。
- `complement_variance_max`: 用于归一化互补性分数的方差上限。
- `coverage_floor`: 特征覆盖度的下限，用于惩罚数据缺失过多的成分。

### 优化算法参数

- **GA/NSGA-II**: `generations`, `pop_size`
- **PSO**: `n_particles`, `iters`, `options`

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
