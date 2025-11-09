# 配置说明（config/algorithms.json）

本文件汇总并细化项目的算法配置项说明，涵盖默认值、类型、语义与调优建议。

- 配置文件路径：`config/algorithms.json`
- 加载策略：在运行时读取后与内置默认值做“深合并”（仅覆盖你提供的键），其余沿用默认。
- 生效范围：GA/NSGA-II/PSO 进化算法与 ACI 指标计算。

## 默认配置（参考）

```json
{
  "ga": {
    "generations": 25,
    "population_size": 32,
    "toxicity_weight": 1.0
  },
  "nsga2": {
    "generations": 30,
    "population_size": 40,
    "max_candidates": null
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
  },
  "viz": {
    "enabled": true,
    "output_dir": "artifacts",
    "dpi": 200,
    "colormap": "viridis"
  }
}
```

说明：若 `config/algorithms.json` 缺失，上述即为运行时使用的内置默认。

## 配置项详解

### 1) GA（单目标，Pymoo GA）
- `generations`（int，默认 25）
  - 进化代数；越大搜索越充分，但耗时更长。建议 20–200。
- `population_size`（int，默认 32）
  - 种群规模；影响解的多样性与计算量。建议 16–128。
- `toxicity_weight`（float，默认 1.0）
  - 单目标适应度公式 `ACI - toxicity_weight * 肝毒性` 中的惩罚权重；范围 `[0, +∞)`。

### 2) NSGA-II（多目标，Pymoo NSGA-II）
- `generations`（int，默认 30）
- `population_size`（int，默认 40）
  - 与 GA 类似，数值越大越充分但更耗时。
- `max_candidates`（int|null，默认 null）
  - 将最终候选解下采样到最多 N 个。保持毒性维度的等间距采样以近似均匀分布；CLI 的 `--max-candidates` 优先于配置。

### 3) PSO（单目标，PySwarms PSO）
- `iterations`（int，默认 40）
  - 迭代轮数。建议 20–200。
- `swarm_size`（int，默认 24）
  - 粒子数量。建议 16–128。
- `options`（object）
  - `c1`（float，默认 1.2）：个体学习因子（cognitive）；偏大更易早收敛，偏小更偏探索。
  - `c2`（float，默认 1.2）：群体学习因子（social）；与 `c1` 共同决定收敛/发散趋势。
  - `w`（float，默认 0.6）：惯性权重；较大偏全球搜索，较小偏局部细化。常见 0.3–0.9。

### 4) ACI 指标相关
- `lambda_weight`（float，默认 0.3）
  - ACI = (1-λ)*平均可用性 + λ*互补性 中的权重 λ，范围 `[0,1]`。
- `complement_variance_max`（float，默认 0.25）
  - 多样性方差的归一化上限；越小越容易得到较高互补性分数。
- `coverage_floor`（float，默认 0.0）
  - 可用性分的“特征覆盖度”下限 `[0,1]`；用于惩罚特征缺失严重的样本。

## 调优建议
- 快速预览/小数据：下调 `generations/iterations` 与 `population_size/swam_size`。
- 强化探索（更充分）：上调 `population_size/swam_size` 与 `generations/iterations`；NSGA-II 更适合获取折衷解集合。
- 更重安全性（压低毒性）：上调 `ga.toxicity_weight`；或在最终加权选择时提升权重。
- 更重多样性：上调 `aci.lambda_weight`；或提升 `pso.options.w` 加强全局搜索。

## 局部覆盖示例
仅调整 GA 的种群规模：

```json
{
  "ga": { "population_size": 64 }
}
```

## 生效原理与回退
- 由 `algorithms._load_algorithm_config()` 读取并深合并。
- 当配置键不存在时，自动回退至内置默认值；当文件不存在时，整体回退至默认配置。
- 非法类型（例如将对象写成字符串）会引发运行期错误，请保持 JSON 结构与类型正确。
### 5) 可视化（viz）
- `enabled`（bool，默认 true）是否启用可视化产物。
- `output_dir`（str，默认 `artifacts`）输出目录。
- `dpi`（int，默认 200）图像分辨率。
- `colormap`（str，默认 `viridis`）帕累托前沿图的配色。
