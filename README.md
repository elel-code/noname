# noname ACI 模板

该项目演示如何基于 pymoo/pyswarms 在中药成分集合上优化「可用性 + 肝毒性」双目标。

## ACI 特征窗口口径

- **Caco-2**：默认使用 `log10(Papp, cm/s)`。若数据为原始 Papp，可先取对数；无法提供时该字段视为缺失。
- **BBB**：默认使用 logBB；若源数据为 0/1 分类，将自动映射为 -0.5/0.5 以兼容连续窗口。
- 其它特征（MW、ALogP、HDon/HAcc、OB、DL、FASA-、HL）均按 FeatureWindow 中的 `note` 说明提供单位，缺失值会触发覆盖度惩罚。

## 配置要点

`config/algorithms.json` 现新增 `aci` 区块，可调节：

- `lambda_weight`：互补项权重（默认 0.3）。
- `complement_variance_max`：互补性归一化时的方差上限（默认 0.25）。
- `coverage_floor`：特征覆盖度的下限，避免极端缺失时得分被完全抹平。

GA/NSGA-II/PSO 的 generations/swarm 等参数仍与原项目一致，单目标的 `toxicity_weight` 会被 `select_candidate_by_weighted_sum` 复用，用于从 Pareto 集挑选目标解。