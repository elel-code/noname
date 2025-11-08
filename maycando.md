把 aci 与“母药单一性惩罚”分离暴露

现在 evaluate_metrics 里把 ACI 直接乘了 (1 - penalty) 返回（aci_with_penalty）。建议同时返回“纯 ACI”与“惩罚后 ACI”，便于做消融分析/可视化；多目标时也能把“惩罚”做第三目标或硬约束（如 majority < 0.8）。
（在 evaluate_metrics 里多返回一个 aci_raw 即可。）

algorithms

让 ACI 的互补项权重可配置

ACI_LAMBDA=0.3 写死在代码里；建议让 compute_aci_score(..., lambda_weight=...) 的默认值读自配置（和 GA/NSGA 参数一处管理），方便做灵敏度分析。

algorithms

缺失特征的“覆盖度惩罚”

现在 compute_desirability_score 对缺失值直接跳过，几何平均会抬高只有少量特征的成分得分。建议乘上一个覆盖度系数：
coverage = (#有效特征)/(#总特征)，最后返回 desirability ** 1.0 * coverage（或把 coverage 线性混入）。

algorithms

互补性方差的尺度参数暴露

互补性里用 variance / 0.25 截断到 1，其中 0.25 是 [0,1] 变量的最大方差；建议把 0.25 作为 COMPLEMENT_VAR_MAX 常量暴露到配置，便于调平滑度。

algorithms

Caco-2、BBB 的取值约定

你的窗口预设 caco2: 0.4–1.5、bbb: -1~1 很合理，但不同数据源口径差异大。建议：

在 FeatureWindow 上增加 note / unit 字段或在 README 里写明口径（logPapp? logBB?）。

若 BBB 为 0/1 分类，当前 _window_score 会把 0/1 当连续值处理；可以在读入数据时先映射成 -0.5 和 +0.5，或对 BBB 单独写个离散映射函数。

algorithms

单目标的权重管理

single_objective_score = ACI - toxicity_weight * toxicity；建议把 toxicity_weight 也写进 config/algorithms.json（你已支持），并提供一个从多目标帕累托前沿到加权和的便捷小函数，便于“赛后”按不同权重挑解。

algorithms

命名与注释的小修正

Ingredient.hepatotoxicity_score 的注释写了 “SwissADME 肝毒性打分”，但 SwissADME 本身不出肝毒性；建议把注释改成“外部模型的肝毒性打分（越高越差）”，保持来源中立。