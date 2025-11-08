# TODO

- [ ] 替换 `sample_ingredients` 为源自 TCMSP/SwissADME 的真实成分数据，保持 `synergy_baseline` 与 `hepatotoxicity_score` 有据可依。
- [ ] 将 NSGA-II 的非支配候选导出（CSV/JSON）并可视化其 Pareto 前沿，以便与 RAW264.7 IL-6 与 LD50 实验结果对照。
- [ ] 如果引入新的安全/协同指标（如 ADMET 其它字段或机器学习预测），在 `CombinationProblem` 中扩展目标数并调整 `diversity_penalty` 逻辑。
- [ ] 为 `algorithms.py` 增加更完整的单元测试（至少覆盖 `vector_to_candidate` 与评估函数），确保后续算法切换在 CI 中稳定。