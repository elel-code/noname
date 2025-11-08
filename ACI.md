用仅基于理化/ADME 特征的“ADME 互补指数（ACI）”来替代“协同指数”。它不声称药效学协同,只用列的字段（MW, AlogP, HDon, HAcc, OB, Caco-2, BBB, DL, FASA-, HL），计算每个成分的“可用性/潜力”分，再在组合层面叠加一个“互补性（多样性）”项作为类协同替代。

指标定义（简单好实现）

- 单体可用性分（Desirability，D_i\in[0,1]）：把每个理化/ADME特征映射到 0–1 的“宜居区”分值，然后取几何平均（避免某项极差被平均掩盖）。
- 组合互补项（Complementarity，C\in[0,1]）：看组合在这些特征上是否“互补/不全都挤在同一端”。实现上用Shannon 熵或特征 z 分的均方差做一个 0–1 归一化的多样性分。
- ACI（越大越好）：

\boxed{\text{ACI}= (1-\lambda)\cdot \overline{D}\ +\ \lambda\cdot C}

其中 \overline{D} 是（配比加权）单体分的平均；\lambda\in[0,1] 默认 0.3（保留可用性为主，互补为辅）。

把 ACI 当“替代协同指数”用；另一指标固定为你自算的肝毒性得分即可组成双目标（max ACI, min HepatoRisk）。
