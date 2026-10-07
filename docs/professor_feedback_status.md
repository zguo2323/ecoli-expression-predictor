# 教授评语与当前回应进度

**记录日期：** 2026-10-06（根据用户最新情况更新）
**项目分支：** `feat/biological-context-phase-1`

## 教授的总体评价

教授肯定项目不仅是 wrapper 或 prompt 工具：它基于真实实验数据训练预测模型，提供 −10/−35 motif、spacer、Shine–Dalgarno score、SD spacing、GC 和 mRNA folding 等可解释特征，也有面向设计的 RBS 排名输出。主要批评是生物学解释及结果表征还不充分：旧结构特征采用 promoter 3′ 端加 RBS，没有纳入 RBS–CDS junction 和早期 CDS；仅报告 Spearman 不足以支撑强结论，还需要 baseline、feature importance、error analysis，以及整组留出 promoter 或 RBS 的困难拆分。

## 逐项回应

| 教授指出的问题 | 当前实现或证据 | 状态与剩余工作 |
|---|---|---|
| folding 上下文缺少 RBS–CDS junction 和早期 CDS | 已依据测得的 `TSS.best` 重建转录前缀至起始 ATG，并拼接暂定 sfGFP 起始 90 nt；模型计算 SD/start 未配对概率与 SD 至 start opening ΔG。B6 在 `log(prot/RNA)` 上相对纯序列 B4 的 5 折增益见 [`summary.md`](../reports/biological_context/summary.md)：promoter、RBS、double-unseen 的 Spearman 配对差 95% CI 分别为 `[0.088, 0.148]`、`[0.041, 0.125]`、`[0.059, 0.178]`；log-MAE 差区间均低于 0。 | **实现并有分组证据；确切 CDS 核验因客观注册限制跳过。** 用户已毕业，Addgene 注册流程要求填写其目前无法提供/选择的组织信息，因此无法创建可访问序列的账号。按用户决定，不再执行该项。保留论文来源的暂定 CDS，不能称为实验质粒的精确序列；在解释与报告中持续披露该限制。 |
| 缺少 baseline 和算法选择证据 | 已有训练集常数/元件 baseline、B3/B4 消融、OSTIR 外部机制 baseline 和同特征/同折 Ridge；XGBoost/Ridge 配对 bootstrap 也已完成。 | **已明显补强。** XGBoost 在 promoter/RBS 留出翻译 proxy 指标优于 Ridge；双新元件差值 CI 跨 0，不足以判定胜者。OSTIR 使用 ViennaRNA 2.7.2、超出 OSTIR 1.1.1 最后验证版本 2.4.18，需保留版本限制。 |
| 缺少 feature importance | 已用 B6 在 4 种场景 × 5 折 × 2 个目标上做 biological-unit held-out permutation：promoter、RBS 和 promoter×RBS 单位分别置换对应特征组，每折每组 30 次。定向测试 [`test_grouped_feature_group_permutation.py`](../tests/test_grouped_feature_group_permutation.py) 检查组覆盖、单位内恒定、置换一致性和重要性方向；全量输出 3,600 条 repeat-level 记录。 | **明显补强。** 结果在 [`grouped_feature_group_permutation_summary.csv`](../reports/biological_context/grouped_feature_group_permutation_summary.csv)，解释和主结果见 [`summary.md`](../reports/biological_context/summary.md) 的 feature-group permutation 小节，图为 [`grouped_feature_group_permutation.png`](../reports/biological_context/grouped_feature_group_permutation.png)。在 `log(prot/RNA)` 上，promoter holdout 最依赖 RBS/SD（Spearman drop / MAE increase **0.346±0.066 / 0.277±0.041**）；RBS holdout 最依赖 promoter（**0.347±0.027 / 0.278±0.031**）；double-unseen 两类序列特征接近（promoter **0.168±0.093 / 0.113±0.070**，RBS **0.168±0.094 / 0.119±0.060**）。Accessibility 组置换后各场景误差也都上升，说明模型使用了这组结构特征。数值是固定五折上的平均模型依赖，不能解释为因果或分子机制。 |
| 缺少 error analysis | 已完成 B6、同特征 Ridge、OSTIR 的 OOF 残差分析：按 translation-proxy 实测五分位、4 种 held-out 场景、promoter/RBS 及 double-unseen pair 汇总；并保存每场景/目标/模型绝对残差最大的 20 个构建。定向测试 [`test_grouped_oof_error_analysis.py`](../tests/test_grouped_oof_error_analysis.py) 验证共享分箱、component/pair 分组、配对比较、eta-squared 和极端案例方向。 | **明显补强了错误位置和模式的记录。** 可查看 [`grouped_oof_error_extreme_cases.csv`](../reports/biological_context/grouped_oof_error_extreme_cases.csv)（400 行，含 scenario、target、model、构建 ID、promoter/RBS ID、actual、predicted、residual、绝对误差和方向），以及 [`grouped_oof_predictions.csv`](../reports/biological_context/grouped_oof_predictions.csv)（完整 OOF 逐构建预测）。具体例子：RBS holdout 的 salis-3-3 上 B6 MAE 1.096、Ridge 0.679（77 个构建）；salis-3-11 为 1.675 vs 1.423（112 个构建）。图表见 [`grouped_oof_error_analysis.png`](../reports/biological_context/grouped_oof_error_analysis.png)。这些是可追溯的失败案例记录；当前目标是让误差透明、便于复查，不要求逐例证明某个确定生物机制。 |
| 验证拆分不够困难 | 已完成 5 折随机诊断、整组 promoter holdout、整组 RBS holdout 和 promoter×RBS 双重未见拆分；同时做了 1,000 次 group bootstrap。double-unseen 训练行排除测试 promoter 和 RBS，交叉配对按 unused 处理。 | **已明显补强。** 结果仍来自同一实验数据和一个固定 seed；bootstrap 区间以这些固定外层预测为条件，不涵盖重训或新实验数据集的变异。 |
| 目前结果不足以支持顶尖评分 | 生物学结构假设、外部/线性 baseline、严格组留出、配对不确定性、跨折分组特征重要性、OOF 残差诊断和 4 个外层 split seed 的 grouped validation 已形成更完整证据链。 | **明显补强但仍有边界。** B6 相对 B4 的翻译 proxy 改善在 promoter/RBS/double-unseen 三类生物 holdout 均 4/4 seeds 同方向；B6 与 Ridge 在双新元件场景仍无稳定胜者。仍是同一实验数据集，不能代替独立数据验证；Addgene 序列核验按用户决定跳过，reporter CDS 暂定。详情见 [`multi-seed summary`](../reports/biological_context/multiseed/multiseed_cross_seed_summary.csv) 和阶段日志阶段 8。 |

## pJ251-GERC 序列核验：需要用户协助的操作

本项目原计划比对 Addgene plasmid **#47441, pJ251-GERC** 的序列。Addgene 当前页面显示 1 条 depositor full sequence 和 3 条 Addgene-verified partial sequences；页面提示登录后才可下载或复制。用户尝试注册时遇到组织信息要求；用户已毕业，目前无法提供/选择符合要求的组织信息，因此不能注册访问。**用户决定将这一步作为客观限制跳过，不再继续尝试。** 相关公开页面记录见 [Addgene sequence page](https://www.addgene.org/47441/sequences/) 和 [Addgene help: sequence access](https://help.addgene.org/hc/en-us/articles/44210206209549-Are-there-any-requirements-for-viewing-plasmid-sequence-data-on-the-Addgene-website)。

**当前状态：跳过。** 原操作方案仅作为将来出现可访问文件时的处理记录：若用户之后通过其他合法渠道获得 GenBank/FASTA，可放到 `data/reference/addgene_47441/` 或上传；再定位 sfGFP ORF、方向与 ATG，并将前 90 nt 与暂定序列比对。Addgene depositor full sequence 不等同于 Addgene 全质粒验证，partial verified sequence 也只能核对其覆盖区域；与论文实验中实际使用的物理质粒是否一致仍需单独证据链。

## OOF 残差分析：已完成描述性诊断

在相同 grouped OOF 预测上比较 B6、Ridge 与 OSTIR，按实测目标五分位统计 signed residual/MAE；按 scenario、promoter、RBS 和 double-unseen promoter×RBS pair 汇总 component errors，并导出最大残差案例。特别检查了 RBS holdout 中 promoter-recognition feature-group permutation effect 较大的现象：B6 的组内 MAE 优于 Ridge 的比例为 108/112 个已知 promoter 背景（96.4%）和 91/111 个留出 RBS（82.0%），说明 B6 的总体优势分布较广，并非只来自少量 promoter。promoter holdout 中 B6 对 Ridge 的胜组率降为 promoter 72.3%、RBS 67.6%；double-unseen promoter×RBS pair 为 47.5%，与严格留出下两模型整体接近相符。与此同时，目标边际 eta-squared 在 RBS holdout 的 promoter 上为 0.346、RBS 上为 0.385；双新元件分别为 0.370、0.386。因此 permutation 排名不能简单归结为目标在 promoter 间变化更大。Eta-squared 是单因素边际摘要，未分解交叉设计和交互；重要性依然表示模型依赖，不是生物因果解释。

残差还显示所有比较方法在翻译 proxy 低值段偏高估、高值段偏低估，double-unseen 两端最明显。该规律可由范围压缩、measurement noise、未建模生物学等多个因素产生，目前没有足够证据确定来源。下一步应优先多随机种子/重复 grouped split，或基于有来源的序列和实验元数据复核具体高误差组；不应将当前描述性分层写成机制结论。

详细审计记录见 [`CDS source audit`](../reports/biological_context/cds_source_audit.md)。

## 预测区间字段与独立数据验证方向

旧 `predict_expression` 曾返回点预测 ±25% 的 `confidence_interval`。这是早期接口规格中的固定占位值，没有统计校准依据，容易让使用者误以为是可靠不确定性区间；阶段 9 已从实现、JSON schema、SPEC、Skill 说明和测试中移除。现在 API 返回 `predicted_prot` 点预测、`features_used` 与 `translation_context_available`。模型比较用的 group-bootstrap 95% CI 是另一项由重采样得出的评测结果，仍保留。

可考虑外部验证的数据候选已记录在 [`biological_modeling_stage_log.md` 阶段 9](biological_modeling_stage_log.md)：Wen et al. 2024 的 520 个调控序列×目标基因 GFP-fusion fluorescence 优先检查序列边界和数据列；Dash et al. 2024 Dryad 的 16 个 native promoter::GFP 数据作为 promoter-only 备选；Bonde et al. EMOPEC 的估计表达表需区分实测与估计值。外部验证尚未运行。
