# 生物学解释优化与评测方案

## 目标

把项目的 mRNA 结构特征从“promoter 末端 + RBS 的 MFE”改成与翻译起始机制相对应的 5′ mRNA 起始区可及性，并用严格的分组留出和有机制含义的 baseline 证明它带来可重复的提升。主要科学问题应表述为：在控制转录量后，序列特征能否预测翻译效率（protein/RNA proxy），特别是对训练时没见过的 promoter 或 RBS 是否成立？

## 仓库现状与需要纠正的解释

- `modules/features/features.py::compute_mrna_folding_energy` 把 `promoter_seq` 最后 30 nt 和 `rbs_seq` 拼起来，再用 `RNA.fold` 的单一最低自由能（MFE）。但 promoter DNA 不等于转录出的 RNA；只有 TSS 之后的实际转录区才能作为 RNA 输入。全序列 MFE 也不能直接表示 SD 和起始密码子的可及性。
- `modules/pipeline/pipeline.py::load_promoter_table` 读入 `TSS.best`，而 `build_dataset` 只保留 promoter ID 和 sequence，因而特征层拿不到 TSS 信息。
- `sd02.xls` 的 RBS sequence 包含 RBS/leader 和末端的 `CATATG` 起始区序列，但不含 reporter CDS 的下游上下文。该研究用这些调控元件控制同一个 sfGFP reporter，因此可尝试从论文的 Addgene 构建（质粒 47441）或补充材料取回精确 reporter 序列；必须核对拼接边界、方向及 TSS 坐标后才能用。
- ViennaRNA 当前不可用时 folding 特征静默返回 `0.0`。这会把“没算出来”伪装成一个有效数值，需改成显式缺失/报错状态，并在结果中记录 folding 引擎和参数。
- 当前普通随机拆分会让相同 promoter 和 RBS 身份同时出现在训练集、测试集。这个结果衡量的是已见元件组合的插值能力，不能证明对新 promoter/RBS 的泛化。

## 优化方案

### 阶段 1：重建真实转录本上下文

1. 保留 `TSS.best`，同时保留原始表行和构建 ID，追溯它相对 promoter 序列的坐标定义。用论文、补充材料及 Addgene 47441 的序列/构建设计核验坐标和序列边界；不要仅凭字段名假定正负号方向。
2. 逐个 promoter 构建 TSS 后的 5′ leader RNA。把 leader、RBS 序列中的实际转录部分、起始密码子以及同一 sfGFP ORF 的开头拼接起来。对至少若干有代表性的构建做手工坐标核对并记录来源。
3. 明确边界约定：RBS 序列中是否含起始密码子、是否包含 BioBrick/限制性位点 scar、起始密码子是哪 3 nt、CDS 从哪里开始。拆掉非转录 DNA 与重复碱基，不对 promoter −10/−35 区做 RNA 折叠。
4. 若无法从权威来源恢复准确 TSS 或 reporter CDS，就先提供“只基于已确认 RBS 区序列”的结构特征并标注局限；不要将猜测的 promoter tail/CDS 输入描述为真实 mRNA。

### 阶段 2：结构特征改为起始区可及性

以翻译起始位点为坐标中心，至少实现下列可解释特征：

- SD 区域（按真实 RBS 注释坐标）未配对概率/平均可及性；
- 起始密码子未配对概率；
- SD + spacer + start 区域的联合开放能量（打开指定区段所需能量）；
- 可选：起始位点前后局部窗口（例如 −30 到 +30 nt）的 MFE，作为与现有特征对照的描述量，而不是把它当作可及性的唯一指标。

优先采用 ViennaRNA partition function/`RNAplfold`/opening energy 相关接口；MFE 可保留作 ablation。特征需固定并记录温度、RNA 参数版本、窗口和最大 base-pair span。单序列 thermodynamic prediction 是近似值，不应称作体内真实结构。比较 OSTIR 的 initiation-rate prediction 作为外部机制 baseline；OSTIR 会把起始区结构与翻译起始相互作用纳入模型。老版 Salis RBS Calculator 可用于方法参考，但其可运行性、许可证和软件版本需先核实，不把不可复现网页结果当主要 benchmark。

### 阶段 3：把预测目标分成转录与翻译问题

项目现有 `prot` 目标是总蛋白输出，包含启动子控制的 RNA 丰度和 RBS 控制的翻译效率。建议主要报告两个相互补充的任务：

1. **端到端表达**：预测 `log(prot)`，对应当前产品目标。
2. **翻译效率 proxy**：预测 `log(prot) - log(RNA)`，或等价的 `log(prot/RNA)`，仅在 RNA、protein 都为正且数据可用时计算。明确这是稳态 protein/RNA proxy，不等同于纯化学的瞬时翻译速率；对零值/低计数值预先规定过滤或检测下限处理，并做敏感性分析。

另外可单独预测 `log(RNA)`，确认 promoter 表征主要作用在转录任务，RBS 与起始区结构特征是否主要改善翻译任务。所有从原始表得到的 `mean.*`、`model.*`、`deltaG` 列先查清其来源；除明确作为“已测元件特性”baseline 外，不得将可能使用了测试构建结果的字段喂给模型。

## Baseline 与模型变体

所有方法共用完全相同的外层数据划分和 target 处理。

| ID | 方法 | 用途 |
|---|---|---|
| B0 | 常数/训练集均值或中位数 | 无序列信息下限 |
| B1 | promoter 与 RBS 的独立边际效应相加（只用训练集拟合元件效应） | 检验组合矩阵是否只需元件身份记忆 |
| B2 | 数据集内独立测得的 promoter `mean_RNA` × RBS `mean.xlat` 组合规则（仅作有测量数据的 oracle/工程 baseline，不能用于未知元件预测） | 对照 Kosuri 式 composability |
| B3 | 现有 10 个特征 + XGBoost | 当前代码 baseline；修复输入缺失/不可用值后重跑，锁定实现版本 |
| B4 | 序列模型，不含 folding（promoter motif/spacing、SD/spacing、GC 等） | 隔离结构贡献 |
| B5 | B4 + 当前 MFE（promoter tail + RBS） | 显示当前结构假设是否有用，预期仅作为负/对照实验 |
| B6 | B4 + 真实 5′ mRNA 起始区的 accessibility/opening 特征 | 核心改进 |
| B7 | OSTIR 对同一完整起始区序列的 initiation rate | 机制型外部 baseline；另可将 log(OSTIR) 与 promoter RNA 预测组合为总蛋白预测 |

B1 的元件效应使用训练集拟合，未见元件的预测回退到训练均值；同时分别报告 seen/unseen 情况。B2 的测量值不适用于部署预测，只用于回答“如果独立 characterization 已知，简单组合能达到什么水平”。

## 分组拆分与主要指标

不要把普通随机 split 当主要结果。固定 5-fold 外层 group CV 或预先冻结的 group train/validation/test，并按 `promoter_id`、`rbs_id` 分别成组。超参数与任何阈值只能在训练/验证组内选。

| 测试场景 | 拆分方法 | 回答的问题 |
|---|---|---|
| 随机组合留出（仅诊断） | 行级随机 split | 已见元件的组合插值，便于复现旧分数 |
| 新 promoter | 整个 promoter 留出 | 对未见启动子能否泛化，RBS 多为已知 |
| 新 RBS | 整个 RBS 留出 | 对未见 RBS 设计能否泛化，promoter 多为已知 |
| 双新元件 | promoter 和 RBS 均不与训练身份重叠 | 最严格的组合外推；单独呈现，避免与前两类平均掉 |

主指标：翻译 proxy 的 Spearman ρ、log-scale MAE、R²；总蛋白也报告相同指标。RBS 设计任务增加每个 promoter 内的 Kendall τ/Spearman、top-k 命中率（实测前 10% 是否被选入预测前 k）、top-1 regret（实测最优值与被选序列实测值的差，建议在 log scale）。置信区间按 promoter/RBS group bootstrap，而不是把相关构建行当独立样本 bootstrap。报告各折分布和每类测试组的 n。

## “生物学提升”的直观证据

需要同时看到预测变好和机制行为符合预期，不能只凭 XGBoost feature importance 宣称解释性：

1. **结构特征增益**：逐折比较 B4、B5、B6 在 held-out RBS 及双新元件上的 translation proxy 表现；报告 ΔSpearman、ΔMAE 及 group-bootstrap 置信区间。B6 相对 B4 的增益是结构特征价值的主要证据。
2. **方向性关联**：在控制 promoter/RBS 身份或使用组内差异后，SD/start 更开放是否对应较高的 `prot/RNA`。画按 accessibility 分箱的实测翻译 proxy 曲线与置信区间，并显示样本量。
3. **配对反事实检查**：对相同 promoter 的 RBS 变体，检查模型能否预测实测翻译效率的排序；再按预测 accessibility 的升降分组，比较实际 `prot/RNA` 变化。此分析利用了组合矩阵的成对设计。
4. **消融和归因稳定性**：按 promoter、RBS、结构、组成等生物学特征组做 drop-group ablation/置换重要性，跨外层折汇总。高度相关的结构特征不能按单个树分裂次数排名解释。
5. **误差分析**：列出每个外层折最差预测对，按低计数/测量噪声、异常高表达、feature 边界、强结构与高表达不一致等类别复核；对 `prot` 与 `RNA` 分别画预测-实测图及残差分布。
6. **模型适用范围**：报告数据集固定的 sfGFP、宿主/培养/载体条件。没有跨 CDS 数据时，结论限于该 reporter 起始上下文；不能声称已验证任意 CDS/RBS 泛化。

建议最终展示一张主图：四个留出场景下 B3/B4/B6/B7 的翻译 proxy Spearman 点图（含 group-bootstrap CI），旁边配一张实测 `prot/RNA` 随起始区 accessibility 变化图。这样能一眼看出增益来自对生物起始区的更好刻画，并展示严格留出下是否还成立。

## 推荐执行顺序与交付

1. **数据语义审计**：把 TSS 坐标、promoter/RBS 边界、reporter CDS 来源写成可追溯的数据字典；若不能无歧义重建，暂停结构特征扩展。
2. **基线重跑**：冻结原始代码、feature 可用性和随机/分组 split；建立同一套折分配与 target 构造。
3. **结构特征实现**：添加 transcript-context builder 和 accessibility 特征，保存 RNA 序列及区间坐标方便人工核查。
4. **OSTIR 对照和消融**：相同转录本输入，不同数据划分；记录版本、参数和可复现命令。
5. **误差分析与报告**：输出可重复的 CSV/JSON 指标、逐构建预测和图表，README 报告旧随机 split 与新 group split 的差异，不只保留最好的分数。

建议交付文件：`data_dictionary.md`、可复用的 split ID 表、每个构建的 held-out prediction 表、模型/特征消融指标表、误差分析图，以及包含环境版本与重跑命令的实验说明。

## 后续验证清单（阶段 4）

阶段 1–3 已完成初步实现和固定拆分 benchmark。阶段 4 的 OSTIR 全量对照、5 折 grouped CV 和 group bootstrap 已运行；CDS 权威序列核验尝试受 Addgene 登录限制，仍未解决。结果和失败记录见 [`biological_modeling_stage_log.md`](biological_modeling_stage_log.md) 与 [`benchmark summary`](../reports/biological_context/summary.md)。

1. **OSTIR/ViennaRNA CLI 环境：已完成并留版本提示**。`ostir-b7` 环境内依赖可用，已输出全量 B7 结果；ViennaRNA CLI 2.7.2 高于 OSTIR 1.1.1 最后验证的 2.4.18，严格版本复现仍可后续检查。
2. **分组交叉验证：已完成**。保存行随机、promoter、RBS、double-unseen 的五折分配与逐折预测，目标为 `log(prot/RNA)` 和 `log(prot)`。
3. **分组 bootstrap：已完成**。1,000 次、固定 seed；保存绝对指标和 B6−B4/B7−B6 配对区间。
4. **reporter CDS 来源：核验未闭环**。公开 Addgene 页列出 full depositor sequence 和 verified partials，但序列数据要求登录；暂定 90 nt 未被逐碱基确认，见 [CDS source audit](../reports/biological_context/cds_source_audit.md)。获得有授权的 GenBank/FASTA 后应完成比对和敏感性分析。

阶段 4 的交付材料：环境版本记录、OSTIR 全量预测文件、固定 fold ID 文件、逐折及 group-bootstrap 指标表、CDS 来源核验记录，以及更新后的主图/summary。当前 CDS 来源限制和 OSTIR 版本提示均已记入阶段日志。

## 算法选择验证：Ridge 已完成，其他模型为备选

阶段 5 已完成 XGBoost 与 Ridge 的同条件比较：使用相同 B3/B4/B6 特征、完全相同的外层 5 折及两个目标；Ridge 在各外层训练折内使用 StandardScaler，并用与外层分组结构匹配的 3 折 inner CV 选择 alpha。训练、测试身份隔离和配对 group bootstrap 见阶段日志及 benchmark summary。结果支持 XGBoost 在随机、单组留出上明显更强；在严格 double-unseen 上，B6 Ridge 与 B6 XGBoost 的指标区间重叠，不能把模型优劣泛化到该场景。

阶段 6 已完成 B6 的 held-out biological feature-group permutation：promoter 与 RBS 特征分别在各自 ID 单位置换，转录本可及性特征按 promoter×RBS 对联合置换；每个外折、目标和组 30 次。主要翻译 proxy 结果显示三个特征组在各生物学 holdout 均有正向模型依赖，但重要性排名随场景变化；详细数值见 [`summary.md`](../reports/biological_context/summary.md) 与阶段日志阶段 6。阶段 7 已完成描述性 OOF 残差分析：比较 B6、Ridge、OSTIR 在多种生物留出中的误差，检查 RBS holdout promoter-feature permutation pattern，按元件身份、实测值区间汇总残差并保存逐构建和极端案例记录。目标是透明呈现误差，不要求对每例证明确定的生物成因。阶段 8 已完成 seed 42、7、123、2026 的 grouped CV split-stability 检查；B6 相对 B4 在三个生物 holdout 均 4/4 seeds 同方向改善，B6 对 Ridge 的 double-unseen 优势仍不稳定。当前数据集相同，外层 split seed 检查不能代替独立实验验证。pJ251-GERC 权威序列核验因 Addgene 注册条件客观无法完成，按用户决定跳过，当前 CDS 仍为 provisional。

产品接口阶段 9 移除了未经校准的 ±25% `confidence_interval` 字段，只返回点预测、使用特征和 TSS 上下文是否可用；详情及测试见阶段日志。下一项可选加强方向是独立外部数据验证：优先检查 Wen et al. 2024 的 520 个 regulatory sequence×gene GFP-fusion measurements 是否能恢复序列边界和实验标签；另有 Dash et al. 2024 Dryad 的 16 个 E. coli promoter::GFP reporter 适合较窄的 promoter 排序验证。当前只是确认了公开候选，没有把外部测试结果冒充为已经完成。

后续可将 ElasticNet、Random Forest 或 ExtraTrees 纳入同一实验框架；保持特征、外层折、目标和训练内调参规则一致，并继续报告多折 Spearman、log-scale MAE、折间波动和配对组 bootstrap。Ridge 标准化系数可作线性模型描述，但受相关特征影响，不应单独作因果解释。

## 文献与开源实现

- Kosuri et al. 2013，项目原始组合数据（12,563 promoter×RBS combinations；同时测 RNA 和 protein，研究也讨论 mRNA 结构对 translation 的影响）：[PNAS/PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3752251/)，[Addgene plasmid 47441 信息由论文记录](https://pubmed.ncbi.nlm.nih.gov/23924614/)。
- Salis, Mirsky & Voigt 2009，RBS Calculator 的热力学翻译起始模型将起始位点周围 mRNA 上下文和 RNA/核糖体相互作用纳入预测：[论文全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2782888/)。
- Roots et al. 2022, OSTIR，开放的 translation initiation rate 预测实现与论文：[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9518832/)，[GitHub](https://github.com/barricklab/ostir)。
- ViennaRNA 官方仓库与 Python API 支持 partition function、base-pairing probability 等 ensemble 计算，不限于单个 MFE：[GitHub](https://github.com/ViennaRNA/ViennaRNA)，[Python API 文档](https://viennarna.readthedocs.io/en/latest/api_python.html)。
- de Smit & van Duin 1990/1994 对 E. coli 研究指出，起始区域结构稳定性与翻译效率相关，强结构会降低表达：[1994 分析](https://www.sciencedirect.com/science/article/pii/S0022283684717141)。
- Seo, Yang & Jung 2009 研究起始密码子下游区域结构对表达的影响，支持把早期 CDS 纳入起始区结构考量：[PubMed](https://pubmed.ncbi.nlm.nih.gov/19579224/)。
