# 生物学建模优化阶段记录

## 阶段 1：重建测量 TSS 到起始密码子的转录前缀

**日期：** 2026-09-30
**分支：** `feat/biological-context-phase-1`

### 修改内容

- promoter 数据加载保留 `TSS.best` 和 `TSS.pct_best`；合并后的构建表也保留这两个字段。
- 新增 `build_transcript_prefix_to_start`：按 TSS 相对 promoter/RBS junction 的偏移拼接 promoter 与 RBS 序列，清理外层引号和展示空格，截取至末端 ATG，并将 DNA 的 T 转换为 RNA 的 U。
- 构建数据集时新增 `transcript_prefix_to_start`。该字段明确只到 ATG，后续 CDS 尚未拼接。
- 新增 [transcript_context_data_dictionary.md](transcript_context_data_dictionary.md)，说明源字段、坐标定义、序列边界及未解决限制。
- 管线集成测试改为写入 pytest 临时目录，避免测试覆盖仓库中的 processed parquet。

### 重点验证

1. 负偏移从 promoter/RBS junction 向 promoter 上游计数，零偏移从 junction 开始。
2. 返回序列必须从 TSS 对应索引开始、以末端 `AUG` 结束，且只含 RNA 字母表 `ACGUN`。
3. promoter/RBS 序列中的引号和展示空格不进入输出。
4. 空序列、非整数/越界 TSS、非 ATG 末端和非 DNA 字符要明确失败。
5. 原始表的每个可用构建都能生成非空前缀，TSS 元数据和输出 parquet 列完整保留。

### 测试设计

- 单元测试用人工序列验证 `tss_best=-2` 和 `tss_best=0` 两种坐标边界，并用参数化输入覆盖无效输入。
- 管线集成测试加载完整原始 xls 表，在临时目录生成 parquet，检查行数、必需列、TSS 与置信度完整性、前缀终止密码子和 RNA 字母表。
- 暂不测试 folding 数值，因为本阶段只构造序列上下文，未把它接入结构特征。

### 运行结果

- 阶段 1 定向测试：`pytest -q tests/test_features.py tests/test_pipeline.py -k 'transcript_prefix or build_dataset'` → **10 passed, 18 deselected**（1.53 s）。
- 相关测试文件整体：`pytest -q tests/test_features.py tests/test_pipeline.py` → **27 passed, 1 failed**（2.64 s）。唯一失败是既有的 `test_compute_mrna_folding_energy_context_dependent`：当前环境没有可导入的 ViennaRNA，旧实现对两条序列都静默返回 `0.0`，所以无法满足该 folding 行为断言。该失败与本阶段新增上下文构建无关，未通过修改断言或 folding 实现隐藏。
- `git diff --check`：通过。

### 本阶段结论与限制

TSS 到起始密码子的序列拼接已由定向单元测试和完整原始数据管线测试覆盖。`transcript_prefix_to_start` **不含 ATG 后 CDS**，所以现在还不能用于评价包含早期编码区的 RBS 可及性；必须先取得并核对 pJ251-GERC 中的 sfGFP reporter 序列。当前旧 folding 特征仍未接入新列，项目现有结构预测仍有原始实现的上下文缺陷。

## 阶段 2：基于 TSS 与早期 CDS 的翻译起始可及性

**日期：** 2026-09-30
**分支：** `feat/biological-context-phase-1`

### 修改内容

- 添加 `data/reference/sfgfp_first90nt.fasta`，使用 Gilliot 博士论文附录 C.1.1 报告的 sfGFP 起始 90 nt。因 Addgene pJ251-GERC 精确序列下载当前不可公开核对，FASTA header 与数据字典均将其标记为暂定来源，不宣称已与质粒 GenBank 逐碱基验证。
- `build_transcript_context` 将实测 TSS 到 RBS 末端 ATG 的 RNA 前缀与 ATG 后 87 nt CDS 拼接（起始 ATG 只保留一次），返回 0-based、半开区间的 SD 与 start 坐标。
- `build_dataset` 保留前缀并新增 `transcript_context`；模型特征改为使用 `TSS_best` 与 reporter CDS 计算 SD 平均未配对概率、start 平均未配对概率和 SD 至 start 整段 opening ΔG。GC、promoter motif、spacer、SD score/spacing 仍保留。
- 用 ViennaRNA partition function 计算同一结构 ensemble 内每个碱基的未配对概率，并以 SD 到 start 全区段硬约束/无约束 ensemble free-energy 差值计算 opening ΔG。固定 37°C、最大配对跨度 120 nt。ViennaRNA 缺失时 folding 函数显式报错，不再静默回传 0。
- 预测/特征 MCP 接受可选 `tss_best`；未提供时可及性特征标记为缺失，返回 `translation_context_available=false`。排名 MCP 同样可传目标 promoter 的 measured TSS。
- `load_and_featurize` 在同一 Python 进程内按输入 parquet 文件版本缓存特征矩阵，避免重复模型评估重新折叠相同数据。
- 更新 README 和数据字典，纠正旧 promoter-tail MFE 的生物学表述，记录 CDS 来源与 ViennaRNA 近似条件。

### 重点验证

1. 负 TSS 偏移、RBS SD motif 和末端 ATG 在含 CDS 的 RNA 上下文内坐标无错位，且 ATG 不重复。
2. SD/start 未配对概率都在 `[0,1]`，opening ΔG 非负；无 TSS 时明确标记结构信息缺失。
3. reporter 参考长度为 90 nt、以 ATG 开始且没有同框 stop codon；pipeline 中每条上下文都比前缀长 87 nt。
4. 未提供 ViennaRNA 时折叠功能显式失败；与环境无关的模型输入与序列特征仍可工作。
5. 新特征接入模型后预测、RBS 排名和评测入口仍能运行。

### 测试设计

- 特征单测以人工 promoter/RBS/CDS 构造序列，并直接断言 SD 与 AUG 区间；另用 ViennaRNA 计算可及性范围和 opening ΔG。
- pipeline 集成测试从仓库的原始 xls 构建完整表，在临时路径断言上下文字符、空值和长度关系。
- 模型集成测试固定随机抽取 600 条构建训练并检查预测/排名/评测 API；这样验证训练调用链，同时避免每轮接口测试都对全数据做数分钟的 RNA 折叠。
- 额外验证 benchmark 的 promoter、RBS、double-unseen 折分无身份泄漏且拆分可复现。

### 运行结果

- 初次模型全数据集成测试在 147 秒仍停留于 ViennaRNA 每构建重复运行多次 partition function；中断后 pytest 报告 **no tests ran**。此运行记为性能失败，不计入通过结果。
- 默认 pytest capture 在该虚拟环境即使只执行 GC 单测也以退出码 139 段错误（栈位于 pytest capture 的 `_readline_workaround`，发生在测试收集前）；使用 `-p no:capture` 后 pytest 正常运行。此环境问题保留记录，没有把捕获插件故障当成代码通过。
- 最初改用 `pfl_fold_up` 并做过一次全量 benchmark；数值核查发现部分长 SD-start 区域返回 0，与硬约束 partition function 不符。该轮 B6 结果已作废。
- 改为 base-pair probability 推导 SD/start 未配对概率，以及硬约束与无约束 partition function 的 free-energy 差，另加数值回归测试（指定构建 SD-start opening ΔG 在 1.3–1.6 kcal/mol）。模型测试数据集缩小到预先固定的 600 行 smoke 样本。
- 特征、pipeline 和 benchmark split 定向测试：`.venv/bin/python -m pytest -p no:capture -q tests/test_features.py tests/test_pipeline.py tests/test_biological_benchmark.py` → **34 passed, 3.77 s**。
- 模型集成测试：`.venv/bin/python -m pytest -p no:capture -q tests/test_model.py` → **8 passed, 15.58 s**。 smoke 样本 600 条，内部测试集 60 条；该轮初次结构实现的 Spearman 0.5418，仅为接口与拟合运行检查，不视作项目主结果。
- 最终全测试（含 partition-function 数值回归测试）：`.venv/bin/python -m pytest -p no:capture -q tests/` → **44 passed, 21.07 s**。 Smoke Spearman 0.5636；仍不作为项目效果结论。
- MCP server 导入/工具注册烟测 → 成功注册 12 个工具，包含新增可选 `tss_best` 参数的 feature/predict/rank 工具。
- `py_compile` 对 feature、pipeline、model、benchmark 和 MCP server 文件：通过。
- 全量模型特征预计算在约 11,696 条构建上耗时约 3–4 分钟；本轮直接计时 100 行约 1.94 秒。benchmark 会全量预计算一次，然后缓存当前进程中的矩阵供多个目标/拆分复用。

### 本阶段限制

- reporter 早期 CDS 仍为论文二手来源、未经 pGERC GenBank 核验；应取得可验证序列后做替换敏感性分析。
- ViennaRNA 的 base-pair/partition-function 概率基于固定参数的热力学 ensemble，不是细胞内实测结构。开放能量是模型转换出的特征，不能解读成实测核糖体解链能。
- 全量可及性计算约 2–4 分钟；用户可复用生成的 parquet，测试中对小样本验证模型集成。

## 阶段 3：有生物含义 baseline 与分组留出评测

**日期：** 2026-09-30
**状态：** 评测脚本、拆分测试与一次全量固定 holdout 基准已完成。

### 修改内容

- 新增 `scripts/biological_benchmark.py`，统一评估 `log(prot)`、`log(RNA)` 和 `log(prot/RNA)` steady-state proxy。
- 四种共用拆分：行级随机（诊断）、整组 promoter 留出、整组 RBS 留出、promoter 与 RBS 均未见的双重留出。double-unseen 的训练集同时排除测试 promoter 与测试 RBS，测试集只取两种 ID 都未见的组合。
- 共同比较训练集 median（B0）、只对已表征元件有效的 measured promoter-RNA×RBS-translation 组合 oracle（B2）、训练集 Ridge 元件身份加性模型（B1）、纯序列特征 XGBoost（B4）、原始 promoter-tail MFE XGBoost ablation（B3）与新 TSS accessibility XGBoost（B6）。所有模型复用同一拆分和目标值。
- 导出逐拆分/逐目标指标、每条 held-out 预测和 B6 特征重要性 CSV，便于组别误差复核。

### 重点验证与测试设计

- 人工 5×5 与 6×6 组合矩阵检查 group split 身份互斥、double-unseen 测试行非空、固定种子拆分相同。
- 全量运行后检查 CSV 行数、每场景的 train/test 规模和 B4→B6 的 Spearman/MAE 差值。结构特征收益只有在 held-out RBS 与 double-unseen 也改善时，才支持跨元件生物泛化。
- 比较总蛋白与翻译 proxy：若可及性只改善 translation proxy 而不稳定改善总 protein，应将结论限定为 translation-related 特征贡献。
- B2 是已测元件特性 oracle，仅作可组合性参考，不能作为未知 promoter/RBS 的部署预测能力。

### 运行结果

- split 单测包含在最终全仓 **44 passed** 中。
- 全量基准修正版运行：`.venv/bin/python -u scripts/biological_benchmark.py --data data/processed/constructs.parquet --output-dir reports/biological_context` → 成功，**11,696** 个构建、4 个 split × 3 targets × 6 个方法；结果已覆盖之前作废的输出。首次执行脚本暴露直接调用时根目录未加入 import path，之后已修复。
- `log(prot/RNA)` 翻译 proxy 上，B6 相比 B3 旧 promoter-tail MFE：random ρ **0.759→0.824**、MAE **0.622→0.542**；promoter holdout ρ **0.482→0.589**、MAE **0.868→0.758**；RBS holdout ρ **0.677→0.739**、MAE **0.682→0.612**；double-unseen ρ **0.410→0.496**、MAE **0.920→0.841**。
- `log(prot)` 总蛋白上，B6 的 Spearman 在四种拆分略有提高（double-unseen **0.552→0.555**），但 double-unseen MAE 从 **1.321** 增至 **1.354**。因此没有宣称所有指标都改善。
- 查看 [`benchmark summary`](../reports/biological_context/summary.md)、全指标 [`CSV`](../reports/biological_context/baseline_metrics.csv)、逐构建残差 [`CSV`](../reports/biological_context/heldout_predictions.csv) 和 B6 importance [`CSV`](../reports/biological_context/b6_feature_importance.csv)。
- 最终全仓测试：`.venv/bin/python -m pytest -p no:capture -q tests/` → **44 passed, 21.07 s**。源码编译检查与 `git diff --check`：通过。

### 评测限制

- `log(prot/RNA)` 是稳态蛋白/RNA proxy，不等同于瞬时翻译起始速率；当前数据 protein、RNA 均为正数，若未来出现非正数需预先制定检测限处理规则。
- B2 使用独立元件测量值，属于 characterization oracle；B1 对未见 ID 通过训练均值路径回退。
- 初版输出当前直接给出 XGBoost impurity importance 与逐行残差，不替代跨外层折的 group permutation importance、bootstrap 置信区间、RBS top-k/regret 或最终论文图。
- 这是单种子、单一 20% holdout，不含跨折置信区间，且没有把 OSTIR 加入同一 benchmark。虚拟环境安装了 OSTIR 1.1.3 并用一个序列验证 Python API 可返回 initiation-rate 预测；其输出伴随 ViennaRNA CLI 依赖缺失警告，因此没有将这次 smoke 结果并入正式 B7 baseline。翻译 proxy 上 B6 在四种场景均优于序列/MFE 对照，是初步生物学证据；总蛋白双新元件 MAE 仍稍差。后续需用 grouped CV、group bootstrap 和完整、版本锁定的 OSTIR baseline 检查稳定性。

## 阶段 4：OSTIR、分组交叉验证与 group bootstrap

**记录日期：** 2026-09-30

本阶段运行了 OSTIR 全量 baseline、5 折 grouped CV 和 1,000 次 group bootstrap，并完成 reporter CDS 来源审计尝试。具体方案见 [`biological_modeling_roadmap.md`](biological_modeling_roadmap.md) 的阶段 4 清单；CSV/JSON 输出见 [`reports/biological_context`](../reports/biological_context/summary.md)。

### 验证重点与测试设计

- 阶段目标：同一个 TSS-derived 5′ UTR + 暂定 sfGFP CDS 起始上下文用于 OSTIR；B3/B4/B6/B7 共用确定性外层折；以不泄漏测试折的训练集线性校准把 `log(OSTIR rate)` 映射到 `log(prot/RNA)`，原始排序由 Spearman 反映。
- 分别按行随机、promoter ID、RBS ID 生成 5 个外层折。double-unseen 的测试集是各折 held-out promoter 与 held-out RBS 的交集；训练集排除任何触及这两类 held-out ID 的行，其余交叉组合标为 unused。这样确保测试元件身份不出现在训练集。
- 外层每折为 `log(prot/RNA)` 与 `log(prot)` 重训 B3、B4、B6；B7 仅评估翻译 proxy。所有方法使用同一 test rows。调参不进行，固定现有 XGBoost 参数。
- bootstrap 对 promoter holdout 重采样 promoter 组，对 RBS holdout 重采样 RBS 组，对 double-unseen 使用 promoter × RBS 交叉 cluster multiplicity；随机行拆分只作为诊断，使用行级 bootstrap。记录模型间配对差值的百分位 95% CI。
- 新增定向测试覆盖合成 5×5 组合矩阵的组间无泄漏和 double-unseen 未使用边、fold 可复现性、OSTIR 正确定位 CDS AUG、OSTIR 缺失/非法输出校验，以及成对 bootstrap 区间。

### 环境与 OSTIR 执行结果

- `conda run -n ostir-b7` 中 `RNAfold`、`RNAsubopt`、`RNAeval` 都可找到；ViennaRNA Python/CLI 为 2.7.2，OSTIR 为 1.1.1。该 Conda 环境与项目 `.venv` 分离，因此 OSTIR 脚本通过 `conda run` 调用 CLI，并按 `benchmark_row_id` 写回预测。
- 单序列 smoke test 使用数据集中真实转录本，OSTIR 在报告 AUG 起始位置 34 返回 initiation-rate `93.0445`。10 行 CSV smoke test 成功。
- 全量运行输入 11,696 个转录本，输出 11,672 个可评分预测，覆盖率 **99.7948%**。24 个无 binding/initiation 记录的构建不伪装为零值，保存在 `ostir_unscored_rows.csv`；它们来自 BBa_J61124（3）、BBa_J61130（2）、BBa_J61131（19）。为让方法在同一测试集比较，所有方法都在 OSTIR 可评分的 11,672 个共同构建上评估。
- OSTIR 全量运行完成，但发出版本提示：ViennaRNA 2.7.2 高于 OSTIR 1.1.1 最后验证的 2.4.18。它不是 CLI 缺失错误；smoke 和全量均完成。报告保留当前版本并披露该兼容性风险，严格版本复现可后续另建 2.4.18 环境验证。

### 分组交叉验证和置信区间结果

- 翻译 proxy 五折 mean ± SD，B4 → B6 Spearman / log-MAE：promoter **0.530±0.073 → 0.643±0.076 / 0.865 → 0.752**；RBS **0.601±0.114 → 0.688±0.100 / 0.759 → 0.679**；double-unseen **0.331±0.118 → 0.461±0.071 / 0.982 → 0.879**。B6−B4 的 Spearman 与 MAE 配对 group-bootstrap 95% CI 在三个生物学 holdout 均不跨 0。`reports/biological_context/group_bootstrap_ci.csv` 保存各法区间和配对差值。
- B7 校准后的 OSTIR mean Spearman / MAE：promoter **0.461±0.078 / 0.879**、RBS **0.441±0.087 / 0.878**、double-unseen **0.433±0.115 / 0.868**。B7−B6 Spearman CI 在 promoter、RBS holdout 偏向 B6；double-unseen CI 跨 0，不能断言任一方更好。
- 总蛋白五折 mean Spearman / MAE，B4 → B6：random **0.853±0.006 → 0.896±0.007 / 0.730 → 0.640**；promoter **0.638±0.109 → 0.700±0.098 / 1.112 → 1.030**；RBS **0.735±0.072 → 0.784±0.061 / 0.935 → 0.863**；double-unseen **0.471±0.236 → 0.542±0.198 / 1.261 → 1.193**。四类中 B6−B4 的 group-bootstrap Spearman CI 均为正、MAE CI 均为负。
- 所有 1,000 次 bootstrap 均得到有效指标；fold seed 42、bootstrap seed 43。随机行诊断使用 row bootstrap，不作为主要生物学证据。

### CDS 来源审计、问题记录和剩余限制

- 查阅 Addgene pJ251-GERC sequence page 后确认该记录列有 1 条 depositor full sequence 和 3 条 Addgene-verified partial sequences，但页面要求登录才能读取/下载序列。故无法将当前从 Gilliot thesis 得来的 90 nt 与实际测量质粒做逐碱基比对；当前序列仍标为 provisional。证据、checksum 和关闭条件记录于 [`cds_source_audit.md`](../reports/biological_context/cds_source_audit.md)。
- 首次全量 OSTIR 运行后，原解析器因假定每行都必须有 OSTIR score 而失败；检查发现确实有 24 条 OSTIR 未评分。随后将解析逻辑改为显式允许缺失、保存缺失 row ID，并在所有方法共同完整案例上比较。此失败与输入或 ViennaRNA CLI 安装无关。
- 最初 OSTIR smoke 脚本不必要地调用全量特征提取，10 行测试也等待所有行 folding。已改为直接使用处理后 parquet 中保存的 `transcript_context`，避免 smoke 测试重复计算；group CV 生成的 `validation_dataset.parquet` 则保存完整模型特征供重跑复用。
- 这仍是同一来源数据集和单组预定 seed 的 5 折评估；bootstrap CI 描述组抽样不确定性，并不覆盖跨数据集、实验室、CDS 或宿主泛化。reporter CDS 权威核验和 sensitivity run 仍未解决；RBS top-k/regret 和跨折 group permutation importance 也未实现。

### 测试运行结果

- 新增测试设计：用 5×5 合成 promoter×RBS 矩阵验证 5 折无身份泄漏、double-unseen 交集和未用交叉行；复用 seed 检查折确定性；以 mock OSTIR CSV 验证 AUG 位置、缺失/非正值输出检查；用合成有组标签预测验证 paired group-bootstrap CI 生成。
- 定向测试：`.venv/bin/python -m pytest -p no:capture -q tests/test_grouped_biological_validation.py tests/test_biological_benchmark.py` → **7 passed**。之后解析器改为允许显式缺失输出，新增脚本仍由全套回归覆盖。
- 全量回归：`.venv/bin/python -m pytest -p no:capture -q tests/` → **49 passed, 20.60 s**。
- 最终全量回归（含 OSTIR 缺失输出允许列表的断言）：`.venv/bin/python -m pytest -p no:capture -q tests/` → **49 passed, 21.21 s**。
- `py_compile` 对 `scripts/biological_benchmark.py`、`scripts/run_ostir_baseline.py`、`scripts/grouped_biological_validation.py` 和 `scripts/plot_grouped_biological_validation.py` 通过；`git diff --check` 通过。
- 用 Conda base 中 Matplotlib 3.7.2 运行绘图脚本，生成 `reports/biological_context/grouped_validation_summary.png`；已目视检查四面板指标与误差区间。
- 产物完整性审计：5 折 × 4 场景的 assignment 可通过 `benchmark_row_id` 追溯；OOF 主键无重复，promoter/RBS/double-unseen 训练和测试身份无重叠；80 个 bootstrap 指标/配对差值条目均有 1,000 个有效 replicate。

### 备选后续方向

把同特征、同分组折、同 target 下的 XGBoost 与 Ridge/ElasticNet、Random Forest、ExtraTrees 比较，报告 Spearman、log-MAE 和折间波动；仅在该实验完成后再判断算法选择优势。该方向已记录在路线图“备选后续方向：验证算法选择”。

## 环境配置准备：Gemini `.env.example`

**日期：** 2026-09-30
**修改：** 新增 `.env.example`，列出代码实际读取的 `GEMINI_API_KEY`；README 改为先复制模板到 `.env` 再填入个人密钥。现有 `.gitignore` 忽略 `.env`，模板本身保持可跟踪且只含占位符。

**验证重点与结果：** 对照 `client_gemini.py` 的 `load_dotenv()` 与 Gemini 客户端默认密钥读取方式，环境变量名一致；人工确认模板不含真实密钥，`.env` 仍在忽略列表中。之后创建 `.venv` 并依 `requirements.txt` 安装项目依赖（含 ViennaRNA）；`.venv` 已被 `.gitignore` 忽略。

## 阶段 5：XGBoost 与 Ridge 同特征、同分组折算法比较

**日期：** 2026-09-30
**状态：** 实现、全量 5 折评估及配对 group bootstrap 已完成。

### 修改内容

- 在 `scripts/grouped_biological_validation.py` 中加入 Ridge，对 B3（历史 MFE）、B4（纯序列）、B6（TSS accessibility）各特征集和 `log(prot/RNA)`、`log(prot)` 两个目标分别训练。与 XGBoost 共用完全相同的外层拆分、训练行、测试行和测试目标。
- Ridge 使用 `StandardScaler + Ridge` pipeline。每个外层训练集内用 3 折 GridSearchCV 按 log-scale MAE 选 alpha；random 场景用行级折，promoter/RBS 场景按相应元件分组，double-unseen 内层验证要求 promoter 与 RBS 均未出现在该折训练集。测试折不参与预处理拟合或调参。
- 发现初始 alpha 网格有多次选中较大的 `1e3` 后，将 15 点网格扩展为 25 个 log-spaced 值 `[1e-5, 1e7]`，并完整重跑 grouped CV 和 1,000 次 group bootstrap。扩展后没有模型选择上端 `1e7`；部分 `log(prot)` 模型选择下端 `1e-5`，即近似无正则化。因此 alpha 下界选择已披露，若 Ridge 后续作为主要模型需增加 OLS/更低 alpha 的敏感性分析。
- 导出逐折和逐构建 Ridge 预测（含 `selected_alpha`）、标准化系数 `grouped_ridge_coefficients.csv`、同条件指标与 XGBoost−Ridge 配对 bootstrap 区间；总览图加入 Ridge B6。
- 路线图更新为：XGBoost/Ridge 算法对照已完成，ElasticNet、Random Forest、ExtraTrees 保留为候选扩展。

### 验证重点与测试设计

1. Ridge scaler 与 alpha 搜索只能接触各外层训练数据；内层拆分需遵循外层生物学分组规则，double-unseen 内层训练行不能含验证 promoter 或 RBS。
2. Ridge 每种特征集应与 XGBoost 在同一外层 fold、目标和测试行比较，并输出逐行预测和选择的 alpha。
3. bootstrap 对两个模型使用共同 held-out 组抽样，直接计算同一重采样样本的指标差，避免不配对区间比较。
4. alpha 网格扩展后检查是否仍有选择值卡在搜索端点；输出 Ridge 跨折均值/标准差、配对 95% CI 和模型对照图。

### 测试及运行结果

- 定向测试：`tests/test_grouped_biological_validation.py` 覆盖 promoter、RBS、double-unseen 内层拆分身份隔离，以及 Ridge 标准化、alpha 搜索范围、有限预测与系数。扩展 alpha 网格后的首次测试因测试断言仍写死旧网格而失败（**1 failed, 8 passed**）；改为导入实现中的 `RIDGE_ALPHAS` 常量后：`.venv/bin/python -m pytest -p no:capture -q tests/test_grouped_biological_validation.py` → **9 passed, 1.85 s**。
- 全量评估：`.venv/bin/python scripts/grouped_biological_validation.py --folds 5 --bootstrap-replicates 1000 --random-state 42` → 成功。XGBoost 每折参数保持既有设置；Ridge 的调参只发生在外层训练折。输出更新至 [`reports/biological_context`](../reports/biological_context/summary.md)。
- 翻译 proxy 的 B6 比较，XGBoost 对 Ridge 的五折 mean ± SD（Spearman / log-MAE）：random **0.827±0.007 / 0.542±0.009 vs 0.538±0.015 / 0.833±0.015**；promoter **0.643±0.076 / 0.752±0.051 vs 0.519±0.042 / 0.851±0.036**；RBS **0.688±0.100 / 0.679±0.048 vs 0.501±0.096 / 0.855±0.041**；double-unseen **0.461±0.071 / 0.879±0.063 vs 0.467±0.077 / 0.863±0.046**。
- 配对 group-bootstrap 中，XGBoost−Ridge 的翻译 proxy Δρ/ΔMAE 95% CI：promoter **[+0.071,+0.171] / [−0.142,−0.055]**，RBS **[+0.141,+0.245] / [−0.218,−0.138]**，double-unseen **[−0.090,+0.076] / [−0.047,+0.073]**。因此在 promoter 或 RBS 单组留出时支持 XGBoost 更好；双新元件时两者差值区间均跨 0，不能宣称有胜者。总蛋白双新元件差值区间同样跨 0。
- Ridge 标准化系数只作为线性模型描述，不解释为因果效应；特征相关性可能影响单个系数。
- 重绘 `grouped_validation_summary.png` 并目视检查四面板趋势、Ridge 图例与 95% bootstrap 区间：通过。
- 全仓测试：`.venv/bin/python -m pytest -p no:capture -q tests/` → **53 passed, 37.10 s**。其中既有 ViennaRNA folding 行为测试仍在本机依赖缺失时被 monkeypatch 跳过外部计算，相关限制沿用阶段 1 记录。
- `py_compile`（grouped validation、plotter、定向测试）及 `git diff --check`：通过。
- 产物完整性：260 条 fold×target×method 指标；485,589 条 OOF 预测；1,240 条 Ridge 系数/alpha 记录；176 个 bootstrap 指标/差值全部各有 **1,000** 个有效 replicate。alpha 实际选择范围 `1e-5` 到 `1e3`，未触及扩展后上界 `1e7`。

## 教授评语回应状态记录

**日期：** 2026-10-06

已将教授关于转录起始区/CDS 上下文、baseline、feature importance、error analysis 和按 promoter/RBS 整组留出的意见逐条映射到当前代码与评测证据。现状、指标和仍未闭环的限制见 [`professor_feedback_status.md`](professor_feedback_status.md)。

用户已说明毕业后无法满足 Addgene 注册要求的组织信息条件，并决定将 pJ251-GERC 序列访问作为客观限制跳过。当前 90 nt 继续标记为论文来源的 provisional CDS，不能称为已核验实验质粒序列。held-out 分组置换重要性已完成，下一项优先补强工作改为对 OOF 残差做系统误差分析；详情见 [`professor_feedback_status.md`](professor_feedback_status.md) 与 [`cds_source_audit.md`](../reports/biological_context/cds_source_audit.md)。

## 阶段 6：B6 生物学特征组的 held-out permutation importance

**日期：** 2026-10-06
**状态：** 代码、测试、全量置换分析和图表完成。

### 设计与修改

- 新增 `scripts/grouped_feature_group_permutation.py`，在既有 4 种场景 × 5 折 × 2 个目标中重训 B6 XGBoost。每一折只在该折 held-out rows 上置换，并比较置换前后 Spearman 和 log-scale MAE。
- 将 B6 的 12 个输入分为三组：`promoter_recognition`（promoter GC、−10/−35 motif、spacer）；`rbs_recognition`（RBS GC、SD score 与 spacing）；`transcript_accessibility`（SD/start 未配对概率和 SD-start opening energy）。每个折、目标和特征组做 30 次置换。
- 先运行 row-wise 探索版后，检查数据发现 promoter 特征在 promoter ID 内恒定、RBS 特征在 RBS ID 内恒定、accessibility 特征在 promoter×RBS pair 内恒定。row-wise 置换会制造同一生物实体有多个相互冲突特征值的输入，故作废该探索版输出。正式方案改为把特征向量分别在 promoter、RBS、promoter×RBS 单位间置换，单位内所有重复行使用同一组置换值，并在运行时检查单位内不变量。
- 正值 delta 定义为原始模型表现减去/优于置换后表现：`delta_spearman = rho_original - rho_permuted`；`delta_mae_log = MAE_permuted - MAE_original`。每折先取 30 次置换均值，再等权汇总五折 mean ± SD；这表示折间变异，不是置信区间或显著性检验。
- 输出 repeat-level、fold-level、aggregate CSV 与四面板图；README summary 与教授评语状态同步更新。

### 重点验证与测试设计

1. 三个特征组必须无重复且完整覆盖 B6 的 12 个特征。
2. 对 promoter/RBS/pair 单位置换时，特征必须在实体内恒定；同组特征使用同一个来源实体，保留组内相关性；实体重复行收到一致的新值。
3. 合成实测/预测数据验证重要性符号：打乱后 Spearman 下降、MAE 上升时两个 delta 应为正。
4. 汇总函数需先在折内平均随机置换，再对外层折等权计算均值和 SD。
5. 全量结果检查四场景、两目标、三组、五折、每组 30 次均存在，预测行覆盖与既有 OSTIR complete-case 外折 assignments 一致。

### 运行结果与解释

- 定向测试：`.venv/bin/python -m pytest -p no:capture -q tests/test_grouped_feature_group_permutation.py` → **5 passed, 1.57 s**。
- 正式全量分析：`.venv/bin/python scripts/grouped_feature_group_permutation.py --permutations 30 --random-state 42` → 成功；3,600 条 repeat-level 记录（4 scenarios × 5 folds × 2 targets × 3 feature groups × 30 permutations），fold-summary 120 行、aggregate-summary 24 行。每个 biological holdout 的 OOF 测试覆盖 11,672 个 complete-case 构建；double-unseen 每折交集共 2,337 行。
- `log(prot/RNA)` 上，各生物学场景中三组都在每个折均给出正向平均 importance。promoter holdout 的 Spearman drop / MAE increase：promoter **0.144±0.102 / 0.110±0.077**、RBS/SD **0.346±0.066 / 0.277±0.041**、accessibility **0.185±0.062 / 0.173±0.038**。RBS holdout 对应 **0.347±0.027 / 0.278±0.031**、**0.225±0.106 / 0.186±0.097**、**0.151±0.044 / 0.158±0.052**。double-unseen 对应 **0.168±0.093 / 0.113±0.070**、**0.168±0.094 / 0.119±0.060**、**0.123±0.023 / 0.109±0.025**。
- 解释：accessibility 特征组在 promoter、RBS、双新元件留出中均显示额外模型依赖，但平均影响低于当前占优的序列组或与其相近。重要性相对顺序随留出场景变化，promoter 和 double-unseen 的部分 feature groups 折间 SD 较大。该结果支持结构信息有补充预测价值，不能证明 RNA 结构的因果效应。
- 结构化误差分析成为下一项补强工作；pJ251-GERC 权威序列核验按用户要求标记为因 Addgene 注册条件客观限制跳过，当前 reporter CDS 仍为 provisional。

### 失败与最终检查

- 第一版完整分析曾采用 row-wise permutation，已因不符合 promoter/RBS/pair 内特征恒定结构而明确作废，并由 biological-unit 版覆盖。
- biological-unit 版第一次执行在预测前尝试把含序列字符串与 ID 的完整 dataframe 转成浮点数，运行失败；修正为只传入 B6 数值特征后全量运行成功。
- `.venv/bin/python -m pytest -p no:capture -q tests/` → **58 passed, 20.41 s**；`py_compile` 与 `git diff --check` 通过。
- `conda run -n base python scripts/plot_grouped_feature_group_permutation.py` 生成图表，并已目视检查 4 面板、三组图例和折间误差线。

## 阶段 7：Grouped OOF 残差与元件分层误差分析

**日期：** 2026-10-06
**状态：** 代码、定向测试、OOF 诊断表和图表完成；全仓测试结果记录于本节末。

### 问题与设计

- 重点追查阶段 6 发现的现象：RBS holdout 中 promoter-recognition 特征组的 permutation effect 高于 RBS/SD 组。分析同一批 OOF 行的 B6、同特征 Ridge 和 OSTIR 预测，检查 B6 相对 Ridge 的误差优势是否只集中在少数 promoter/RBS。
- 对每个 scenario/target/model，按实测 target 构造共享的五分位箱，计算样本数、signed residual（actual−predicted）、MAE、median absolute error 和 RMSE。分箱仅用于事后诊断，不参与模型训练、调参或校准。
- 按 held-out promoter ID、RBS ID、以及 double-unseen 的 promoter×RBS pair 汇总 residual metrics；同一 OOF 行配对比较 B6 与 Ridge、B6 与 OSTIR 的绝对误差，报告每个元件组的 MAE 差和组间分布。
- 新增目标值的单因素边际 eta-squared，量化 target 按 promoter 或 RBS ID 分组均值所关联的总平方和比例，专门用于检验“promoter 特征重要”是否可简单归因于 promoter 组间目标变异更大。该量不是交叉设计的方差分解，不能解析交互或因果机制。
- 导出极端残差行，保留 construct、序列身份、fold、真实值、预测值与误差方向，便于后续人工检查。绘制 RBS holdout 的组件级模型配对误差与 translation-proxy 分位残差图。

### 重点验证与测试设计

1. 合成数据验证同一 benchmark row 在不同模型之间获得相同 target 分位；分位计数完整，且改变分位数参数可测试。
2. 合成 double-unseen 数据验证 promoter、RBS 及 promoter×RBS 三种 component error 分组，配对比较仅对相同 row 求误差差；正的 Ridge−B6 MAE 差表示 B6 更好。
3. 可控的二因子合成目标验证 eta-squared 按 ID 聚合：完全由 promoter 决定时 promoter eta²=1、RBS eta²=0。
4. 极端案例测试验证每种 scenario/target/model 限制输出数，并正确识别 underprediction/overprediction。
5. 全量产物核对 scenario、target、model 覆盖与 CSV 行数；查看 RBS-holdout 和 double-unseen 主要结论，目视检查图表标注、图例、误差方向和组件 scatter 对角线。

### 运行结果与解释

- 定向测试 `.venv/bin/python -m pytest -p no:capture -q tests/test_grouped_oof_error_analysis.py` → **6 passed, 0.48 s**。早期测试迭代中，分箱测试初始假设与默认五分位数量不符，另一个配对误差测试使用精确浮点相等导致失败；改用显式分箱数与浮点容差断言后通过。之后新增 eta-squared 合成测试，并将分析输出契约测试扩展为检查 `target_component_variance`。
- 复现命令 `.venv/bin/python scripts/grouped_oof_error_analysis.py --extreme-cases 20` 成功。输出：actual-bin 100 行、component error 16,145 行、component summary 45 行、target component variance 16 行、paired component comparison 9,687 行、paired summary 27 行、极端残差 400 行（每 scenario×target×主模型各 20 条）。
- RBS holdout、`log_translation_proxy` 中，B6 对 Ridge 的 component-level MAE 更低：已知 promoter 背景 **108/112（96.4%）**，held-out RBS **91/111（82.0%）**。跨组 median MAE：promoter 组 **0.679 vs 0.807**、RBS 组 **0.584 vs 0.766**（B6 vs Ridge）。提示优势较广泛，但不能解释 promoter 组的 permutation effect 为什么最大。
- 继续检查其它严格场景：promoter holdout 下 B6 组 MAE 胜过 Ridge 的 promoter 比例 **72.3%**、RBS 比例 **67.6%**；double-unseen 下 pair 胜组率 **47.5%**（1,110/2,337 pairs）。RBS-holdout 的广泛优势不能推广为模型在所有生物留出场景都一致占优。
- 对应 target marginal eta-squared：RBS holdout promoter **0.346**、RBS **0.385**；double-unseen promoter **0.370**、RBS **0.386**。因此目标的单因素组间变异没有显示 promoter 明显大于 RBS。该指标忽略 promoter×RBS 交互和交叉结构，不能将其与 permutation effect 视作同一种量。
- 翻译 proxy residual 有明显范围压缩：RBS holdout B6 的最低实测五分位平均 residual **−0.948**（过预测，MAE **1.078**, n=2,335），最高五分位 **+0.713**（低预测，MAE **0.782**, n=2,335）；double-unseen 最低/最高五分位分别为 **−1.209 / +1.051**（MAE **1.340 / 1.112**, n=468）。Ridge 与 OSTIR 同方向，因此现有证据不能把该偏差归为某一个算法或某项生物机制。
- held-out salis-3-3 仍是局部失败案例：77 个构建的 B6 MAE **1.096**，Ridge **0.679**；salis-3-11 为 **1.675 vs 1.423**（112 构建）。需要配合原始测量与序列背景进一步核验，当前不赋予机制解释。
- 图表由 `conda run -n base python scripts/plot_grouped_oof_error_analysis.py` 生成并目视检查；标注仅保留一个最大模型差异 RBS，防止多标签重叠。
- `.venv/bin/python -m pytest -p no:capture -q tests/` → **64 passed, 21.12 s**；包括本阶段 6 个定向测试与此前已有测试。
- `.venv/bin/python -m py_compile scripts/grouped_oof_error_analysis.py scripts/plot_grouped_oof_error_analysis.py tests/test_grouped_oof_error_analysis.py` 与 `git diff --check` → 通过。

### 限制与后续

- 该步骤是 OOF 预测误差的描述性分层，没有对误差类别（低计数测量噪声、序列边界、预测结构不符等）做自动或人工的生物成因归类。
- actual-quintile 的误差趋势不能作为校准后效果；不应用这些 OOF test labels 来修改模型或做事后校正。
- 当前发现支持更透明地说明模型适用范围和失败组，但没有证明结构可及性或 promoter feature 的因果作用。优先后续为重复外层分组/多个随机种子、或可核实独立数据验证；Addgene 序列访问仍按用户决定跳过。

## 阶段 8：多随机种子 grouped validation 稳健性检查

**日期：** 2026-10-06
**状态：** 4 个外层 split seeds（42、7、123、2026）已完成并汇总；最终全仓测试结果见本节末。

### 问题与设计

- 目标是看主要 grouped-CV 结论是否依赖 seed 42 刚好抽到的那一组外层折，而不是新增一轮难以解释的机制归因。
- 对 seeds **7、123、2026** 重跑现有 5-fold grouped CV，与既有 seed **42** 逐折结果并列。每个 seed 都使用相同的 OSTIR complete-case 构建（11,672 行）、四种外层场景（random、promoter holdout、RBS holdout、double-unseen）、两个目标（`log_translation_proxy`、`log_prot`），及相同 B3/B4/B6、Ridge、OSTIR 实现与指标。OSTIR 未评分的 24 个行仍一致排除。
- 改变的是外层分组和 Ridge 内层 CV 的随机划分；XGBoost 本身 `random_state` 保持固定 42。因此本实验主要测量**拆分敏感性**，不涵盖 XGBoost 初始化/训练随机性。
- 新增 `--skip-bootstrap` 参数以避免每个 seed 重复计算已有 seed42 的 1,000 次组 bootstrap。seed 42 的原 bootstrap/配对区间保持作为固定切分的不确定性结果；多 seed 结果另外报告每个 seed 的五折均值、seed 间 SD/range 和 B6−B4、B6−Ridge、B6−OSTIR 的配对折差方向。
- 新增 `scripts/summarize_multiseed_grouped_validation.py` 和合成数据测试，生成拼接 fold metrics、seed-level metrics、跨 seed 汇总及配对差异 CSV。

### 重点验证与测试设计

1. 每个 seed 必须保存 4 场景×5 折×2 目标的相同模型结果；确认预测覆盖和 complete-case 样本数一致。
2. 汇总测试用小型合成数据检查 seed 均值、跨 seed 统计、配对比较和指标符号（Spearman 增大有利，MAE 减小有利）。
3. 逐 seed 检查主要 biological holdout 的 B6−B4、B6−Ridge 配对差方向；分别展示严格 double-unseen，不以较容易的 random split 掩盖它。
4. 元数据明确记录 multi-seed run 跳过 bootstrap，且原 seed42 bootstrap 文件未被覆盖。

### 运行过程、失败记录与结果

- 首次 seed7 命令沿用默认 1,000 次 bootstrap。分组 OOF 预测、fold metrics 等文件已经写完，但重复 bootstrap 运行耗时较长；为避免每个新 seed 都重做已存在的固定 seed bootstrap，主动中断该冗余步骤。中断栈位于 `_resample_group_rows`，是 `KeyboardInterrupt`，没有发现模型或分折异常；保留其 OOF 文件后，以新增 `--skip-bootstrap` 正式重跑并完整写入 metadata。该中断不计作有效评估结果。
- 正式命令分别为：`.venv/bin/python scripts/grouped_biological_validation.py --folds 5 --skip-bootstrap --random-state 7 --output-dir reports/biological_context/multiseed/seed_7`，seed 123 和 2026 同理替换参数/目录。三个 seed 均成功，metadata 的 `bootstrap_skipped=true`，complete-case `n_common_scored=11672`。
- 汇总命令：`.venv/bin/python scripts/summarize_multiseed_grouped_validation.py --run 42=reports/biological_context --run 7=reports/biological_context/multiseed/seed_7 --run 123=reports/biological_context/multiseed/seed_123 --run 2026=reports/biological_context/multiseed/seed_2026`。输出 1,040 条 fold metrics、208 条 seed-level summaries、52 条 model summaries、800 条 paired fold deltas、160 条 paired seed deltas 和 40 条 paired-delta summaries，均位于 `reports/biological_context/multiseed/`。
- `log_translation_proxy` 的 B6 4-seed mean ± SD across seed-level 5-fold means（Spearman / log-MAE）：promoter holdout **0.630±0.009 / 0.759±0.008**；RBS holdout **0.705±0.017 / 0.672±0.006**；double-unseen **0.476±0.021 / 0.877±0.013**。
- B6 相对 B4 sequence-only 在 promoter、RBS、double-unseen 三种生物学 holdout 中，4/4 seeds 的平均 CV Spearman 更高且 MAE 更低。B6 相对同特征 Ridge，在 promoter 和 RBS holdout 中也都 4/4 seeds 指标更好；double-unseen 仍没有稳定胜者（B6−Ridge Spearman 差仅 2/4 seeds 为正，MAE 差仅 1/4 seeds 有利 B6）。这与原有 seed42 的双新元件结果一致：核心结构特征模型在单组未见时稳健优于 sequence-only / Ridge，但双新元件外推依然困难，不能靠 seed 变化消除。
- `test_summarize_multiseed_grouped_validation.py` 与 `test_grouped_biological_validation.py` 定向测试 → **10 passed, 1.65 s**。最终 `.venv/bin/python -m pytest -p no:capture -q tests/` → **65 passed, 22.79 s**。新增/修改脚本 `py_compile` 与 `git diff --check` → 通过。

### 产物位置与范围

- 每 seed 的独立目录：`reports/biological_context/multiseed/seed_{7,123,2026}/`，含 OOF predictions、fold assignments、fold metrics、B6 feature importance、Ridge coefficients 和 run metadata。seed42 原始产物仍在 `reports/biological_context/`。
- 汇总表：`multiseed_cross_seed_summary.csv`（主汇总）、`multiseed_paired_delta_summary.csv`（同折配对模型差）、`multiseed_seed_level_metrics.csv`（各 seed 明细），以及其他完整长表。
- 多 seed 稳健性是相同数据、相同设计下的重分组检查，不是独立数据集/新实验复现。小 seed 数量的 SD 与胜出比例用于描述稳定性，不作为正式显著性检验或总体泛化保证。

## 阶段 9：移除未校准的预测区间，并确定外部数据候选

**日期：** 2026-10-06
**状态：** 固定 ±25% 的伪置信区间已从 API 移除；独立公开数据源已检索并评估适配性。

### 修改与验证

- 原 `predict_expression` 返回 `confidence_interval=[prediction×0.75, prediction×1.25]`。代码和 [`SPEC.md`](../SPEC.md) 已明确这只是 placeholder，没有由 held-out residual、bootstrap 或 conformal calibration 推导，因此不能作为统计置信区间。现在接口仅返回 `predicted_prot`、`features_used`、`translation_context_available`；RBS 排名输出原本只传递点预测，没有区间字段。
- 同步更新 `modules/model/predict_expression.json` 的 outputs/example、`modules/model/SKILL.md`、`tests/test_model.py`。测试现在要求点预测及特征字段存在，并明确断言 API 不再输出 `confidence_interval`。
- 定向测试 `.venv/bin/python -m pytest -p no:capture -q tests/test_model.py` → **7 passed, 18.96 s**。代码/JSON/SPEC/Skill 搜索确认旧 placeholder 只留有测试“字段不存在”的断言；另外的 grouped bootstrap 置信区间属于有效、独立计算的评测统计量，并未移除。
- 修改后的全仓测试 `.venv/bin/python -m pytest -p no:capture -q tests/` → **64 passed, 21.34 s**；`git diff --check` 和 `modules/model/model.py` 的 `py_compile` 均通过。外部数据只完成检索和适用性审查，尚未导入或计算外部验证分数。

### 外部测量数据候选

公开文献检索发现可作为后续独立验证的真实测量数据：

1. **Wen et al., 2024（优先评估）**：在 *E. coli* 中测试 15 个人工 regulatory sequences、41 个代谢基因，报告 520 组调控序列×基因的 GFP fusion fluorescence；Supplementary Table S3 有荧光读数，mmc1.xlsx/mmc2.docx 为开放补充材料。方法包含 promoter、RBS、目标基因起始 180 bp 与 GFP fusion context。宿主/培养条件、载体、reporter/CDS context 不同；并且所选 regulatory sequences 来自既有 promoter-RBS library，可能与 Kosuri 训练集中的 regulatory identities 重叠。必须先做精确 sequence overlap audit；若有重叠，称为独立实验条件/编码上下文验证，不称为新 promoter/RBS 泛化。不能把原始荧光数值与 Kosuri `prot` 直接校准比较，但可在确定序列边界后，以冻结模型做外部 rank/相关性分析，并按 gene/CDS 分层检查排序变化。论文/数据入口：[ScienceDirect article](https://www.sciencedirect.com/science/article/pii/S2405805X24000851)，[PMC full text and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC11137365/)。
2. **Dash et al., 2024 Dryad**：16 个 E. coli global-regulator full-length promoter::GFP reporters，Dryad 提供 plasmid FASTA/图谱以及单细胞流式和 bulk fluorescence 文件，荧光单位为 a.u.。可验证固定 GFP/RBS 框架下的 promoter 排序或相关性；没有 RBS 变化，压力/生长条件也与训练集不同：[Dryad dataset](https://datadryad.org/dataset/doi:10.5061/dryad.b2rbnzsm8)。
3. **Bonde et al., 2016 EMOPEC**：补充材料有 3,864 个 GFP expression estimate、106 个经单克隆验证的样本及序列/寡核苷酸；更适合作为 RBS 设计比较的次级候选。因关键表包含 estimated/desired-level quantities，正式采用前要区分模型估计值与真实实验测量，不能把全部 3,864 行都当独立实测标签：[Nature Methods article and supplements](https://www.nature.com/articles/nmeth.3727)，[PubMed record](https://pubmed.ncbi.nlm.nih.gov/26752768/)。

推荐先核对第一项补充材料中能否逐一恢复 promoter/RBS/leader/CDS 边界和实际测量列；若边界明确，用它作为独立数据集，不参与重训或校准，优先报告 rank Spearman/Kendall 与样本数，并分 gene/CDS 报告。若无法可靠拆分 regulatory sequence，则采用 Dryad promoter reporter 数据做范围较窄的 promoter-only 检查。外部数据 benchmark **尚未执行**，不能将文献中他人报告的结果说成本项目验证结果。

## 阶段 10：收敛对话 MCP 工具并更新 README

**日期：** 2026-10-07

### 审查与修改

- 逐一核对 MCP JSON 定义、server 注册逻辑、Python 实现与模块 Skill。保留 7 个序列解释工具和 2 个面向设计的模型工具：`extract_all_features`、单特征计算工具、`predict_expression`、`rank_rbs_for_promoter`。
- 从 MCP 注册面移除 3 个不适合会话调用的入口：旧 promoter-tail MFE（会被误读为真实转录本结构）、`evaluate_model`（依赖本地 parquet 与随机拆分，不能替代严格分组评估）、`build_dataset`（离线数据准备）。离线 Python 实现/评测流程保留。
- MCP 返回特征包含 NaN 时转换为标准 JSON `null`，避免未知 TSS 时结构特征以非标准 JSON NaN 返回。删除预测/特征 schema 中未经真实运行核验的数值样例；同步更新模块 Skill、README 及历史 SPEC 的现行状态提示。
- README 补齐工具目录、离线构建与训练、MCP 客户端配置、实验前未知 TSS 的用法边界、预测尺度/置信区间限制，以及生物验证报告入口。

### 测试设计与结果

- 新增 `tests/test_mcp_tools.py`：用注册清单断言 MCP 只暴露 9 个预期用户工具；用含 NaN/Inf 的嵌套结果断言返回值可解析为严格 JSON，非有限数转成 null。
- 初始 `unittest discover` 没有收集到测试（本仓库测试使用 pytest 风格）；改用项目 pytest 命令运行完整套件：`.venv/bin/python -m pytest -p no:capture -q tests/` → **66 passed, 20.97 s**。
- `git diff --check`、模块 JSON 解析和 `server.py` 编译检查 → 通过。

## 阶段 11：缺失实测 TSS 时先询问用户

**日期：** 2026-10-07

### 修改与验证重点

- 在 `predict_expression`、`rank_rbs_for_promoter` 与 `extract_all_features` 的 MCP 描述中加入交互规则：缺少 TSS 时先问用户是否有匹配构建语境的实测 `tss_best`；提供则使用；没有、不确定或希望继续则省略并继续；不猜测，也不重复询问同一构建语境。
- 同步更新 features/model Skill 和 README。追问由 MCP 客户端中的助手发起；MCP 工具本身仍是接收可选参数的无状态函数。
- 新增回归测试检查三个 schema 描述都要求“ask once”且禁止推断/编造 TSS。

### 测试结果

- 完整套件：`.venv/bin/python -m pytest -p no:capture -q tests/` → **67 passed, 22.96 s**。
- JSON 解析及 `git diff --check` → 通过。
