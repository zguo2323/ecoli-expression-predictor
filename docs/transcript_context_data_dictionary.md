# Kosuri 构建的转录上下文数据字典

## 来源与字段定义

项目使用 Kosuri et al. 2013 的 `sd01.xls`（Promoters）、`sd02.xls`（RBSs）和 `sd03.xls`（Constructs）。补充材料对序列和坐标的说明是：

- promoter `Sequence` 含开头的 cut site；序列中空格用于分隔设计片段，末端 5 nt 是用于 TSS 识别的 promoter barcode。
- `TSS.best` 是最常见的 TSS 相对 promoter/RBS junction 的位置；`TSS.pct_best` 是从该 TSS 开始的 contig 比例。
- RBS `Sequence` 含结尾 cut site。
- 该构建库以 pGERC 表达 superfolder GFP。Addgene 记录确认 pJ251-GERC（#47441）含 sfGFP reporter。

来源：[Kosuri et al. 补充材料](https://pmc.ncbi.nlm.nih.gov/articles/PMC3752251/)，[Addgene pJ251-GERC #47441](https://www.addgene.org/47441/)。

## 本分支构建的数据列

| 列名 | 含义 |
|---|---|
| `promo_seq` | `sd01.xls` 中的 promoter `Sequence`，去掉外层引号，保留原始内部空格与片段边界 |
| `rbs_seq` | `sd02.xls` 中的 RBS `Sequence`，去掉外层引号，保留原始内部空格与末端 cut site |
| `TSS_best` | 测量的优势 TSS 相对 promoter/RBS junction 的带符号偏移；负值表示在 junction 上游 |
| `TSS_pct_best` | 从优势 TSS 开始的 RNA contig 比例 |
| `transcript_prefix_to_start` | 从优势 TSS 起，到 RBS 末端 ATG 起始密码子（含 ATG）为止的 RNA 序列；DNA 的 T 被转为 U |
| `transcript_context` | `transcript_prefix_to_start` 加上 ATG 后的 87 nt reporter CDS 上下文；RNA 字母表 |

`transcript_prefix_to_start` 的拼接规则是：移除序列引号和用于展示的空白；将 `TSS_best` 转成相对 promoter/RBS junction 的序列索引；从拼接后的 promoter+RBS 构建序列该索引处开始截取，到末端 ATG 结束。RBS 表中的 111 条序列均以 `CATATG` 结束，因此末端 `ATG` 可作为当前数据的起始密码子边界；ATG 前的 `CAT` 保留为转录区序列。

## 早期 CDS 参考序列

仓库中的 [`sfgfp_first90nt.fasta`](../data/reference/sfgfp_first90nt.fasta) 保存了一条 90 nt 的 sfGFP reporter 序列。它来自 [Gilliot 博士论文附录 C.1.1](https://research-information.bris.ac.uk/files/432056523/Final_Copy_2024_01_23_Gilliot_P_PhD_Redacted.pdf) 对 Kosuri 构建的描述，包含起始 ATG 和后续 87 nt。由于 [Addgene pJ251-GERC 页面](https://www.addgene.org/47441/)上的质粒序列下载需要登录，目前无法将其与 Addgene 的 GenBank 记录逐碱基核对。故该序列标记为**暂定参考**，可用于方法开发和敏感性分析，不能表述为经质粒记录验证的精确 CDS。

## 当前边界与限制

- `transcript_context` 已用于提取 SD/start ensemble unpaired probability 与 SD 至 start 的 opening ΔG；CDS 上下文仅为上述暂定的前 90 nt。
- 当前 `compute_mrna_folding_energy` 保留为历史 promoter-tail MFE 计算函数，模型训练特征不再使用该值。它仍不是按实测 TSS 定位的生理转录本结构特征。
- ViennaRNA 的可及性使用 `pfl_fold_up` 的单序列 partition-function ensemble 概率，37°C、全转录上下文窗口、最大配对跨度 120 nt；opening ΔG 按 `−RT ln(Popen)` 换算。不模拟核糖体动力学、RNA 降解、离子/细胞条件或转录耦联。[ViennaRNA Python API](https://viennarna.readthedocs.io/en/latest/api_python.html)
- 这些测量使用单一 sfGFP reporter 与特定质粒/实验条件。即使补齐 CDS，推论也首先限于该 reporter 及相同构建语境。
