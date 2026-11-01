# GSE132903 阿尔茨海默病 vs 对照：从 GEO 到火山图（赛博医学PhD 视频配套资料）

## 数据
- GEO 编号：GSE132903（Illumina HumanHT-12 V4.0，GPL10558）
- 样本：颞中回脑组织，阿尔茨海默病 (AD) 97 例 vs 非痴呆对照 (ND) 98 例
- 原始数据未打包，运行脚本会自动从 GEO 官方 FTP 下载到 input/

## 分析流程（scripts/01_volcano_pipeline.R）
1. 下载 series matrix 和平台注释 GPL10558
2. 探针转基因名，重复探针保留平均表达最高的（42,179 探针 → 28,859 基因）
3. 按 GEO2R 规则判断 log2：数据已在 log 尺度，不再转换
4. limma：AD − Control，BH 校正
5. 阈值 |log2FC| > 0.585（1.5 倍）且 adj.P < 0.05：上调 202 个，下调 279 个

注：本示例为教学演示，未校正年龄、性别等协变量。

## 运行
需要 R ≥ 4.3 和包 GEOquery、limma、ggplot2、ggrepel：
```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma")); install.packages(c("ggplot2","ggrepel"))
setwd("GSE132903_AD_volcano")   # 解压后的文件夹
source("scripts/01_volcano_pipeline.R")   # 下载、清洗、limma、默认火山图
source("scripts/02_volcano_academic.R")   # 美化版学术风格火山图
```

## 文件
- scripts/：01 为完整分析流程 + 默认火山图，02 为美化版学术风格火山图
- metadata/：样本分组表、探针-基因对照表
- output/：全部基因 limma 结果、显著差异基因表、清洗后的表达矩阵 (rds)、运行日志、sessionInfo
  - 注：`clean_expr_gene.rds`（清洗后的表达矩阵）未存放在本仓库中，由运行 `scripts/01_volcano_pipeline.R` 重新生成；`input/`（原始 GEO 下载）在脚本运行时自动创建。
- figures/：默认火山图 volcano_default (PDF/PNG)、美化版学术风格火山图 volcano_academic (PDF/PNG/TIFF 600 dpi)、样本箱线图 QC

数据来源：NCBI GEO GSE132903 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132903)，原始数据版权归数据提交者，本资料仅供学习。
