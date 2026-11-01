## =====================================================================
## GSE132903 阿尔茨海默病 (AD) vs 非痴呆对照 (ND) — 颞中回 差异分析 + 火山图
## Pipeline: GEO 下载 -> 探针注释 -> 数据清洗 -> limma -> 火山图
## Run:  Rscript scripts/volcano_pipeline.R   (工作目录 = 本文件夹根目录)
## =====================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(ggplot2); library(ggrepel)
})
options(timeout = 1200, stringsAsFactors = FALSE)
## 请把工作目录设为本文件夹根目录 (含 scripts/ 的那一层)，例如 setwd("~/GSE132903_AD_volcano")
for (d in c("input","metadata","scripts","output","figures")) dir.create(d, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

## ---- 1. 下载 series matrix ------------------------------------------
GSE <- "GSE132903"
msg("Step 1: downloading ", GSE, " series matrix")
sm <- file.path("input", paste0(GSE, "_series_matrix.txt.gz"))
if (!file.exists(sm))   # 直接从 GEO FTP 下载 series matrix
  download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/%s_series_matrix.txt.gz",
                        substr(GSE, 1, nchar(GSE) - 3), GSE, GSE), sm, mode = "wb")
gset <- getGEO(filename = sm, getGPL = FALSE)
msg("Title: ", gset@experimentData@title)
msg("Platform: ", annotation(gset), " | features x samples: ", paste(dim(gset), collapse = " x "))

## ---- 2. 平台注释 (GPL10558, Illumina HumanHT-12 V4.0) ----------------
msg("Step 2: downloading platform annotation ", annotation(gset))
gf <- file.path("input", paste0(annotation(gset), ".soft.gz"))
gpl <- if (file.exists(gf)) getGEO(filename = gf) else getGEO(annotation(gset), destdir = "input")
msg("GPL title: ", Meta(gpl)$title)
anno <- Table(gpl)[, c("ID", "Symbol")]
anno$Symbol <- trimws(anno$Symbol)
write.csv(anno, "metadata/GPL10558_probe2symbol.csv", row.names = FALSE)

## ---- 3. 分组 (来自 GEO pheno 'diagnosis:ch1') ------------------------
pd <- pData(gset)
pheno <- data.frame(gsm = pd$geo_accession, title = pd$title,
                    tissue = pd$`tissue:ch1`, diagnosis = pd$`diagnosis:ch1`,
                    sex = pd$`Sex:ch1`, age = pd$`expired_age (years):ch1`)
write.csv(pheno, "metadata/GSE132903_pheno.csv", row.names = FALSE)
group <- factor(ifelse(pheno$diagnosis == "AD", "AD", "Control"), levels = c("Control", "AD"))
msg("Step 3: tissue = ", paste(unique(pheno$tissue), collapse = ","))
msg("diagnosis table: ", paste(names(table(pheno$diagnosis)), table(pheno$diagnosis), sep = "=", collapse = ", "))
msg("Group sizes: AD = ", sum(group == "AD"), ", Control (ND) = ", sum(group == "Control"))

## ---- 4. log2 判断 (GEO2R 标准检查) -----------------------------------
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
msg("Step 4: quantiles = ", paste(round(qx, 3), collapse = ", "), " | needs log2: ", LogC)
if (LogC) { ex[ex <= 0] <- NaN; ex <- log2(ex); msg("log2 transform applied") } else
  msg("No log2 transform (data already VST/log2-scale, quantile normalized by submitters)")
png("figures/QC_boxplot_samples.png", width = 2400, height = 900, res = 150, type = "cairo")
boxplot(ex, col = ifelse(group == "AD", "#ff9a3c", "#2ee6e6"), outline = FALSE, las = 2,
        names = FALSE, main = "GSE132903 expression distribution per sample (orange=AD, cyan=Control)")
dev.off()

## ---- 5. 探针 -> 基因 (多探针取平均表达最高者) ------------------------
sym <- anno$Symbol[match(rownames(ex), anno$ID)]
keep <- !is.na(sym) & sym != ""
ex <- ex[keep, ]; sym <- sym[keep]
mm <- rowMeans(ex)
o <- order(sym, -mm)
ex <- ex[o, ]; sym <- sym[o]
first <- !duplicated(sym)
expr_gene <- ex[first, ]; rownames(expr_gene) <- sym[first]
msg("Step 5: probes with symbol = ", sum(keep), " -> unique genes = ", nrow(expr_gene),
    " (duplicates collapsed by max mean expression)")
saveRDS(list(expr = expr_gene, group = group, pheno = pheno), "output/clean_expr_gene.rds")

## ---- 6. limma 差异分析 ----------------------------------------------
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
fit <- lmFit(expr_gene, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(AD - Control, levels = design)))
deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg <- data.frame(gene = rownames(deg), deg, row.names = NULL)

## 阈值：|logFC| > 0.585 (1.5 倍) 且 adj.P < 0.05
## (AD 脑组织 bulk 芯片效应量小，|logFC|>1 的基因很少，故用 0.585)
FC_CUT <- 0.585; P_CUT <- 0.05
msg("For reference: genes with |logFC|>1 & adj.P<0.05 = ", sum(abs(deg$logFC) > 1 & deg$adj.P.Val < P_CUT),
    "; adj.P<0.05 any FC = ", sum(deg$adj.P.Val < P_CUT))
deg$change <- ifelse(deg$adj.P.Val < P_CUT & deg$logFC >  FC_CUT, "Up",
              ifelse(deg$adj.P.Val < P_CUT & deg$logFC < -FC_CUT, "Down", "NS"))
deg$change <- factor(deg$change, levels = c("Up", "Down", "NS"))
write.csv(deg, "output/GSE132903_AD_vs_Control_limma_all_genes.csv", row.names = FALSE)
write.csv(subset(deg, change != "NS"), "output/GSE132903_AD_vs_Control_DEGs_sig.csv", row.names = FALSE)
msg("Step 6: thresholds |logFC|>", FC_CUT, " & adj.P<", P_CUT,
    " -> Up = ", sum(deg$change == "Up"), ", Down = ", sum(deg$change == "Down"))
U <- deg[deg$change == "Up", ];   up <- head(U[order(-U$logFC), ], 5)
D <- deg[deg$change == "Down", ]; dn <- head(D[order(D$logFC), ], 5)
cat("\nTop5 Up (by logFC):\n");   print(up[, c("gene","logFC","adj.P.Val")], row.names = FALSE)
cat("\nTop5 Down (by logFC):\n"); print(dn[, c("gene","logFC","adj.P.Val")], row.names = FALSE)
cat("\nTop10 by adj.P:\n"); print(head(deg[, c("gene","logFC","adj.P.Val")], 10), row.names = FALSE)

## ---- 7. 火山图 -------------------------------------------------------
deg$nlp <- -log10(deg$adj.P.Val)
sig <- deg[deg$change != "NS", ]
sig$score <- abs(sig$logFC) * sig$nlp      # 标注：每侧取 |logFC| x -log10P 最大的 5 个
SU <- sig[sig$change == "Up", ]; SD <- sig[sig$change == "Down", ]
lab <- rbind(head(SU[order(-SU$score), ], 5), head(SD[order(-SD$score), ], 5))
nU <- sum(deg$change == "Up"); nD <- sum(deg$change == "Down")
xlim <- max(abs(deg$logFC)) * 1.08

## 默认风格火山图
cols_pub <- c(Up = "#D7301F", Down = "#2166AC", NS = "grey75")
p_pub <- ggplot(deg, aes(logFC, nlp, colour = change)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 3.2, fontface = "italic",
                  max.overlaps = 50, box.padding = 0.4, segment.colour = "grey50") +
  scale_colour_manual(values = cols_pub,
                      labels = c(Up = paste0("Up (", nU, ")"), Down = paste0("Down (", nD, ")"), NS = "NS")) +
  coord_cartesian(xlim = c(-xlim, xlim)) +
  labs(x = expression(log[2]~"Fold Change (AD vs Control)"), y = expression(-log[10]~"adjusted P"),
       colour = NULL, title = "GSE132903: AD vs non-demented control (middle temporal gyrus)",
       subtitle = paste0("AD n=", sum(group=="AD"), ", Control n=", sum(group=="Control"),
                         " | limma, |log2FC|>", FC_CUT, ", BH adj.P<", P_CUT)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12))
ggsave("figures/volcano_default.png", p_pub, width = 7, height = 6.5, dpi = 300, type = "cairo")
ggsave("figures/volcano_default.pdf", p_pub, width = 7, height = 6.5, device = cairo_pdf)

msg("Done. Figures in figures/, tables in output/")
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
