## =====================================================================
## 美化版火山图（学术期刊风格）
## 输入：output/GSE132903_AD_vs_Control_limma_all_genes.csv（由 01_volcano_pipeline.R 生成）
## 风格要点：NPG 配色、点大小随显著性变化、阈值虚线、上下调数量标注、
##          Top 基因带框标签、无网格的边框主题，适合直接放进论文 Figure
## 运行：工作目录设为资料包根目录后 source("scripts/02_volcano_academic.R")
## =====================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

deg <- read.csv("output/GSE132903_AD_vs_Control_limma_all_genes.csv")
FC_CUT <- 0.585; P_CUT <- 0.05                     # 与 01 脚本保持一致
deg$nlp <- -log10(deg$adj.P.Val)
deg$change <- factor(deg$change, levels = c("Up", "Down", "NS"))
nU <- sum(deg$change == "Up"); nD <- sum(deg$change == "Down")

## 需要标注的基因：每侧按 |logFC| × -log10(adj.P) 综合得分取前 8 个
sig <- subset(deg, change != "NS"); sig$score <- abs(sig$logFC) * sig$nlp
lab <- rbind(head(sig[sig$change == "Up", ][order(-sig$score[sig$change == "Up"]), ], 8),
             head(sig[sig$change == "Down", ][order(-sig$score[sig$change == "Down"]), ], 8))

cols <- c(Up = "#E64B35", Down = "#3C5488", NS = "#D9D9D9")   # Nature Publishing Group 配色
xlim <- ceiling(max(abs(deg$logFC)) * 10) / 10 + 0.1
ymax <- max(deg$nlp) * 1.05

p <- ggplot(deg, aes(logFC, nlp)) +
  # 先画不显著的点（灰色、小、半透明），再画显著点，避免被遮挡
  geom_point(data = subset(deg, change == "NS"), colour = cols["NS"], size = 0.8, alpha = 0.6) +
  geom_point(data = subset(deg, change != "NS"),
             aes(fill = change, size = nlp), shape = 21, colour = "white", stroke = 0.25, alpha = 0.9) +
  geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_label_repel(data = lab, aes(label = gene, colour = change),
                   size = 3, fontface = "italic", fill = "white", label.size = 0.25,
                   label.padding = unit(0.15, "lines"), box.padding = 0.5, max.overlaps = Inf,
                   segment.colour = "grey40", segment.size = 0.3, min.segment.length = 0,
                   show.legend = FALSE, seed = 42) +
  annotate("text", x =  xlim * 0.97, y = ymax, hjust = 1, vjust = 1, size = 4, fontface = "bold",
           colour = cols["Up"], label = paste0("Up: ", nU)) +
  annotate("text", x = -xlim * 0.97, y = ymax, hjust = 0, vjust = 1, size = 4, fontface = "bold",
           colour = cols["Down"], label = paste0("Down: ", nD)) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_size_continuous(range = c(0.8, 3.5), guide = "none") +   # 点越大越显著，与纵轴一致，不单列图例
  scale_x_continuous(limits = c(-xlim, xlim), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(log[2]~"fold change (AD / Control)"),
       y = expression(-log[10]~"(adjusted "*italic(P)*")"),
       title = "AD vs non-demented control",
       subtitle = "GSE132903, middle temporal gyrus (AD n = 97, Control n = 98)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey30", size = 10),
        legend.position = "none")

ggsave("figures/volcano_academic.pdf", p, width = 6.5, height = 6, device = cairo_pdf)
ggsave("figures/volcano_academic.png", p, width = 6.5, height = 6, dpi = 600)
ggsave("figures/volcano_academic.tiff", p, width = 6.5, height = 6, dpi = 600, compression = "lzw")
cat("Done: figures/volcano_academic.(pdf/png/tiff)\n")
