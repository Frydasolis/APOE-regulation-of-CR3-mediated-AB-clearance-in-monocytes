library(Seurat)
library(Seurat)
library(ggplot2)
library(ggpubr)
monocytes <- subset(seurat_obj, subset = predicted.celltype.l1 == "Mono")

cr3_genes <- c("ITGB2", "ITGAM")
modules <- list(CR3 = cr3_genes)

monocytes <- AddModuleScore(
  object = monocytes,
  features = modules,
  name = "ModuleScore"
)

colnames(monocytes@meta.data)[(ncol(monocytes@meta.data))] <- "CR3_Score"

monocytes$condition <- factor(monocytes$condition, levels = c("HD", "AD"))
monocytes$genotype <- factor(monocytes$genotype, levels = c("E3E3", "E3E4", "E4E4"))

table(monocytes$condition, monocytes$genotype)

df <- monocytes@meta.data

comparisons <- list(
  c("E3E3", "E3E4"),
  c("E3E3", "E4E4"),
  c("E3E4", "E4E4")
)

ggplot(df, aes(x = genotype, y = CR3_Score, fill = genotype)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif"
  ) +
  facet_wrap(~condition, scales = "free_x") +
  theme_classic() +
  xlab("Genotype") +
  ylab("CR3 Module Score") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


  ggsave(
  filename = "CR3_ModuleScore_by_APOE_condition.pdf",
  plot = last_plot(),
  width = 8, height = 5
)



