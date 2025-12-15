library(Seurat)
library(Seurat)
library(ggplot2)
library(ggpubr)


cr3_genes <- c("ITGB2", "ITGAM")
modules <- list(CR3 = cr3_genes)

celltypes <- unique(seurat_obj$predicted.celltype.l1)

for (ct in celltypes) {

  message("Procesando: ", ct)

  cells_ct <- subset(
    seurat_obj,
    subset = predicted.celltype.l1 == ct
  )

 
  if (length(unique(cells_ct$condition)) < 2) {
    message("  -> Se omite (solo un grupo)")
    next
  }


  cells_ct <- AddModuleScore(
    object = cells_ct,
    features = modules,
    name = "ModuleScore"
  )

  colnames(cells_ct@meta.data)[ncol(cells_ct@meta.data)] <- "CR3_Score"

  cells_ct$condition <- factor(
    cells_ct$condition,
    levels = c("HD", "AD")
  )

  df <- cells_ct@meta.data

  p <- ggplot(df, aes(x = condition, y = CR3_Score, fill = condition)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif"
    ) +
    theme_classic() +
    xlab("Condition") +
    ylab("CR3 Module Score") +
    scale_fill_brewer(palette = "Set2") +
    ggtitle(paste("CR3 score –", ct))

 
  ggsave(
    filename = paste0(
      "CR3_ModuleScore_HD_vs_AD_", gsub(" ", "_", ct), ".pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
}


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



