
# ─────────────────────────────────────────────

library(Seurat)
library(dplyr)
library(ggplot2)
mono_obj <- subset(alzheimerfiltercomplete, predicted.celltype.l1 == "Mono")
Idents(mono_obj) <- "condition"

deg_ADvsHD <- FindMarkers( 
  mono_obj,
  ident.1 = "AD",
  ident.2 = "HD",
  group.by = "condition",
  test.use = "MAST",
  logfc.threshold = 0
)
deg_ADvsHD$gene <- rownames(deg_ADvsHD)


# DEG 2: APOE3/4 vs APOE3/3
Idents(mono_obj) <- "APOE"
deg_APOE34vs33 <- FindMarkers(
  subset(mono_obj, subset = APOE %in% c("E3E3", "E3E4")),
  ident.1 = "E3E4",
  ident.2 = "E3E3",
  group.by = "APOE",
  test.use = "MAST",
  logfc.threshold = 0
)
deg_APOE34vs33$gene <- rownames(deg_APOE34vs33)


DEG 3: APOE4/4 vs APOE3/3

deg_APOE44vs33 <- FindMarkers(
  subset(mono_obj, subset = APOE %in% c("E4E4", "E3E3")),
  ident.1 = "E4E4",
  ident.2 = "E3E3",
  group.by = "APOE",
  test.use = "MAST",
  logfc.threshold = 0
)
deg_APOE44vs33$gene <- rownames(deg_APOE44vs33)

deg_APOE44vs43 <- FindMarkers(
  subset(mono_obj, subset = APOE %in% c("E4E4", "E3E4")),
  ident.1 = "E4E4",
  ident.2 = "E3E4",
  group.by = "APOE",
  test.use = "MAST",
  logfc.threshold = 0
)
deg_APOE44vs43$gene <- rownames(deg_APOE44vs43)


merged_deg3433 <- inner_join(
  deg_ADvsHD[, c("gene", "avg_log2FC", "p_val_adj")],
  deg_APOE34vs33[, c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene",
  suffix = c("_ADvsHD", "_APOE")
)

merged_deg4433 <- inner_join(
  deg_ADvsHD[, c("gene", "avg_log2FC", "p_val_adj")],
  deg_APOE34vs33[, c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene",
  suffix = c("_ADvsHD", "_APOE")
)
merged_deg4434 <- inner_join(
  deg_ADvsHD[, c("gene", "avg_log2FC", "p_val_adj")],
  deg_APOE44vs43[, c("gene", "avg_log2FC", "p_val_adj")],
  by = "gene",
  suffix = c("_ADvsHD", "_APOE")
)

merged_deg3433 <- merged_deg3433 %>%
  filter(p_val_adj_ADvsHD < 0.05 | p_val_adj_APOE < 0.05)
merged_deg4433 <- merged_deg4433 %>%
  filter(p_val_adj_ADvsHD < 0.05 | p_val_adj_APOE < 0.05)
merged_deg4434 <- merged_deg4434 %>%
  filter(p_val_adj_ADvsHD < 0.05 | p_val_adj_APOE < 0.05)


P <- ggplot(merged_deg, aes(x = avg_log2FC_ADvsHD, y = avg_log2FC_APOE)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    x = "Log2FC: AD vs HD",
    y = "Log2FC: APOE3/4 vs APOE3/3",
    title = "Comparación de DEGs en Monocitos"
  ) +
  theme(text = element_text(size = 14))
ggsave("APOE34E33ADHD.pdf", plot = P , width = 5, height = 6)
P <- ggplot(merged_deg4433, aes(x = avg_log2FC_ADvsHD, y = avg_log2FC_APOE)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    x = "Log2FC: AD vs HD",
    y = "Log2FC: APOE4/4 vs APOE3/3",
    title = "Comparación de DEGs en Monocitos"
  ) +
  theme(text = element_text(size = 14))
ggsave("APOE44E33ADHD.pdf", plot = P , width = 5, height = 6)


top_labels_34 <- Merged_DEG_ADvsHD_APOE34vs33 %>%
  slice_max(order_by = abs(avg_log2FC_ADvsHD), n = 20) %>%
  bind_rows(Merged_DEG_ADvsHD_APOE34vs33 %>% slice_max(order_by = abs(avg_log2FC_APOE), n = 20)) %>%
  distinct(gene, .keep_all = TRUE)

P1 <- ggplot(Merged_DEG_ADvsHD_APOE34vs33, aes(x = avg_log2FC_ADvsHD, y = avg_log2FC_APOE)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_text_repel(data = top_labels_34, aes(label = gene), size = 2, max.overlaps = 30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "A) APOE3/4 vs 3/3",
    x = "Log2FC: AD vs HD",
    y = "Log2FC: APOE3/4 vs APOE3/3"
  ) +
  theme(text = element_text(size = 7))

top_labels_44 <- Merged_DEG_ADvsHD_APOE44vs33 %>%
  slice_max(order_by = abs(avg_log2FC_ADvsHD), n = 20) %>%
  bind_rows(Merged_DEG_ADvsHD_APOE44vs33 %>% slice_max(order_by = abs(avg_log2FC_APOE), n = 20)) %>%
  distinct(gene, .keep_all = TRUE)

P2 <- ggplot(Merged_DEG_ADvsHD_APOE44vs33, aes(x = avg_log2FC_ADvsHD, y = avg_log2FC_APOE)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_text_repel(data = top_labels_44, aes(label = gene), size = 2, max.overlaps = 30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "B) APOE4/4 vs 3/3",
    x = "Log2FC: AD vs HD",
    y = "Log2FC: APOE4/4 vs APOE3/3"
  ) +
  theme(text = element_text(size = 7))

top_labels_4434 <-Merged_DEG_ADvsHD_APOE44vs34 %>%
  slice_max(order_by = abs(avg_log2FC_ADvsHD), n = 20) %>%
  bind_rows(Merged_DEG_ADvsHD_APOE44vs34 %>% slice_max(order_by = abs(avg_log2FC_APOE), n = 20)) %>%
  distinct(gene, .keep_all = TRUE)

P3 <- ggplot(Merged_DEG_ADvsHD_APOE44vs34, aes(x = avg_log2FC_ADvsHD, y = avg_log2FC_APOE)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_text_repel(data = top_labels_4434, aes(label = gene), size = 2, max.overlaps = 30) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "C) APOE4/4 vs 3/4",
    x = "Log2FC: AD vs HD",
    y = "Log2FC: APOE4/4 vs APOE3/4"
  ) +
  theme(text = element_text(size = 7))


combined_plot <- P1 | P2 | P3 + plot_layout(guides = "collect")


ggsave("DEG_comparison_APOE_panels.pdf", combined_plot, width = 15, height = 6)
write.csv(merged_deg3433, "Merged_DEG_ADvsHD_APOE34vs33.csv", row.names = FALSE)
write.csv(merged_deg4433, "Merged_DEG_ADvsHD_APOE44vs33.csv", row.names = FALSE)
write.csv(merged_deg4434, "Merged_DEG_ADvsHD_APOE44vs34.csv", row.names = FALSE)


merged_deg <- Merged_DEG_ADvsHD_APOE44vs33 %>% mutate(gene = rownames(Merged_DEG_ADvsHD_APOE44vs33))

q1_genes <- merged_deg %>%
  filter(avg_log2FC_ADvsHD > 0, avg_log2FC_APOE > 0) %>%
  pull(gene)
q1_indices <- as.numeric(q1_genes)
q1_gene_symbols <- rownames(mono_obj)[q1_indices]

# Convertir a Entrez IDs para enrichment
gene_df <- bitr(q1_gene_symbols,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

entrez_ids <- gene_df$ENTREZID
ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",             
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
barplot(ego, showCategory = 15, title = "GO BP Enrichment")

ggsave("GO_BP_enrichment_4433AD.pdf", width = 10, height = 6)

ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")
pdf("KEGG_barplot4434.pdf", width = 7, height = 6)
barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")
dev.off()
merged_deg <- Merged_DEG_ADvsHD_APOE34vs33 %>% mutate(gene = rownames(Merged_DEG_ADvsHD_APOE34vs33))

q1_genes <- merged_deg %>%
  filter(avg_log2FC_ADvsHD > 0, avg_log2FC_APOE > 0) %>%
  pull(gene)

q1_indices <- as.numeric(q1_genes)

q1_gene_symbols <- rownames(mono_obj)[q1_indices]

gene_df <- bitr(q1_gene_symbols,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

entrez_ids <- gene_df$ENTREZID
ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",             
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

barplot(ego, showCategory = 15, title = "GO BP Enrichment")
ggsave("GO_BP_enrichment_3433AD.pdf", width = 10, height = 6)

ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")

pdf("KEGG_barplot3433.pdf", width = 8, height = 6)
barplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")
dev.off()



