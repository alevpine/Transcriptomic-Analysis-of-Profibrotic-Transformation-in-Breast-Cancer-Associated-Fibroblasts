library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
library(tximport); library(DESeq2)
library(AnnotationDbi); library(org.Hs.eg.db)
library(ComplexHeatmap); library(circlize); library(grid)

DATA_DIR   <- "data"
RESULT_DIR <- "results"

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
COL_NF      <- "#2166AC"
COL_CAF     <- "#B2182B"
condition_colors <- c(NF = COL_NF, CAF = COL_CAF)

theme_clean <- theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = "grey30", linewidth = 0.3),
        legend.position = "bottom")

refseq_keys   <- keys(org.Hs.eg.db, keytype = "REFSEQ")
refseq2symbol <- AnnotationDbi::select(
  org.Hs.eg.db, keys = refseq_keys,
  columns = c("REFSEQ", "SYMBOL"), keytype = "REFSEQ")
tx2gene <- refseq2symbol %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(REFSEQ, .keep_all = TRUE) %>%
  dplyr::rename(TXNAME = REFSEQ, GENEID = SYMBOL)

patients <- paste0("P", 1:12)
sample_info <- data.frame(
  sample_id = c(paste0(patients, "_CAF"), paste0(patients, "_NF")),
  condition = factor(rep(c("CAF", "NF"), each = 12), levels = c("NF", "CAF")),
  patient   = factor(rep(patients, 2)), stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample_id
sample_info$files <- file.path(DATA_DIR, paste0(sample_info$sample_id, "_abundance.tsv"))
stopifnot(all(file.exists(sample_info$files)))

txi <- tximport(files = sample_info$files, type = "kallisto",
                tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "no")
colnames(txi$counts) <- sample_info$sample_id
colnames(txi$abundance) <- sample_info$sample_id
colnames(txi$length) <- sample_info$sample_id

dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ patient + condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]; dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "CAF", "NF"), alpha = PADJ_CUTOFF)
vsd <- vst(dds, blind = TRUE)
res_df <- as.data.frame(res) %>% mutate(gene = rownames(res)) %>% dplyr::filter(!is.na(padj))
write.csv(res_df %>% arrange(padj), file.path(RESULT_DIR, "DE_results.csv"), row.names = FALSE)
expr_mat <- assay(vsd)

profibrotic_genes <- list(
  "Myofibroblast" = c("ACTA2","TAGLN","MYL9","CNN1","FAP","PDPN","THY1","MYLK",
                      "TPM1","TPM2","DES","VIM","S100A4","PDGFRA","PDGFRB","CAV1","CTHRC1","POSTN"),
  "ECM" = c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL6A1","COL11A1","COL12A1",
            "FN1","LOX","LOXL2","TNC","SPARC","BGN","DCN","VCAN","LAMA4","LAMB1","ELN","FBLN1","THBS1","THBS2"),
  "MMPs / TIMPs" = c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP11","MMP13","MMP14",
                     "MMP16","ADAMTS1","ADAMTS2","ADAMTS5","TIMP1","TIMP2","TIMP3"),
  "TGF-b / SMAD" = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                     "SERPINE1","CTGF","LTBP1","LTBP2","LTBP4","BMP2","BMP4","INHBA","ACVR1","TGFBI","NOG","FST"),
  "IL-6 / JAK / STAT" = c("IL6","IL6R","IL6ST","JAK1","JAK2","STAT3","STAT1","SOCS1","SOCS3","IL11","LIF","OSM","IL10"),
  "NF-kB / Inflammation" = c("NFKB1","NFKB2","RELA","IKBKB","IL1A","IL1B","IL1R1","TNFAIP3",
                             "CCL2","CCL5","CCL7","CXCL1","CXCL2","CXCL5","CXCL8","CXCL12","CXCL14","PTGS2","ICAM1","VCAM1"),
  "EMT" = c("SNAI1","SNAI2","TWIST1","TWIST2","ZEB1","ZEB2","CDH1","CDH2","CDH11",
            "FOXC2","SOX4","SOX9","PRRX1","HMGA2","ITGA5","ITGB1"),
  "Wnt" = c("WNT2","WNT2B","WNT3A","WNT5A","WNT5B","WNT7A","WNT11","FZD1","FZD2",
            "FZD7","LRP5","LRP6","CTNNB1","AXIN2","LEF1","TCF7L2","DKK1","SFRP1","SFRP2","RSPO3"),
  "Hedgehog" = c("SHH","IHH","DHH","PTCH1","PTCH2","SMO","GLI1","GLI2","GLI3","SUFU","HHIP"),
  "YAP / Hippo" = c("YAP1","WWTR1","TEAD1","TEAD2","TEAD4","CYR61","CTGF","LATS1","LATS2","ANKRD1","AMOTL2","AJUBA"),
  "Angiogenesis" = c("VEGFA","VEGFB","VEGFC","PGF","FLT1","KDR","FLT4","PDGFA","PDGFB","PDGFC",
                     "FGF1","FGF2","FGF7","FGFR1","FGFR2","ANGPT1","ANGPT2","TEK","HIF1A","EPAS1","NRP1"),
  "Pre-met. niche" = c("SPP1","S100A4","HGF","MET","CXCR4","CSF1","LGALS3","COMP","SRGN"))

pathway_colors <- c("Myofibroblast"="#DC143C","ECM"="#4169E1","MMPs / TIMPs"="#9932CC",
                    "TGF-b / SMAD"="#FF8C00","IL-6 / JAK / STAT"="#DAA520","NF-kB / Inflammation"="#FF6347",
                    "EMT"="#8B008B","Wnt"="#228B22","Hedgehog"="#008080","YAP / Hippo"="#2F4F4F",
                    "Angiogenesis"="#00CED1","Pre-met. niche"="#8B4513")

subtype_sigs <- list(
  myCAF  = c("ACTA2","TAGLN","MYL9","MYLK","CNN1","TPM1","TPM2"),
  matCAF = c("CTHRC1","POSTN","COL11A1","COL12A1","COL10A1","COMP","SFRP4"),
  iCAF   = c("IL6","CXCL12","CCL2","CXCL14","IL1B","CXCL1","CXCL2"))
subtype_colors <- c(myCAF = "#DC143C", matCAF = "#4169E1", iCAF = "#DAA520")

core_fibrotic <- c("COL1A1","COL1A2","COL3A1","FN1","LOX","LOXL2","ACTA2","POSTN",
                   "CTHRC1","SERPINE1","CTGF","FAP","TNC","SPARC","BGN","THBS1")
core_fibrotic <- core_fibrotic[core_fibrotic %in% rownames(expr_mat)]

calc_score <- function(genes, mat) {
  genes <- genes[genes %in% rownames(mat)]
  if (length(genes) < 2) return(setNames(rep(NA, ncol(mat)), colnames(mat)))
  colMeans(t(scale(t(mat[genes, ]))), na.rm = TRUE)
}

scores_wide <- data.frame(
  sample = colnames(expr_mat),
  myCAF = calc_score(subtype_sigs$myCAF, expr_mat),
  matCAF = calc_score(subtype_sigs$matCAF, expr_mat),
  iCAF = calc_score(subtype_sigs$iCAF, expr_mat),
  condition = sample_info$condition, patient = sample_info$patient)
scores_wide$dominant <- apply(scores_wide[, c("myCAF","matCAF","iCAF")], 1,
                              function(x) c("myCAF","matCAF","iCAF")[which.max(x)])

# Fig 1a: PCA
pca <- prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
pct_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     condition = sample_info$condition, patient = sample_info$patient)
ggsave(file.path(RESULT_DIR, "Fig1a_PCA.png"),
       ggplot(pca_df, aes(PC1, PC2, color = condition)) + geom_point(size = 3) +
         stat_ellipse(level = 0.95, linewidth = 0.7, linetype = "dashed") +
         scale_color_manual(values = condition_colors, name = "") +
         labs(title = "PCA: CAF vs NF", x = sprintf("PC1 (%s%%)", pct_var[1]),
              y = sprintf("PC2 (%s%%)", pct_var[2])) +
         theme_clean + theme(aspect.ratio = 1), width = 6, height = 6, dpi = 300)

# Fig 1b: top DEGs
gene_pw_lookup <- function(g) {
  for (pw in names(profibrotic_genes)) if (g %in% profibrotic_genes[[pw]]) return(pw)
  return("Other")
}
top_up <- res_df %>% dplyr::filter(padj < PADJ_CUTOFF, log2FoldChange > LFC_CUTOFF) %>% arrange(padj) %>% head(15)
top_down <- res_df %>% dplyr::filter(padj < PADJ_CUTOFF, log2FoldChange < -LFC_CUTOFF) %>% arrange(padj) %>% head(15)
top_deg <- rbind(top_up, top_down) %>%
  mutate(pathway = sapply(gene, gene_pw_lookup),
         sig = ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", "*"))) %>%
  arrange(log2FoldChange) %>% mutate(gene = factor(gene, levels = gene))
ggsave(file.path(RESULT_DIR, "Fig1b_top_DEGs.png"),
       ggplot(top_deg, aes(x = gene, y = log2FoldChange, color = pathway)) +
         geom_segment(aes(xend = gene, y = 0, yend = log2FoldChange), linewidth = 0.7) +
         geom_point(size = 3) + geom_hline(yintercept = 0, linewidth = 0.4) +
         geom_text(aes(label = sig, y = log2FoldChange + ifelse(log2FoldChange >= 0, 0.2, -0.2)),
                   size = 3, color = "black") + coord_flip() +
         scale_color_manual(values = c(pathway_colors, "Other" = "grey50"), name = "Pathway") +
         labs(title = "Top differentially expressed genes: CAF vs NF",
              subtitle = "15 most significant up + 15 down", x = NULL, y = "log2 Fold Change") +
         theme_clean + theme(axis.text.y = element_text(size = 8, face = "italic"),
                             legend.position = "right", legend.text = element_text(size = 7)),
       width = 9, height = 8, dpi = 300)

# Fig 2: heatmap
gene_order <- c(); gene_to_pathway <- c()
for (pw in names(profibrotic_genes)) {
  pr <- profibrotic_genes[[pw]][profibrotic_genes[[pw]] %in% rownames(vsd)]
  pr <- pr[!pr %in% gene_order]
  if (length(pr) > 0) { gene_order <- c(gene_order, pr); gene_to_pathway <- c(gene_to_pathway, rep(pw, length(pr))) }
}
row_order <- c(paste0("P", 1:12, "_NF"), paste0("P", 1:12, "_CAF"))
row_order <- row_order[row_order %in% colnames(vsd)]
mat_t <- t(scale(t(expr_mat[gene_order, row_order]))); mat_t <- t(mat_t)
ht <- Heatmap(mat_t, name = "z-score",
              col = colorRamp2(c(-2, 0, 2), c(COL_NF, "#F7F7F7", COL_CAF)),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = factor(sample_info[row_order, "condition"], levels = c("NF","CAF")),
              column_split = factor(gene_to_pathway, levels = names(profibrotic_genes)),
              row_gap = unit(3,"mm"), column_gap = unit(1,"mm"),
              row_title_gp = gpar(fontsize = 11, fontface = "bold"), column_title = NULL,
              row_labels = gsub("_(NF|CAF)$", "", row_order),
              row_names_gp = gpar(fontsize = 9), row_names_side = "left",
              column_names_gp = gpar(fontsize = 5.5, fontface = "italic"), column_names_rot = 90,
              top_annotation = HeatmapAnnotation(Pathway = gene_to_pathway, col = list(Pathway = pathway_colors),
                                                 show_legend = TRUE, show_annotation_name = FALSE, height = unit(4,"mm"),
                                                 annotation_legend_param = list(Pathway = list(direction = "horizontal", nrow = 2))),
              left_annotation = rowAnnotation(Condition = sample_info[row_order, "condition"],
                                              col = list(Condition = condition_colors), show_legend = TRUE, show_annotation_name = FALSE,
                                              width = unit(4,"mm"), annotation_legend_param = list(Condition = list(direction = "horizontal", nrow = 1))),
              heatmap_legend_param = list(direction = "horizontal", legend_width = unit(4,"cm")),
              width = ncol(mat_t) * unit(2.5,"mm"), height = nrow(mat_t) * unit(6,"mm"))
png(file.path(RESULT_DIR, "Fig2_heatmap.png"), width = 30, height = 12, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
dev.off()

# Fig 3a: fibrotic score
fib_score <- calc_score(core_fibrotic, expr_mat)
fib_df <- data.frame(sample = names(fib_score), score = fib_score,
                     condition = sample_info$condition, patient = sample_info$patient)
pv_fib <- wilcox.test(fib_df$score[fib_df$condition == "CAF"],
                      fib_df$score[fib_df$condition == "NF"], paired = TRUE)$p.value
star_fib <- ifelse(pv_fib < 0.001, "***", ifelse(pv_fib < 0.01, "**", ifelse(pv_fib < 0.05, "*", "ns")))
ggsave(file.path(RESULT_DIR, "Fig3a_fibrotic_score.png"),
       ggplot(fib_df, aes(condition, score)) +
         geom_boxplot(aes(fill = condition), alpha = 0.3, outlier.shape = NA, width = 0.5) +
         geom_line(aes(group = patient), color = "grey55", linewidth = 0.35) +
         geom_point(aes(color = condition), size = 2.5) +
         scale_fill_manual(values = condition_colors) + scale_color_manual(values = condition_colors) +
         annotate("text", x = 1.5, y = max(fib_df$score) * 1.05, label = star_fib, size = 6, fontface = "bold") +
         labs(title = "Composite fibrotic score: CAF vs NF",
              subtitle = paste0("Mean z-score of ", length(core_fibrotic), " core fibrotic genes"),
              x = NULL, y = "Fibrotic score") +
         theme_clean + theme(legend.position = "none"), width = 5, height = 6, dpi = 300)

# Fig 3b: subtype scores
subtype_scores <- do.call(rbind, lapply(names(subtype_sigs), function(st)
  data.frame(subtype = st, sample = colnames(expr_mat), score = calc_score(subtype_sigs[[st]], expr_mat),
             condition = sample_info$condition, patient = sample_info$patient, stringsAsFactors = FALSE)))
subtype_scores$subtype <- factor(subtype_scores$subtype, levels = c("myCAF","matCAF","iCAF"))
st_pvals <- subtype_scores %>% group_by(subtype) %>%
  summarise(pval = wilcox.test(score[condition == "CAF"], score[condition == "NF"], paired = TRUE)$p.value,
            ymax = max(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(star = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns"))))
ggsave(file.path(RESULT_DIR, "Fig3b_subtype_scores.png"),
       ggplot(subtype_scores, aes(condition, score)) +
         geom_boxplot(aes(fill = condition), alpha = 0.3, outlier.shape = NA, width = 0.5) +
         geom_line(aes(group = patient), color = "grey55", linewidth = 0.35) +
         geom_point(aes(color = condition), size = 2.2) +
         geom_text(data = st_pvals, aes(x = 1.5, y = ymax * 1.06, label = star),
                   size = 5, fontface = "bold", inherit.aes = FALSE) +
         facet_wrap(~ subtype, scales = "free_y", nrow = 1) +
         scale_fill_manual(values = condition_colors) + scale_color_manual(values = condition_colors) +
         labs(title = "CAF subtype scores: CAF vs NF", subtitle = "Mean z-score of subtype marker genes",
              x = NULL, y = "Subtype score") +
         theme_clean + theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12)),
       width = 10, height = 5, dpi = 300)

# Fig 4: waterfall
all_profib <- unique(unlist(profibrotic_genes))
gene_pw_map <- data.frame(gene = character(), pathway = character(), stringsAsFactors = FALSE)
added <- c()
for (pw in names(profibrotic_genes)) {
  gs <- profibrotic_genes[[pw]]; gs <- gs[!gs %in% added]
  if (length(gs) > 0) { gene_pw_map <- rbind(gene_pw_map, data.frame(gene = gs, pathway = pw, stringsAsFactors = FALSE)); added <- c(added, gs) }
}
wf_df <- res_df %>% dplyr::filter(gene %in% all_profib) %>%
  left_join(gene_pw_map, by = "gene") %>% dplyr::filter(!is.na(pathway)) %>%
  mutate(sig = ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", ""))),
         pathway = factor(pathway, levels = names(profibrotic_genes))) %>%
  group_by(pathway) %>% mutate(gene_ordered = reorder(gene, log2FoldChange)) %>% ungroup()
ggsave(file.path(RESULT_DIR, "Fig4_waterfall.png"),
       ggplot(wf_df, aes(gene_ordered, log2FoldChange, fill = pathway)) +
         geom_col(width = 0.75, show.legend = FALSE) + geom_hline(yintercept = 0, linewidth = 0.3) +
         geom_text(aes(label = sig, y = log2FoldChange + ifelse(log2FoldChange >= 0, 0.15, -0.25)),
                   size = 2.2, color = "black") +
         coord_flip() + facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
         scale_fill_manual(values = pathway_colors) +
         labs(title = "Profibrotic genes: log2 fold change (CAF vs NF)",
              subtitle = "* padj<0.05  ** padj<0.01  *** padj<0.001", x = NULL, y = "log2 Fold Change") +
         theme_clean + theme(axis.text.y = element_text(size = 7, face = "italic"),
                             strip.text = element_text(face = "bold", size = 9),
                             strip.background = element_rect(fill = "grey95", color = NA)),
       width = 14, height = 16, dpi = 300)

# Fig 5: pathway activity
pathway_activity <- do.call(rbind, lapply(names(profibrotic_genes), function(pw) {
  gs <- profibrotic_genes[[pw]]; gs <- gs[gs %in% rownames(expr_mat)]
  if (length(gs) < 2) return(NULL)
  act <- colMeans(t(scale(t(expr_mat[gs, ]))), na.rm = TRUE)
  data.frame(pathway = pw, sample = names(act), activity = act,
             condition = sample_info[names(act), "condition"],
             patient = sample_info[names(act), "patient"], stringsAsFactors = FALSE)
}))
pathway_activity$pathway <- factor(pathway_activity$pathway, levels = names(profibrotic_genes))
pathway_activity$dominant <- scores_wide[pathway_activity$sample, "dominant"]
pw_pvals <- pathway_activity %>% group_by(pathway) %>%
  summarise(pval = wilcox.test(activity[condition == "CAF"], activity[condition == "NF"], paired = TRUE)$p.value,
            ymax = max(activity), .groups = "drop") %>%
  mutate(star = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns"))),
         pathway = factor(pathway, levels = names(profibrotic_genes)))
ggsave(file.path(RESULT_DIR, "Fig5_pathway_activity.png"),
       ggplot(pathway_activity, aes(condition, activity)) +
         geom_boxplot(aes(fill = condition), alpha = 0.3, outlier.shape = NA, width = 0.5) +
         geom_line(aes(group = patient), color = "grey55", linewidth = 0.3) +
         geom_point(aes(color = condition), size = 1.5) +
         geom_text(data = pw_pvals, aes(x = 1.5, y = ymax * 1.08, label = star),
                   size = 4, fontface = "bold", inherit.aes = FALSE) +
         facet_wrap(~ pathway, scales = "free_y", ncol = 4) +
         scale_fill_manual(values = condition_colors) + scale_color_manual(values = condition_colors) +
         labs(title = "Pathway activity: CAF vs NF", subtitle = "Wilcoxon signed-rank paired test",
              x = NULL, y = "Pathway activity (mean z-score)") +
         theme_clean + theme(legend.position = "none", strip.text = element_text(face = "bold", size = 9)),
       width = 12, height = 10, dpi = 300)

# Fig 6: pathway x subtype
pw_caf <- pathway_activity %>% dplyr::filter(condition == "CAF")
pw_summary <- pw_caf %>% group_by(pathway, dominant) %>%
  summarise(mean_act = mean(activity, na.rm = TRUE), .groups = "drop") %>%
  mutate(pathway = factor(pathway, levels = rev(names(profibrotic_genes))),
         dominant = factor(dominant, levels = c("myCAF","matCAF","iCAF")))
pw_kw <- pw_caf %>% group_by(pathway) %>%
  summarise(kw_p = tryCatch(kruskal.test(activity ~ dominant)$p.value, error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(kw_star = ifelse(is.na(kw_p), "", ifelse(kw_p < 0.001, "***",
                                                  ifelse(kw_p < 0.01, "**", ifelse(kw_p < 0.05, "*", "")))),
         pathway = factor(pathway, levels = rev(names(profibrotic_genes))))
ggsave(file.path(RESULT_DIR, "Fig6_pathway_by_subtype.png"),
       ggplot(pw_summary, aes(x = dominant, y = pathway, fill = mean_act)) +
         geom_tile(color = "white", linewidth = 1.2) +
         geom_text(aes(label = sprintf("%.2f", mean_act)), size = 3.5, color = "black") +
         geom_text(data = pw_kw %>% dplyr::filter(kw_star != ""),
                   aes(x = 3.6, y = pathway, label = kw_star), size = 3.5, color = "grey20", inherit.aes = FALSE) +
         scale_fill_gradient2(low = COL_NF, mid = "white", high = COL_CAF, midpoint = 0, name = "Mean\nactivity") +
         scale_x_discrete(position = "top") + coord_cartesian(xlim = c(0.5, 4.2), clip = "off") +
         labs(title = "Pathway activity in CAFs by dominant subtype",
              subtitle = "Kruskal-Wallis significance on the right", x = NULL, y = NULL) +
         theme_clean + theme(axis.text.y = element_text(size = 9, face = "bold"),
                             axis.text.x.top = element_text(size = 10, face = "bold"),
                             axis.line = element_blank(), legend.position = "right",
                             plot.margin = margin(5, 30, 5, 5)),
       width = 10, height = 7, dpi = 300)
