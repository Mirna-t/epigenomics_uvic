# JQ1-treated SCC methylation practical
# Ready-to-run template based on the source files provided in the chat.
#
# Main source mapping:
# - Raw import / QC / detectionP / preprocessQuantile:
#   script_2025_EPIC_modif.R (unpaginated), Minfi.pdf pp. 1365-1366,
#   02_Methylation_workshop_1.pdf PDF p. 16
# - Probe filtering rationale:
#   Normalization_review.pdf pp. 930-931; SWAN.pdf p. 7
# - Modelling with limma on M-values:
#   ChAMP.pdf p. 429
# - Promoter / promoter+island summaries:
#   paper.pdf pp. 62548, 62552; script_2025_EPIC_modif.R
#
# IMPORTANT:
# 1) Edit idat_dir to the folder that contains the IDAT files and sample sheet.
# 2) Make sure your sample sheet contains Treatment and CellLine columns.
#    Treatment should have DMSO and JQ1.
#    CellLine should have SCC15 and SCC25.

suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(missMethyl)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ChAMP)         # for rmSNPandCH
  library(tibble)
})

# -----------------------------
# 1. INPUT
# -----------------------------
idat_dir <- "PATH_TO_IDAT_FOLDER"
out_dir  <- "JQ1_practical_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

targets <- read.metharray.sheet(base = idat_dir)

# Harmonize metadata names if necessary
# Edit these lines if your sample sheet uses different column names.
if (!"Treatment" %in% colnames(targets)) stop("Sample sheet must contain a 'Treatment' column.")
if (!"CellLine"  %in% colnames(targets)) stop("Sample sheet must contain a 'CellLine' column.")

targets$Treatment <- factor(targets$Treatment, levels = c("DMSO", "JQ1"))
targets$CellLine  <- factor(targets$CellLine, levels = c("SCC15", "SCC25"))

rgset <- read.metharray.exp(targets = targets, verbose = TRUE)

# -----------------------------
# 2. RAW QC
# -----------------------------
mset_raw <- preprocessRaw(rgset)
qc <- getQC(mset_raw)

pdf(file.path(out_dir, "QC_plotQC.pdf"), width = 7, height = 6)
plotQC(qc)
dev.off()

pdf(file.path(out_dir, "QC_beta_density_raw.pdf"), width = 8, height = 6)
densityPlot(mset_raw, sampGroups = targets$Treatment, main = "Raw beta density by treatment")
dev.off()

detP <- detectionP(rgset)

pdf(file.path(out_dir, "QC_detectionP_barplot.pdf"), width = 10, height = 5)
barplot(colMeans(detP),
        col = as.numeric(targets$Treatment),
        las = 2, cex.names = 0.8,
        main = "Mean detection p-value per sample")
legend("topright", legend = levels(targets$Treatment),
       fill = seq_along(levels(targets$Treatment)), bty = "n")
dev.off()

qc_sample_table <- data.frame(
  Sample = sampleNames(rgset),
  Mean_DetectionP = colMeans(detP),
  Treatment = targets$Treatment,
  CellLine = targets$CellLine
)
write.csv(qc_sample_table, file.path(out_dir, "QC_sample_metrics.csv"), row.names = FALSE)

# Remove poor-quality samples
keep_samples <- colMeans(detP) < 0.01
rgset <- rgset[, keep_samples]
targets <- targets[keep_samples, , drop = FALSE]
detP <- detP[, keep_samples, drop = FALSE]

# Remove poor-quality probes
keep_probes <- rowSums(detP < 0.01) == ncol(rgset)
rgset <- rgset[keep_probes, ]

writeLines(
  c(
    paste("Samples retained:", ncol(rgset)),
    paste("Probes retained after detectionP filter:", nrow(rgset))
  ),
  con = file.path(out_dir, "QC_retention_summary.txt")
)

# -----------------------------
# 3. NORMALIZATION
# -----------------------------
# Source: Minfi.pdf p. 1365
# preprocessQuantile is suitable when strong global methylation shifts are not expected.
# If your exploratory plots suggest strong global distribution shifts, consider preprocessFunnorm instead.
grset <- preprocessQuantile(rgset)

beta <- getBeta(grset)
Mval <- getM(grset)

# -----------------------------
# 4. REMOVE SNP / CROSS-REACTIVE / SEX-CHROMOSOME PROBES
# -----------------------------
# Source: script_2025_EPIC_modif.R (unpaginated)
Mval_filt <- rmSNPandCH(Mval,
                        dist = 2,
                        mafcut = 0.05,
                        rmcrosshyb = TRUE,
                        rmXY = TRUE)

beta_filt <- beta[rownames(Mval_filt), , drop = FALSE]

writeLines(
  paste("Probes retained after rmSNPandCH:", nrow(Mval_filt)),
  con = file.path(out_dir, "Probe_filter_summary.txt")
)

# -----------------------------
# 5. PCA: TREATMENT VS CELL LINE
# -----------------------------
# Source: 02_Methylation_workshop_1.pdf PDF pp. 15-16
var_idx <- order(apply(beta_filt, 1, var), decreasing = TRUE)[1:min(10000, nrow(beta_filt))]
pca <- prcomp(t(beta_filt[var_idx, ]), center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(beta_filt),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Treatment = targets$Treatment,
  CellLine = targets$CellLine
)

var_exp <- round(100 * summary(pca)$importance[2, 1:2], 2)

p1 <- ggplot(pca_df, aes(PC1, PC2, color = Treatment, shape = CellLine)) +
  geom_point(size = 3) +
  labs(title = "PCA of normalized beta values",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw()

p2 <- ggplot(pca_df, aes(PC1, PC2, color = CellLine, shape = Treatment)) +
  geom_point(size = 3) +
  labs(title = "Same PCA, emphasizing CellLine",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw()

ggsave(file.path(out_dir, "PCA_treatment_cellline_1.pdf"), p1, width = 7, height = 5)
ggsave(file.path(out_dir, "PCA_treatment_cellline_2.pdf"), p2, width = 7, height = 5)
write.csv(pca_df, file.path(out_dir, "PCA_scores.csv"), row.names = FALSE)

# Simple PC association checks to justify model choice
pc1_lm <- summary(lm(PC1 ~ CellLine + Treatment, data = pca_df))
pc2_lm <- summary(lm(PC2 ~ CellLine + Treatment, data = pca_df))
capture.output(pc1_lm, file = file.path(out_dir, "PC1_model_summary.txt"))
capture.output(pc2_lm, file = file.path(out_dir, "PC2_model_summary.txt"))

# -----------------------------
# 6. LIMMA DIFFERENTIAL METHYLATION
# -----------------------------
# Source: ChAMP.pdf p. 429; script_2025_EPIC_modif.R
design_treat_only <- model.matrix(~ Treatment, data = targets)
design_with_cellline <- model.matrix(~ CellLine + Treatment, data = targets)

fit_treat_only <- eBayes(lmFit(Mval_filt, design_treat_only))
fit_with_cellline <- eBayes(lmFit(Mval_filt, design_with_cellline))

res_treat_only <- topTable(fit_treat_only, coef = "TreatmentJQ1", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name")

res <- topTable(fit_with_cellline, coef = "TreatmentJQ1", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name")

sig_res <- res %>% filter(adj.P.Val < 0.05)

write.csv(res_treat_only, file.path(out_dir, "DMP_results_treatment_only.csv"), row.names = FALSE)
write.csv(res, file.path(out_dir, "DMP_results_cellline_adjusted.csv"), row.names = FALSE)
write.csv(sig_res, file.path(out_dir, "DMP_results_cellline_adjusted_significant.csv"), row.names = FALSE)

# -----------------------------
# 7. ANNOTATION
# -----------------------------
# Promoter definition based on provided script and paper.pdf:
# TSS1500 / TSS200 / 5'UTR / 1stExon
ann <- as.data.frame(getAnnotation(grset))

res_annot <- res %>%
  left_join(ann %>% dplyr::select(Name, UCSC_RefGene_Name,
                                  UCSC_RefGene_Group, Relation_to_Island),
            by = "Name")

sig_annot <- res_annot %>% filter(adj.P.Val < 0.05)

promoter_sig <- sig_annot %>%
  filter(grepl("TSS1500|TSS200|5'UTR|1stExon", UCSC_RefGene_Group))

promoter_island_sig <- promoter_sig %>%
  filter(Relation_to_Island == "Island")

write.csv(res_annot, file.path(out_dir, "DMP_results_annotated.csv"), row.names = FALSE)
write.csv(promoter_sig, file.path(out_dir, "DMP_promoter_significant.csv"), row.names = FALSE)
write.csv(promoter_island_sig, file.path(out_dir, "DMP_promoter_island_significant.csv"), row.names = FALSE)

# Gene lists
split_genes <- function(x) {
  ux <- unique(unlist(strsplit(paste(na.omit(x), collapse = ";"), ";")))
  ux[nchar(ux) > 0]
}

genes_all <- split_genes(sig_annot$UCSC_RefGene_Name)
genes_promoter <- split_genes(promoter_sig$UCSC_RefGene_Name)
genes_promoter_island <- split_genes(promoter_island_sig$UCSC_RefGene_Name)

write.table(genes_all, file.path(out_dir, "Genes_all_significant.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_promoter, file.path(out_dir, "Genes_promoter_significant.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_promoter_island, file.path(out_dir, "Genes_promoter_island_significant.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

summary_counts <- data.frame(
  Metric = c("All significant DMPs", "Promoter significant DMPs", "Promoter+Island significant DMPs",
             "Unique genes near all significant DMPs", "Unique genes near promoter significant DMPs",
             "Unique genes near promoter+island significant DMPs"),
  Value = c(nrow(sig_annot), nrow(promoter_sig), nrow(promoter_island_sig),
            length(genes_all), length(genes_promoter), length(genes_promoter_island))
)
write.csv(summary_counts, file.path(out_dir, "Summary_counts.csv"), row.names = FALSE)

# -----------------------------
# 8. ENRICHMENT
# -----------------------------
# missMethyl is preferred for array data because it handles probe-number bias.
sig_cpg <- sig_annot$Name
all_cpg <- rownames(Mval_filt)

if (length(sig_cpg) > 0) {
  go_res <- gometh(sig.cpg = sig_cpg, all.cpg = all_cpg, collection = "GO")
  kegg_res <- gometh(sig.cpg = sig_cpg, all.cpg = all_cpg, collection = "KEGG")
  write.csv(go_res, file.path(out_dir, "GO_enrichment.csv"), row.names = FALSE)
  write.csv(kegg_res, file.path(out_dir, "KEGG_enrichment.csv"), row.names = FALSE)
}

# -----------------------------
# 9. OPTIONAL HEATMAP OF TOP DMPs
# -----------------------------
top_heat <- sig_annot %>%
  slice_min(order_by = adj.P.Val, n = min(100, n())) %>%
  pull(Name)

if (length(top_heat) > 1) {
  mat <- beta_filt[top_heat, , drop = FALSE]
  zmat <- t(scale(t(mat)))
  pdf(file.path(out_dir, "Heatmap_top100_DMPs.pdf"), width = 8, height = 10)
  heatmap(zmat, scale = "none", Colv = NA, col = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
}

# -----------------------------
# 10. SHORT REPORTING TEXT
# -----------------------------
report_lines <- c(
  paste("Final sample count:", ncol(Mval_filt)),
  paste("Final probe count:", nrow(Mval_filt)),
  paste("Significant DMPs (FDR < 0.05):", nrow(sig_annot)),
  paste("Promoter DMPs:", nrow(promoter_sig)),
  paste("Promoter+Island DMPs:", nrow(promoter_island_sig)),
  "",
  "Interpretation checklist:",
  "- Does PCA separate samples mainly by CellLine or by Treatment?",
  "- If CellLine drives PC1/PC2, keep CellLine in the final model.",
  "- Discuss whether JQ1 effects remain after CellLine adjustment.",
  "- Link enriched terms to BRD4 / JQ1 biology only after checking the provided background papers."
)
writeLines(report_lines, con = file.path(out_dir, "Report_ready_summary.txt"))

message("Analysis template completed. Outputs written to: ", normalizePath(out_dir))
