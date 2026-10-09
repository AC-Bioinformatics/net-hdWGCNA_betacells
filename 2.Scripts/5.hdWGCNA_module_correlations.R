# hdWGCNA part 3: what do the modules track? Module eigengenes (MEs, computed on the
# v3 ND/T2D single cells in 4.hdWGCNA_network.R) are related to:
#   page 1 - the authors' beta-cell subclusters (C0-C7)
#   page 2 - donor traits (HbA1c, BMI, age, sex; per donor) and cell QC metrics (per cell)
#   page 3 - condition: ND vs T2D, per donor
# Donor-level questions use one value per donor (mean ME over the donor's cells), so
# 13 ND vs 10 T2D donors are compared - not ~49,000 cells, which would make almost
# any difference "significant".
Root = '/home/alan.culligan/network'
Data_dir = file.path(Root, '1.Data', 'GSE221156')
Out_dir = file.path(Root, '3.Results')
setwd(Root)

library(Seurat)
library(tidyverse)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

seu <- readRDS(file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
md  <- seu[[]]

# Module eigengenes (no batch correction - all cells are v3), grey = unassigned genes
MEs  <- GetMEs(seu, harmonized = FALSE)
mods <- setdiff(colnames(MEs), "grey")
MEs  <- as.matrix(MEs[, mods])
modules <- GetModules(seu)
print(table(modules$module))

# Hub genes: connectivity of every gene to every module (kME), top 10 per module.
# reassign_modules = FALSE keeps the modules exactly as the network built them.
seu <- ModuleConnectivity(seu, harmonized = FALSE, reassign_modules = FALSE)
hubs <- GetHubGenes(seu, n_hubs = 10)
write.csv(GetModules(seu), file.path(Out_dir, "10.Modules_kME.csv"), row.names = FALSE)
write.csv(hubs, file.path(Out_dir, "10.Hub_genes.csv"), row.names = FALSE)
hub_lab <- hubs %>% mutate(module = as.character(module)) %>% group_by(module) %>%
  summarise(hubs = paste(head(gene_name, 5), collapse = ", "))
mod_lab <- setNames(paste0(mods, " (", table(modules$module)[mods], ")"), mods)

# Heatmap helper: modules x variables, coloured by correlation, labelled with r (and stars)
heat <- function(r, p = NULL, title, subtitle, xlab) {
  df <- as.data.frame(as.table(r)); colnames(df) <- c("module", "var", "r")
  df$label <- sprintf("%.2f", df$r)
  if (!is.null(p)) {
    stars <- cut(as.vector(p), c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
    df$label <- paste0(df$label, stars)
  }
  df$module <- factor(mod_lab[as.character(df$module)], levels = rev(mod_lab))
  ggplot(df, aes(var, module, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "r") +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

## Page 1: authors' subclusters - correlation of each ME with membership of each cluster (cells)
clusters <- sort(unique(as.character(md$Clusters)))
r_clus <- sapply(clusters, function(k) cor(MEs, as.numeric(md$Clusters == k)))
rownames(r_clus) <- mods
write.csv(r_clus, file.path(Out_dir, "10.Module_cluster_cor.csv"))
p1 <- heat(r_clus, title = "Modules vs the authors' beta-cell subclusters",
           subtitle = paste0("Pearson r between each cell's module eigengene and membership of each cluster (", nrow(md), " v3 ND/T2D cells)\n",
                             "C0 = insulin secretion, C4 = CD63-high / low ribo, C6 = senescence (paper clusters 1, 5, 7)"),
           xlab = "Author cluster")

## Page 2: donor traits (per donor) and cell QC (per cell)
donor_ME <- aggregate(MEs, by = list(donor = as.character(md$donor)), FUN = mean)
donor_md <- md %>% distinct(donor, .keep_all = TRUE) %>%
  transmute(donor = as.character(donor), condition, HbA1c, BMI, age, male = as.numeric(sex == "male"))
dd <- merge(donor_ME, donor_md, by = "donor")
traits <- c("HbA1c", "BMI", "age", "male")
r_tr <- sapply(traits, function(t) sapply(mods, function(m) cor(dd[[m]], dd[[t]], method = "spearman", use = "pairwise")))
p_tr <- sapply(traits, function(t) sapply(mods, function(m) suppressWarnings(cor.test(dd[[m]], dd[[t]], method = "spearman", exact = FALSE)$p.value)))
write.csv(cbind(r_tr, setNames(as.data.frame(p_tr), paste0(traits, "_p"))), file.path(Out_dir, "10.Module_trait_cor_donor.csv"))
p2a <- heat(r_tr, p_tr, title = "Modules vs donor traits",
            subtitle = paste0("Spearman r of donor-mean eigengene with each trait, ", nrow(dd), " donors (* p<0.05, ** p<0.01, *** p<0.001, unadjusted); male = 1"),
            xlab = "Donor trait")
qc <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
r_qc <- sapply(qc, function(q) cor(MEs, md[[q]], method = "spearman"))
rownames(r_qc) <- mods
write.csv(r_qc, file.path(Out_dir, "10.Module_QC_cor_cells.csv"))
p2b <- heat(r_qc, title = "Modules vs cell QC metrics",
            subtitle = "Spearman r per cell - strong correlations flag modules driven by depth / mito / ribo rather than biology",
            xlab = "QC metric (authors' nCount / nFeature / percent.mt)")

## Page 3: condition - donor-mean eigengene, ND vs T2D
long <- dd %>% select(donor, condition, all_of(mods)) %>%
  pivot_longer(all_of(mods), names_to = "module", values_to = "ME")
stats <- long %>% group_by(module) %>%
  summarise(ND_median = median(ME[condition == "ND"]), T2D_median = median(ME[condition == "T2D"]),
            wilcox_p = wilcox.test(ME[condition == "T2D"], ME[condition == "ND"], exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(wilcox_p_adj = p.adjust(wilcox_p, "BH"))
# with age and sex as covariates (linear model on donor means)
stats$lm_p_age_sex <- sapply(stats$module, function(m) {
  d <- merge(dd[, c("donor", m)], donor_md, by = "donor"); colnames(d)[2] <- "ME"
  coef(summary(lm(ME ~ condition + age + male, data = d)))["conditionT2D", "Pr(>|t|)"]
})
stats <- left_join(stats, hub_lab, by = "module")
write.csv(stats, file.path(Out_dir, "10.Module_condition_donor.csv"), row.names = FALSE)
print(stats)
long$module <- factor(long$module, levels = mods,
                      labels = sprintf("%s\nWilcoxon p=%.2g (BH %.2g)", mods,
                                       stats$wilcox_p[match(mods, stats$module)], stats$wilcox_p_adj[match(mods, stats$module)]))
p3 <- ggplot(long, aes(condition, ME, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  facet_wrap(~ module, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(ND = "#F8766D", T2D = "#00BFC4")) +
  labs(title = "Module eigengenes: ND vs T2D", y = "Donor-mean module eigengene", x = NULL,
       subtitle = "One point per donor (13 ND, 10 T2D, v3); age/sex-adjusted p-values in 10.Module_condition_donor.csv") +
  theme_bw() + theme(legend.position = "none", strip.text = element_text(size = 8))

pdf(file.path(Out_dir, "10.Module_correlations.pdf"), width = 14, height = 10)
print(p1)
print(p2a / p2b)
print(p3)
dev.off()
message("Saved: ", file.path(Out_dir, "10.Module_correlations.pdf"))

saveRDS(seu, file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
message("Saved (with kME): ", file.path(Data_dir, "GSE221156_beta_hdwgcna.rds"))
