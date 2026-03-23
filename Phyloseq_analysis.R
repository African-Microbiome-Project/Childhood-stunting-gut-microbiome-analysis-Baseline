# 0. Load libraries
library(phyloseq)
library(magrittr)
library(vegan)

# 1. Set base directory

base_dir <- "C:/Users/kayle/Combined_Sep_data/Kaylee results"

# 2. Load merged sequence & taxonomy tables

seqtab_all <- readRDS(file.path(base_dir, "Merged_seqtab_nochim_ALL.rds"))
taxa_all <- readRDS(file.path(base_dir, "Merged_taxa_assigned_ALL.rds"))

# 3. Load Informal metadata

metadata_informal <- read.csv(file.path(base_dir, "Combined_Metadata_Informal.csv"),
                              row.names = 1, check.names = FALSE)


# 4. Subset seq table and taxa to only Informal samples

common_samples <- intersect(rownames(seqtab_all), rownames(metadata_informal))
seqtab_informal <- seqtab_all[common_samples, , drop = FALSE]
taxa_informal <- taxa_all[colnames(seqtab_informal), , drop = FALSE]

# 5. Build phyloseq object

OTU_informal <- otu_table(seqtab_informal, taxa_are_rows = FALSE)
TAX_informal <- tax_table(as.matrix(taxa_informal))
SAM_informal <- sample_data(metadata_informal)
ps_informal <- phyloseq(OTU_informal, TAX_informal, SAM_informal)

# 6. Save initial merged phyloseq object

saveRDS(seqtab_informal, file.path(base_dir, "Merged_seqtab_nochim_Informal.rds"))
write.csv(seqtab_informal, file.path(base_dir, "Merged_seqtab_nochim_Informal.csv"))

saveRDS(taxa_informal, file.path(base_dir, "Merged_taxa_assigned_Informal.rds"))
write.csv(taxa_informal, file.path(base_dir, "Merged_taxa_assigned_Informal.csv"))

saveRDS(ps_informal, file.path(base_dir, "Merged_phyloseq_object_Informal.rds"))
cat("Merged outputs saved for Informal metadata.\n")

#check phyloseq objects)
ps_informal
table(tax_table(ps_informal)[, "Phylum"], exclude = NULL)

# 7. Clean phyloseq object

ps_informal_clean <- ps_informal %>%
  subset_taxa(!is.na(Phylum) & Phylum != "") %>%
  subset_taxa(Family != "Mitochondria") %>%
  subset_taxa(Order != "Chloroplast")

saveRDS(ps_informal_clean, file.path(base_dir, "Merged_phyloseq_object_Informal_clean.rds"))
cat("Cleaned phyloseq object (Informal) saved.\n")

#check phyloseq objects)
ps_informal_clean
table(tax_table(ps_informal_clean)[, "Phylum"], exclude = NULL)

# 8. Unrarefied rarefaction curves

otu_mat_informal <- as(otu_table(ps_informal_clean), "matrix")
if(taxa_are_rows(ps_informal_clean)) otu_mat_informal <- t(otu_mat_informal)

# Plot to console
rarecurve(otu_mat_informal, step = 1000, col = "steelblue", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Unrarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")

# Save PDF
pdf(file.path(base_dir, "Informal_unrarefied_rarefaction.pdf"), width = 7, height = 5)
rarecurve(otu_mat_informal, step = 1000, col = "steelblue", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Unrarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")
dev.off()

# Save PNG
png(file.path(base_dir, "Informal_unrarefied_rarefaction.png"), width = 7*600, height = 5*600, res = 600)
rarecurve(otu_mat_informal, step = 1000, col = "steelblue", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Unrarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")
dev.off()


# 9. Rarefy phyloseq object

rarefy_depth_informal <- 176000
set.seed(123)
ps_informal_rarefied <- rarefy_even_depth(ps_informal_clean,
                                          sample.size = rarefy_depth_informal,
                                          rngseed = 123,
                                          replace = FALSE)

saveRDS(ps_informal_rarefied, file.path(base_dir, "Merged_rarefied_Informal.rds"))
cat("Rarefied phyloseq object (Informal) saved.\n")

#sanity check ( check phyloseq objects)
ps_informal_rarefied
table(tax_table(ps_informal_rarefied)[, "Phylum"], exclude = NULL)


# 10. Rarefied rarefaction curves

otu_mat_informal_r <- as(otu_table(ps_informal_rarefied), "matrix")
if(taxa_are_rows(ps_informal_rarefied)) otu_mat_informal_r <- t(otu_mat_informal_r)

# Plot to console
rarecurve(otu_mat_informal_r, step = 1000, col = "darkgreen", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Rarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")

# Save PDF
pdf(file.path(base_dir, "Informal_rarefied_rarefaction.pdf"), width = 7, height = 5)
rarecurve(otu_mat_informal_r, step = 1000, col = "darkgreen", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Rarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")
dev.off()

# Save PNG
png(file.path(base_dir, "Informal_rarefied_rarefaction.png"), width = 7*600, height = 5*600, res = 600)
rarecurve(otu_mat_informal_r, step = 1000, col = "darkgreen", cex = 0.6,
          label = FALSE, main = "Rarefaction Curves – Informal (Rarefied)",
          xlab = "Sequencing Depth (Reads)", ylab = "Observed ASVs")
dev.off()
