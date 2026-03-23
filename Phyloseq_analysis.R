# 0. Load required libraries

library(dada2)
library(phyloseq)

# 1. Set paths and read input files

base_dir <- "C:/Users/kayle/Combined_Sep_data/Kaylee results"

seqtab_J <- readRDS(file.path(base_dir, "Jseqtab.nochim.rds"))
seqtab_N <- readRDS(file.path(base_dir, "Nseqtab.nochim.rds"))

taxa_J  <- readRDS(file.path(base_dir, "Jtaxa_assigned.rds"))
taxa_N  <- readRDS(file.path(base_dir, "Ntaxa_assigned.rds"))

track_J <- readRDS(file.path(base_dir, "Jtrack_reads.rds"))
track_N <- readRDS(file.path(base_dir, "Ntrack_reads.rds"))

metadata <- read.csv(file.path(base_dir, "173Combined_MetadataF.csv"), 
                     row.names = 1, check.names = FALSE)

# 2. Merge sequence tables

seqtab_all <- mergeSequenceTables(seqtab_J, seqtab_N)
cat("Merged sequence table dimensions (samples x ASVs):", dim(seqtab_all), "\n")


# 3. Check sample name consistency

cat("Samples in seqtab not in metadata:\n")
print(setdiff(rownames(seqtab_all), rownames(metadata)))
cat("Samples in metadata not in seqtab:\n")
print(setdiff(rownames(metadata), rownames(seqtab_all)))

# 4. Merge taxonomy tables

taxa_all <- taxa_J
new_seqs <- setdiff(rownames(taxa_N), rownames(taxa_all))

if (length(new_seqs) > 0) {
  taxa_all <- rbind(taxa_all, taxa_N[new_seqs, , drop = FALSE])
  cat("Added", length(new_seqs), "new taxa from second run.\n")
}

asv_seqs <- colnames(seqtab_all)
missing_taxa <- setdiff(asv_seqs, rownames(taxa_all))

if (length(missing_taxa) > 0) {
  cat("WARNING:", length(missing_taxa), "ASVs in merged seqtab have no taxonomy. Filling with NA.\n")
  ntax <- ncol(taxa_all)
  na_rows <- matrix(NA, nrow = length(missing_taxa), ncol = ntax,
                    dimnames = list(missing_taxa, colnames(taxa_all)))
  taxa_all <- rbind(taxa_all, na_rows)
}

taxa_all <- taxa_all[colnames(seqtab_all), , drop = FALSE]


# 5. Merge track-reads tables

intersecting_samples <- intersect(rownames(track_J), rownames(track_N))
if (length(intersecting_samples) > 0) {
  cat("WARNING: Duplicate sample names found between runs:\n")
  print(intersecting_samples)
}

track_all_df <- rbind(track_J, track_N)

# 6. Build phyloseq object

OTU <- otu_table(seqtab_all, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxa_all))
SAM <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAM)

# 7. Save merged outputs

saveRDS(seqtab_all, file = file.path(base_dir, "Merged_seqtab_nochim.rds"))
write.csv(seqtab_all, file = file.path(base_dir, "Merged_seqtab_nochim.csv"))

saveRDS(taxa_all, file = file.path(base_dir, "Merged_taxa_assigned.rds"))
write.csv(taxa_all, file = file.path(base_dir, "Merged_taxa_assigned.csv"))

saveRDS(track_all_df, file = file.path(base_dir, "Merged_track_reads.rds"))
write.csv(track_all_df, file = file.path(base_dir, "Merged_track_reads.csv"))

saveRDS(ps, file = file.path(base_dir, "Merged_phyloseq_object.rds"))

cat("\nMerging completed successfully.\n")
cat("Saved files:\n")
cat("- Merged_seqtab_nochim.rds / .csv\n")
cat("- Merged_taxa_assigned.rds / .csv\n")
cat("- Merged_track_reads.rds / .csv\n")
cat("- Merged_phyloseq_object.rds\n")

# ================================================================
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


#check phyloseq objects
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

#check phyloseq objects
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

