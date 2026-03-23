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

# 0. Load libraries
library(magrittr)  


# 1. Set base directory

base_dir <- "C:/Users/kayle/Combined_Sep_data/Kaylee results"

# 2. Load merged phyloseq object

ps <- readRDS(file.path(base_dir, "Merged_phyloseq_object.rds"))

# 3. Remove unwanted taxa

ps_clean <- ps %>%
  subset_taxa(!is.na(Phylum) & Phylum != "") %>%   # remove NA or blank Phylum
  subset_taxa(Family != "Mitochondria") %>%       # remove mitochondria
  subset_taxa(Order != "Chloroplast")             # remove chloroplast

# 4. Quick QC: check remaining samples and ASVs

cat("Samples:", nsamples(ps_clean), "\n")
cat("ASVs:", ntaxa(ps_clean), "\n")

# 5. Save cleaned phyloseq object

saveRDS(ps_clean, file = file.path(base_dir, "Merged_phyloseq_object_clean.rds"))
