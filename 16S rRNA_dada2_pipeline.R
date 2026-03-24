# Load the necessary packages
library(dada2)

# Define the path to your combined fastq files
path <- "/mnt/lustre/users/nkambunga/NRehousing_nodup"

# Define the path to your combined samples
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz$", full.names = TRUE)) # Forawrd reads
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz$", full.names = TRUE)) # Reverse reads

# Check the number of files (R1 and R2 should match)
if (length(fnFs) != length(fnRs)) {
  stop("The number of forward and reverse reads do not match!")
}
print(length(fnFs))  # Should equal length(fnRs)

# Create a directory for filtered reads
filtered_path <- file.path(path, "filtered")
if (!dir.exists(filtered_path)) dir.create(filtered_path)

# File paths for filtered reads
filtFs <- file.path(filtered_path, paste0(basename(fnFs)))
filtRs <- file.path(filtered_path, paste0(basename(fnRs)))

# Filter and trim the reads (adjust truncLen based on read length and quality)
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen = c(220, 200),  # Adjust based on read quality
                     trimLeft = c(19, 20),# Triming from both sides
                maxN = 0,               # No Ns allowed
                     maxEE = c(2, 3),        # Maximum expected errors
                     truncQ = 2,             # Truncate reads with low quality
                     rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

# Learn the error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Dereplicate the reads for both forward and reverse
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Extract sample names from file names
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Name the dereplicated objects with sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Run the DADA2 algorithm for paired-end reads
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired-end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=20, maxMismatch=0)

# Create a sequence table
seq_table <- makeSequenceTable(mergers)

# Remove chimeras
seq_table_no_chim <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE)

# Save the sequence table
saveRDS(seq_table_no_chim, file = file.path(path, "Nseqtab.nochim.rds"))

# Assign taxonomy using the SILVA nr99 database
silva_nr99_w_species <- "/mnt/lustre/users/nkambunga/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxa <- assignTaxonomy(seq_table_no_chim, silva_nr99_w_species, multithread = TRUE)

# Save the taxonomy assignment
saveRDS(taxa, file = file.path(path, "Ntaxa_assigned.rds"))

# Optionally, save the taxonomy table in a CSV file
taxa_df <- as.data.frame(taxa)
write.csv(taxa_df, file = file.path(path, "Ntaxa_assigned.csv"))

# Track the number of reads through each step
track <- cbind(out, rowSums(seq_table_no_chim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample_names

# Save the tracking information
saveRDS(track, file = file.path(path, "Ntrack_reads.rds"))
write.csv(track, file = file.path(path, "Ntrack_reads.csv"))

#---------------------------------------------------------

# Load the necessary packages
library(dada2)

# Define the path to your combined fastq files
path <- "/mnt/lustre/users/nkambunga/JRehousing_nodup"

# Define the path to your combined samples
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz$", full.names = TRUE)) # Forawrd reads
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz$", full.names = TRUE)) # Reverse reads

# Check the number of files (R1 and R2 should match)
if (length(fnFs) != length(fnRs)) {
  stop("The number of forward and reverse reads do not match!")
}
print(length(fnFs))  # Should equal length(fnRs)

# Create a directory for filtered reads
filtered_path <- file.path(path, "filtered")
if (!dir.exists(filtered_path)) dir.create(filtered_path)

# File paths for filtered reads
filtFs <- file.path(filtered_path, paste0(basename(fnFs)))
filtRs <- file.path(filtered_path, paste0(basename(fnRs)))

# Filter and trim the reads (adjust truncLen based on read length and quality)
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     truncLen = c(220, 200),  # Adjust based on read quality
                     trimLeft = c(19, 20),# Triming from both sides
                maxN = 0,               # No Ns allowed
                     maxEE = c(2, 3),        # Maximum expected errors
                     truncQ = 2,             # Truncate reads with low quality
                     rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

# Learn the error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Dereplicate the reads for both forward and reverse
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Extract sample names from file names
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Name the dereplicated objects with sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Run the DADA2 algorithm for paired-end reads
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired-end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=20, maxMismatch=0)

# Create a sequence table
seq_table <- makeSequenceTable(mergers)

# Remove chimeras
seq_table_no_chim <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE)

# Save the sequence table
saveRDS(seq_table_no_chim, file = file.path(path, "Jseqtab.nochim.rds"))

# Assign taxonomy using the SILVA nr99 database
silva_nr99_w_species <- "/mnt/lustre/users/nkambunga/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxa <- assignTaxonomy(seq_table_no_chim, silva_nr99_w_species, multithread = TRUE)

# Save the taxonomy assignment
saveRDS(taxa, file = file.path(path, "Jtaxa_assigned.rds"))

# Optionally, save the taxonomy table in a CSV file
taxa_df <- as.data.frame(taxa)
write.csv(taxa_df, file = file.path(path, "Jtaxa_assigned.csv"))

# Track the number of reads through each step
track <- cbind(out, rowSums(seq_table_no_chim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample_names

# Save the tracking information
saveRDS(track, file = file.path(path, "Jtrack_reads.rds"))
write.csv(track, file = file.path(path, "Jtrack_reads.csv"))

                                                           
