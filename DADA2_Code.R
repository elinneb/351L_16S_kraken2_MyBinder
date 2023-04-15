### To run the following code:

### If you want to run 1 line of code:
### Click once on a line of code you want to run, and click the run button above ^

### If you are directed to run multiple lines of code at a time
### Drag to hightlight all the lines you are instructed to run
### Then click the run button above ^

### When STOP sign symbol appears in the upper bar of the Console below, CODE IS RUNNING
### Code is done when the arrow ">" appears in the Console below.


### Load dada2 program
library(dada2); packageVersion("dada2")


# define directory where seq output files are
path <- "~/"
list.files(path)



### read in sequencing files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)


# visualize the quality profiles of the forward reads
png("Fwd_Read_Quality_Profiles.png")
plotQualityProfile(fnFs[1:2])
dev.off()

# visualize the quality profiles of the reverse reads
png("Rev_Read_Quality_Profiles.png")
plotQualityProfile(fnRs[1:2])
dev.off()


### Filter and Trim
# 16S V4 region is ~254bp

# Assign the filenames and path for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim reads, note: 16S V4 region is ~254bp
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), rm.phix=TRUE)
head(out)

write.table(out, "Read_Filter_In_Out.txt", sep = "\t", quote=F)


# visualize the quality profiles of the FILTERED forward reads
png("FiltFwd_Read_Quality_Profiles.png")
plotQualityProfile(filtFs[1:2])
dev.off()

# visualize the quality profiles of the FILTERED reverse reads
png("FiltRev_Read_Quality_Profiles.png")
plotQualityProfile(filtRs[1:2])
dev.off()



### Generate an error model of our data
# Learns the specific error-signature of our dataset for later steps
# You don't have to worry about the specifics of how this works, it's to
# inform the program on where the data are vulnerable to errors and how to
# account for it.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)



### Dereplication 
# Combines all identical sequencing reads into into “unique sequences” 
# with a corresponding “abundance” equal to the number of reads with that unique sequence 
# I.e. if you have 100 identical sequences, it keeps 1 and attaches the number "100" to it
# Dereplication substantially reduces computation time by eliminating redundant comparisons.
# NOTE the number of unique sequences per total sequence reads in red below after code is run
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
# (attaches sample names to new files for the program's reference)
names(derepFs) <- sample.names
names(derepRs) <- sample.names



### Inferring ASVs - Sample Inference Algorithm
# Incorporates the consensus quality profiles and abundances of each unique sequence
# Then figures out if each sequence is of biological origin or spurious
# Better to run on combined samples to resolve low-abundance sequences
# But for computing reasons each group will run one pair
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]



### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=TRUE)



### Construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods
seqtab <- makeSequenceTable(mergers)
## DIM function - returns the number of elements (ASVs) from seqtab in a one-dimensional array
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



### Remove chimeras
# Sees if there are any lower-abundance sequences that can be made by mixing left and right portions of more-abundant sequences
# If high, might be from ambiguous bases in primer sequences - revisit removal of primers
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# DIM function - returns the number of elements (ASVs) from seqtab.nochim in a one-dimensional array
dim(seqtab.nochim)
# Get percentage of non-chimeric reads
sum(seqtab.nochim)/sum(seqtab)


### Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100,1))
colnames(track) <- c("input", "filtered", "dadaF", "dadaR", "merged", "nonchim", "final_perc_reads_retained")
rownames(track) <- sample.names

head(track)
write.table(track, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)


######################################################
#### WORKS UNTIL HERE, TROUBLESHOOTING NEXT STEPS ####
######################################################


# Save last file as output so we can re-open R and clear memory
write.table(seqtab.nochim, "seqtab.nochim.tsv", quote=FALSE, sep="\t", col.names=NA)


### Go to JuptyerLab tab and right-click to download the seqtab.nochim file locally
### CLOSE BINDER DOWN. Restart binder link. 
### Drag-drop the seqtab.nochim file into the file list in JupyterLab launch page.
### Reopen R-studio. WHEN YOU REOPEN, RESUME HERE!!!

seqtab.nochim <- as.matrix(read.delim("seqtab.nochim.tsv", sep = "\t", header = T, row.names=1, stringsAsFactors=FALSE))

### Load dada2 program
library(dada2); packageVersion("dada2")

### Assign taxonomy with latest Silva database, v.138.1
# Compares ASVs to database of 16S sequences to classify taxa
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

write.table(taxa, "Taxonomic_Table.txt", sep = "\t", quote=F)



