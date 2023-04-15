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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20),
                     maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
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
# Learns the specific error-signature of our dataset
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize the error rates
png("Fwd_Error_Rate_Plots.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()

png("Rev_Error_Rate_Plots.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()


### Dereplication 
# Combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence 
# I.e. if you have 100 identical sequences, it keeps 1 and attaches the number "100" to it
# Dereplication substantially reduces computation time by eliminating redundant comparisons.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names



### Inferring ASVs - Sample Inference Algorithm
# Incorporates the consensus quality profiles and abundances of each unique sequence
# Then figures out if each sequence is of biological origin or spurious
# Better to run on combined samples to resolve low-abundance sequences
# But for computing reasons each group will run one pair - next week will start with full data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]



### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])



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


### Assign taxonomy with latest Silva database, v.138.1
# Compares ASVs to database of 16S sequences to classify taxa
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

write.table(taxa, "Taxonomic_Table.txt", sep = "\t", quote=F)




################
### Phyloseq ###
################

### Install phyloseq and ggplot2 in R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

install.packages("ggplot2")


### Load Phyloseq and ggplot2
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())


##### THIS SECTION IS ADDED BY BECKET TO CREATE RELATIVE ABUNDANCE TABLES #####
##### SKIP TO "Construct sample data.frame" if you want to do normal Phyloseq code ####

### Make matrix table
t_seqtab <- as.data.frame(t(seqtab.nochim))
t_seqtab1 <- t_seqtab
t_seqtab1$abund_sum <- rowSums(t_seqtab1)

taxabun <- merge(taxa, t_seqtab1, by = 'row.names', type = "full", match = "all")

# sort rows by most abundant Phlya
taxabun <- taxabun[with(taxabun, order(abund_sum, decreasing=TRUE)),]

# remove row sums column
taxabun <- taxabun[,-102]

write.csv(taxabun, "Taxa Absolute Abundance Table.csv", row.names = TRUE, quote = F)


### Convert to relative abundance
library(dplyr)

t_seqtab_perc <- t_seqtab %>% mutate(across(where(is.numeric), ~ ./sum(.)))

t_seqtab_perc$abund_sum <- rowSums(t_seqtab_perc)

taxabun_perc <- merge(taxa, t_seqtab_perc, by = 'row.names', type = "full", match = "all")

# sort rows by most abundant Phlya
taxabun_perc <- taxabun_perc[with(taxabun_perc, order(abund_sum, decreasing=TRUE)),]

# remove row sums column
taxabun_perc <- taxabun_perc[,-102]

write.csv(taxabun_perc, "Taxa Relative Abundance Table.csv", row.names = TRUE, quote = F)







### Construct sample data.frame
samples.out <- rownames(seqtab.nochim)
site <- sapply(strsplit(samples.out, "-"), `[`, 1)
dirt <- sapply(strsplit(samples.out, "-"), `[`, 2)
dirt <- sapply(strsplit(dirt, "_"), `[`, 1)
dirt <- substr(dirt,2,2)
sample <- sapply(strsplit(samples.out, "_"), `[`, 2)


samdf <- data.frame(Site=site, Soil_Type=dirt, Sample=sample)
rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps


# Visualize Alpha-diversity by soil type:
pdf("All Samples Alpha Diversity by Soil Type.pdf",width = 8, height = 7)
plot_richness(ps, x="Soil_Type", measures=c("Shannon", "Simpson"), color="Site")
dev.off()

# Visualize Alpha-diversity by site:
pdf("All Samples Alpha Diversity by Site.pdf",width = 8, height = 7)
plot_richness(ps, x="Site", measures=c("Shannon", "Simpson"), color="Soil_Type")
dev.off()



# Beta-Diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

# Beta-diversity by soil type
pdf("All Samples Beta Diversity by Soil Type.pdf",width = 8, height = 7)
plot_ordination(ps.prop, ord.nmds.bray, color="Soil_Type", title="Bray NMDS")
dev.off()


# Beta-diversity by site
pdf("All Samples Beta Diversity by Site.pdf",width = 8, height = 7)
plot_ordination(ps.prop, ord.nmds.bray, color="Site", title="Bray NMDS")
dev.off()


# Bar Plots
top <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:25]
ps.top <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)

# change taxrank to look at different ones
taxrank <- "Genus"

pdf(file=paste("All Samples",taxrank,"Barplot.pdf"),width = 8, height = 7)
plot_bar(ps.top, x="Sample", fill=paste(taxrank)) + geom_bar(position="fill", stat="identity") + facet_wrap(~Soil_Type, scales="free_x")
dev.off()


# #Agglomerate to phylum-level and rename
# library(dplyr)
# ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
# phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
# phyloseq::otu_table(ps_phylum)[1:5, 1:5]
# 
# #Melt and plot
# phyloseq::psmelt(ps_phylum) %>%
#   ggplot(data = ., aes(x = Soil_Type, y = Abundance)) +
#   geom_boxplot(outlier.shape  = NA) +
#   geom_jitter(aes(color = OTU), height = 0, width = .2) +
#   labs(x = "", y = "Abundance\n") +
#   facet_wrap(~ OTU, scales = "free")


### Subset based on site (choose which subset you want to use, don't use both)
# Can change what you choose to subset by (what's in parentheses)

# By 1 site
site <- "C1"
samdf_SS <- subset(samdf, Site == paste(site))

# By 2 sites
site <- "C1"
site2 <- "C2"
samdf_SS <- subset(samdf, Site == paste(site) | Site == paste(site2))



ps_SS <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf_SS), 
               tax_table(taxa))
ps_SS


# Visualize Alpha-diversity by soil type:
pdf(file=paste(site,"Samples Alpha Diversity by Soil Type.pdf"),width = 8, height = 7)
plot_richness(ps_SS, x="Soil_Type", measures=c("Shannon", "Simpson"), color="Site")
dev.off()


# Beta-Diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop_SS <- transform_sample_counts(ps_SS, function(otu) otu/sum(otu))
ord.nmds.bray_SS <- ordinate(ps.prop_SS, method="NMDS", distance="bray")

# Beta-diversity by soil type
pdf(file=paste(site,"Samples Beta Diversity by Soil Type.pdf"),width = 8, height = 7)
plot_ordination(ps.prop_SS, ord.nmds.bray_SS, color="Soil_Type", title="Bray NMDS")
dev.off()

# Beta-diversity by site
pdf(file=paste(site,"Samples Beta Diversity by Site.pdf"),width = 8, height = 7)
plot_ordination(ps.prop_SS, ord.nmds.bray_SS, color="Site", title="Bray NMDS")
dev.off()


# Bar Plots
top <- names(sort(taxa_sums(ps_SS), decreasing=TRUE))[1:50]
ps_SS.top <- transform_sample_counts(ps_SS, function(OTU) OTU/sum(OTU))
ps_SS.top <- prune_taxa(top, ps_SS.top)


taxrank_site <- "Class"

pdf(file=paste(site,"Samples",taxrank_site,"Barplot.pdf"),width = 8, height = 7)
plot_bar(ps_SS.top, x="Sample", fill=paste(taxrank_site)) + geom_bar(position="fill", stat="identity") + facet_wrap(~Soil_Type, scales="free_x")
dev.off()

