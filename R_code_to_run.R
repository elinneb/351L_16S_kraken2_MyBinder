# Load packages
library(phyloseq)
library(microeco)
library(file2meco)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(tidyr)

# Import biom table from kraken2 analysis
ps <- suppressWarnings(import_biom("merged_16S_kraken_ALL.biom"))

# Set new tax table column names
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Clean up column names (only text before 1st "_") of the OTU table
sample_names(ps) <- sub("_.*", "", sample_names(ps))

# Read in metadata table
metadata <- read.table("sample_metadata.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Convert Metadata to a sample_data Object
metadata <- sample_data(metadata)

# Reorder metadata to match phyloseq sample order
metadata <- metadata[sample_names(ps), , drop = FALSE]

# Convert metadata to phyloseq sample_data and assign directly
sample_data(ps) <- sample_data(metadata)

# Filter out chloroplasts and mitochondria
ps <- subset_taxa(ps, Order != "o__Chloroplast" & Family != "f__Mitochondria")

# Force order of sites and months for plotting
ps@sam_data$Month <- factor(ps@sam_data$Month, levels = c("07 (Jul)","08 (Aug)","09 (Sep)", "10 (Oct)", "11 (Nov)", "12 (Dec)", "01 (Jan)", "02 (Feb)"))
ps@sam_data$Site <- factor(ps@sam_data$Site, levels = c("Mission Point","Paradise Point","Fiesta Sunset Beach", "Fanuel Street Park", "Kendall Frost Marsh", "De Anza Cove", "Leisure Lagoon", "South Shores", "Tecolote Creek", "Rose Creek"))


## Create Microtable object (microeco)
meco_dataset <- phyloseq2meco(ps)

############# TO GENERATE COMMUNITY COMPOSITION PLOTS

# Create object that agglomerates to Order
t <- trans_abund$new(dataset = meco_dataset, taxrank = "Order", ntaxa = 25)

g_s <- t$plot_bar(bar_full = TRUE, others_color = "grey70", facet = c("Site"), legend_text_italic = FALSE, clustering = FALSE, x_axis_name = NULL, color_values = RColorBrewer::brewer.pal(10, "Paired")) + theme_classic() + theme(text = element_text(size = 18), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(fill = guide_legend(ncol = 1, title = "Order"))
g_m <- t$plot_bar(bar_full = TRUE, others_color = "grey70", facet = c("Month"), legend_text_italic = FALSE, clustering = FALSE, x_axis_name = NULL, color_values = RColorBrewer::brewer.pal(10, "Paired")) + theme_classic() + theme(text = element_text(size = 18), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + guides(fill = guide_legend(ncol = 1, title = "Order"))

plot_grid(g_s, g_m, ncol = 1, align = 'v', labels = "AUTO")

ggsave("MB_16S_composition_barplot.pdf", width = 24, height = 18, limitsize = FALSE)



############# To screen for pathogen levels across sites

# clean taxonomy
tax <- as.data.frame(tax_table(ps))
tax[] <- lapply(tax, function(x) sub("^[a-z]__", "", x))
tax_table(ps) <- tax_table(as.matrix(tax))

# transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# define a simple pathogen/opportunist list
pathogen_genera <- c(
  "Escherichia-Shigella",
  "Salmonella",
  "Campylobacter",
  "Vibrio",
  "Staphylococcus",
  "Streptococcus",
  "Enterococcus",
  "Klebsiella",
  "Pseudomonas",
  "Legionella",
  "Listeria"
)

# Aggregate to the genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# extract and flag
otu <- as.data.frame(otu_table(ps_genus))
if (!taxa_are_rows(ps_genus)) otu <- t(otu)

tax <- as.data.frame(tax_table(ps_genus))

otu$Genus <- tax$Genus

# sum abundance of flagged genera per sample
path_summary <- otu %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  filter(Genus %in% pathogen_genera) %>%
  summarise(across(where(is.numeric), sum)) %>%
  pivot_longer(cols = everything(),
               names_to = "Sample",
               values_to = "PathogenFraction")

# add metadata
meta_df <- as.data.frame(sample_data(ps_genus))
meta_df$Sample <- rownames(meta_df)

path_summary <- left_join(path_summary, meta_df, by = "Sample")

# plot by site
ggplot(path_summary, aes(x = Site, y = PathogenFraction)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Fraction of Reads in Pathogen-Associated Genera",
    y = "Relative Abundance",
    x = "Site"
  )

ggsave("MB_16S__pathogen_bySite.pdf", width = 10, height = 6, limitsize = FALSE)



############# OPTIONAL
############# TO ANALYZE ALPHA DIVERSITY

plot_richness(ps, x = "Month", measures = c("Shannon", "Simpson")) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = Site), width = 0.2, size = 2) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0)
  )

ggsave("MB_16S__alpha_byMonth.pdf", width = 10, height = 6, limitsize = FALSE)



plot_richness(ps, x = "Site", measures = c("Shannon", "Simpson")) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = Month), width = 0.2, size = 2) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0)
  )

ggsave("MB_16S__alpha_bySite.pdf", width = 10, height = 6, limitsize = FALSE)


