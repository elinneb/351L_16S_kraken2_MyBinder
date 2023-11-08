install.packages("vegan")
install.packages("dendextend")
install.packages("tidyr")
install.packages("viridis")
install.packages("reshape")

install.packages("BiocManager")

source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("phyloseq")
