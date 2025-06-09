# from QIIME2 to r - see qiime2_r_tutorial_script (latest date...)

library(tidyverse)
library(vegan)
library(BiocManager)
library(phyloseq)
library(viridis)
library(pheatmap)
otu <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/heatmap/methanotrophs/feature-table.tsv", sep = ",", header = T, 
                  comment.char = "", check.names = FALSE)

taxonomy <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/heatmap/methanotrophs/taxonomy.tsv", sep = ",", header = T, 
                  comment.char = "", check.names = FALSE)

str(otu)
otu <- otu %>%
  mutate_if(is.integer, as.numeric)

otu_desired_order <- c("8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2",
                       "23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1",
                      "50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6", 
                      "345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6")
otu_rd <- otu[, otu_desired_order]
str(otu_rd)



metadata <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/heatmap/methanotrophs/sample-metadata.tsv", sep = "\t", header = T, 
                       comment.char = "", row.names = 1, check.names = FALSE)

OTU = otu_table(as.matrix(otu_rd), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums <- tax_table(as.matrix(taxonomy))

# merge the data (this is the phyloseq object)
hm_methnotrophs <- phyloseq(OTU, SAMPLE, taxasums)

# Extract the abundance matrix from the phyloseq object
otu_matrix <- as.matrix(otu_table(hm_methnotrophs))




# Optional: Prepare a taxonomy annotation dataframe for the rows (OTUs)
taxa_annotation <- as.data.frame(tax_table(hm_methnotrophs))

# Extract the Genus column from taxa_annotation and create a new data frame for row annotations
genus_annotation <- data.frame(Classification = taxa_annotation$genus)

# Assuming 'Genus' is the column with genus names in your taxa_annotation
rownames(otu_matrix) <- taxa_annotation$genus

# Create an annotation data frame for the Order
combined_annotation <- data.frame(
  Phylum = taxa_annotation$phylum,
  Order = taxa_annotation$order
)
  
# Find duplicate 'genus' names and add a suffix to make them unique
make_unique <- function(names) {
  counts <- table(names)
  duplicated_names <- names(counts[counts > 1])
  for (name in duplicated_names) {
    idx <- which(names == name)
    # Add suffix to duplicates except the first occurrence
    names[idx[-1]] <- paste0(names[idx[-1]], "_", seq_along(idx[-1]))
  }
  names
}

# Apply the function to the 'Genus' column of your taxa_annotation
unique_genus_names <- make_unique(taxa_annotation$genus)
# Now set the row names of otu_matrix and order_annotation using these unique names
rownames(otu_matrix) <- unique_genus_names
rownames(combined_annotation) <- unique_genus_names

# remove _ from sample names
colnames(otu_matrix) <- gsub("_", " ", colnames(otu_matrix))

# Create the heatmap

heatmap_motanotroph <-  pheatmap(otu_matrix,
         annotation_row = combined_annotation,
         scale = "row", # This scales each row (OTU) to have zero mean and unit variance
         clustering_distance_rows = "euclidean",
         cluster_cols=FALSE,
         treeheight_row = 20,
         cellheight = 10,
         cellwidth = 10,
         clustering_method = "complete",
         color = viridis::viridis(100), # Using viridis color scale
         show_rownames = T, # Set to F if there are too many rows
         show_colnames = T # Set to F if there are too many columns
)

heatmap_motanotroph
        pdf("/home/oded/Documents/orit_sivan/articles/oded/NGS/heatmap/methanotrophs/24.11.18_heatmap_methanotroph_yedoma.pdf",  width = 10, height = 4, family="ArialMT")
        heatmap_motanotroph
        dev.off()
