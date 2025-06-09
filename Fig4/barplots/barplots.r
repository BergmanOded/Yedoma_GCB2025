library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# all features (to look at overall diversity)

metadata <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/biom/sample-metadata.tsv", sep = "\t", header = T, 
                               comment.char = "", row.names = 1, check.names = FALSE)
str(metadata)

metadata$depth <- as.numeric(metadata$depth)

ASV <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/biom/feature-table.tsv", sep = "\t", header = T, 
                           comment.char = "", row.names = 1, check.names = FALSE)
# Define the desired order of columns
desired_order <- c("8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2",
                   "23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1",
                   "50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6", "345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6")

# Reorder the columns of otu_table according to the desired order
ASV <- ASV[, desired_order]
str(ASV)

ASV <- ASV %>%
  mutate_if(is.integer, as.numeric) 

# taxonomy file (has more ASVs as it was constructed based on additional experiments)

taxonomy<-read.table("/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/taxonomy/biom/taxonomy.tsv", sep = "\t", header = T, comment.char = "", row.names = 1, check.names = FALSE)

# Filter taxonomy file to include only ASVs present in ASV
filtered_taxonomy <- taxonomy[rownames(taxonomy) %in% rownames(ASV), ]

# Reorder taxonomy to match the order of ASVs in ASV
ordered_taxonomy <- filtered_taxonomy[match(rownames(ASV), rownames(filtered_taxonomy)), ]

# Verify the order (optional)
all(rownames(ordered_taxonomy) == rownames(ASV))  # Should return TRUE

# Save the result (optional)
write.table(ordered_taxonomy, "/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/taxonomy/biom/24.11.18_yedoma_filtered_ordered_taxonomy.tsv", sep = "\t", quote = FALSE)

# Split the Taxon column into separate taxonomic ranks
ordered_taxonomy_parsed <- ordered_taxonomy %>%
  mutate(
    kingdom = sub("^d__([^;]*);.*", "\\1", Taxon),
    phylum = sub("^.*p__([^;]*);.*", "\\1", Taxon),
    class = sub("^.*c__([^;]*);.*", "\\1", Taxon),
    order = sub("^.*o__([^;]*);.*", "\\1", Taxon),
    family = sub("^.*f__([^;]*);.*", "\\1", Taxon),
    genus = sub("^.*g__([^;]*);.*", "\\1", Taxon),
    species = sub("^.*s__([^;]*)$", "\\1", Taxon)
  )

# Remove the original Taxon column
ordered_taxonomy_parsed <- ordered_taxonomy_parsed %>% select(-Taxon, -Confidence)

# ASVs + taxonomy
ASV_taxonomy <- cbind(ASV = rownames(ASV), ASV, ordered_taxonomy_parsed)

write.csv(ASV_taxonomy, "/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/taxonomy/biom/24.11.18_ASV_taxonomy.csv", row.names = FALSE)

ASV_taxonomy <- cbind(ASV, ordered_taxonomy_parsed)

ASV_arch<- dplyr::filter(ASV_taxonomy, grepl('Archaea', kingdom))
ASV_bact <- dplyr::filter(ASV_taxonomy, grepl('Bacteria', kingdom))

# now omit the Kingdom column, for final otu_table for phyloseq (by kingdom)

ASV_arch <- ASV_arch[c(1:27)]

ASV_bact <- ASV_bact[c(1:27)]

# separate taxonomy to Bacteria and Archaea 

taxonomy_arch<- dplyr::filter(ordered_taxonomy_parsed, grepl('Archaea', kingdom))
taxonomy_bact <- dplyr::filter(ordered_taxonomy_parsed, grepl('Bacteria', kingdom))

# phyloseq object

otu_table_bact = otu_table(as.matrix(ASV_bact), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_bact <- tax_table(as.matrix(taxonomy_bact))
BP_bact <- phyloseq(otu_table_bact, SAMPLE, taxasums_bact)

otu_table_arch = otu_table(as.matrix(ASV_arch), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_arch <- tax_table(as.matrix(taxonomy_arch))
BP_arch <- phyloseq(otu_table_arch, SAMPLE, taxasums_arch)


# Transform BP to relative abundances
BP.rel.bact <- transform_sample_counts(BP_bact, function(x) x / sum(x))
BP.rel.arch <- transform_sample_counts(BP_arch, function(x) x / sum(x))

# Aggregate at the Phylum level
BP.rel.bact_phylum <- tax_glom(BP.rel.bact, "phylum")
BP.rel.arch_phylum <- tax_glom(BP.rel.arch, "phylum")

# convert phyla below 1% to others and plot

# first bacteria 
# Calculate the total abundance of each phylum across all samples
phylum_abundances_bact <- psmelt(BP.rel.bact_phylum) %>%
  group_by(phylum) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  ungroup()

# Identify low-abundance phyla
low_abundance_phyla <- phylum_abundances_bact %>%
  filter(TotalAbundance < 0.01) %>%
  pull(phylum)

# Get the taxonomy table as a data frame
tax_table_df_bact <- as.data.frame(tax_table(BP.rel.bact))

# Update phylum names for low-abundance phyla
tax_table_df_bact$phylum <- if_else(tax_table_df_bact$phylum %in% low_abundance_phyla, "Others", tax_table_df_bact$phylum)

# Update the taxonomy table in BP.rel
tax_table(BP.rel.bact) <- tax_table(as.matrix(tax_table_df_bact))

# Aggregate at Phylum level again to combine "Others"
BP.rel.bact_phylum_updated <- tax_glom(BP.rel.bact, "phylum")

# to get same colors as other graphs -

getPalette = colorRampPalette(brewer.pal(11, "Paired"))

speciesList_bact_phyl = unique(tax_table(BP.rel.bact_phylum_updated)[,"phylum"])
speciesPalette_bact_phyl = getPalette(length(speciesList_bact_phyl))
names(speciesPalette_bact_phyl) = speciesList_bact_phyl

# Convert phyloseq object to data frame
  bp_data_bact <- psmelt(BP.rel.bact_phylum_updated)
  
  # Define your custom labels and desired order
  sample_order <- c("711", "695", "655", "564", "535", "345", "285", "198", "162", "100", "50",
                    "305", "295", "259", "186", "166", "106", "86", "68", "56", "23",
                    "225", "195", "146", "84", "38", "8")
  
  
  # Map original sample names to custom labels
  bp_data_bact$sample <- factor(bp_data$Sample, levels = c("8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2",
                                                      "23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1",
                                                      "50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6", "345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6"), 
                           labels = sample_order)                                                   
  
  # Create the plot
  barplot_phyl_bact <- ggplot(bp_data_bact, aes(x = sample, y = Abundance*100, fill = phylum)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.02) +
    ylab("Abundance (%)") +
    scale_x_discrete(labels = sample_order) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Convert to percentages
    scale_fill_manual(values = speciesPalette_bact_phyl) +
    theme(axis.title.y = element_text(color = "black", size = 16, vjust = 3),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    coord_flip() 
  
  barplot_phyl_bact
  pdf("/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/figs/24.11.18_barplot_bact_phylum_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot_phyl_bact
  dev.off()
  
  
  # now archaea 
  # Calculate the total abundance of each phylum across all samples
  phylum_abundances_arch <- psmelt(BP.rel.arch_phylum) %>%
    group_by(phylum) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    ungroup()
  
  # Identify low-abundance phyla
  low_abundance_phyla <- phylum_abundances_arch %>%
    filter(TotalAbundance < 0.01) %>%
    pull(phylum)
  
  # Get the taxonomy table as a data frame
  tax_table_df_arch <- as.data.frame(tax_table(BP.rel.arch))
  
  # Update phylum names for low-abundance phyla
  tax_table_df_arch$phylum <- if_else(tax_table_df_arch$phylum %in% low_abundance_phyla, "Others", tax_table_df_arch$phylum)
  
  # Update the taxonomy table in BP.rel
  tax_table(BP.rel.arch) <- tax_table(as.matrix(tax_table_df_arch))
  
  # Aggregate at Phylum level again to combine "Others"
  BP.rel.arch_phylum_updated <- tax_glom(BP.rel.arch, "phylum")
  
  # to get same colors as other graphs -
  
  getPalette = colorRampPalette(brewer.pal(8, "Paired"))
  
  speciesList_arch_phyl = unique(tax_table(BP.rel.arch_phylum_updated)[,"phylum"])
  speciesPalette_arch_phyl = getPalette(length(speciesList_arch_phyl))
  names(speciesPalette_arch_phyl) = speciesList_arch_phyl

  # Convert phyloseq object to data frame
  bp_data_arc <- psmelt(BP.rel.arch_phylum_updated)
  
  # Map original sample names to custom labels
  bp_data_arc$sample <- factor(bp_data_arc$Sample, levels = c("711_BH6", "695_BH6", "655_BH6", "564_BH6", "535_BH6", "345_BH6", "285_BH6", 
                                                              "198_BH6", "162_BH6", "100_BH6", "50_BH6", 
                                                              "305_BH1", "295_BH1", "259_BH1", "186_BH1", "166_BH1", "106_BH1", "86_BH1", 
                                                              "68_BH1", "56_BH1", "23_BH1", 
                                                              "225_BH2", "195_BH2", "146_BH2", "84_BH2", "38_BH2", "8_BH2"), 
                           labels = sample_order)                                                   
  
  # Create the plot
  barplot_phyl_arc <- ggplot(bp_data_arc, aes(x = sample, y = Abundance*100, fill = phylum)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.02) +
    ylab("Abundance (%)") +
    scale_x_discrete(labels = sample_order) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Convert to percentages
    scale_fill_manual(values = speciesPalette_arch_phyl) +
    theme(axis.title.y = element_text(color = "black", size = 16, vjust = 3),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    coord_flip() 
  
  barplot_phyl_arc
  pdf("/home/oded/Documents/orit_sivan/articles/oded/NGS/barplot/24.11.16/figs/24.12.22_barplot_arc_phylum_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot_phyl_arc
  dev.off()
  

  
  
  
  
  
  

  # Aggregate at the class level
  BP.rel.bact_class <- tax_glom(BP.rel.bact, "class")
  BP.rel.arch_class <- tax_glom(BP.rel.arch, "class")
  
  # first bacteria 
  # Calculate the total abundance of each class across all samples
  class_abundances_bact <- psmelt(BP.rel.bact_class) %>%
    group_by(class) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    ungroup()
  
  # Identify low-abundance class
  low_abundance_class <- class_abundances_bact %>%
    filter(TotalAbundance < 0.01) %>%
    pull(class)
  
  # Get the taxonomy table as a data frame
  tax_table_df_bact <- as.data.frame(tax_table(BP.rel.bact))
  
  # Update class names for low-abundance class
  tax_table_df_bact$class <- if_else(tax_table_df_bact$class %in% low_abundance_class, "Others", tax_table_df_bact$class)
  
  # Update the taxonomy table in BP.rel
  tax_table(BP.rel.bact) <- tax_table(as.matrix(tax_table_df_bact))
  
  # Aggregate at class level again to combine "Others"
  BP.rel.bact_class_updated <- tax_glom(BP.rel.bact, "class")
  
  # to get same colors as other graphs -
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  speciesList_bact_class = unique(tax_table(BP.rel.bact_class_updated)[,"class"])
  speciesPalette_bact_class = getPalette(length(speciesList_bact_class))
  names(speciesPalette_bact_class) = speciesList_bact_class
  
  # geom_bar makes the borders within the bars
  
  #all samples, by month
  barplot_class_bact <- plot_bar(BP.rel.bact_class_updated, x= "sample", y = "Abundance", fill="class") +
    ylab("Abundance (proportion)")+ xlab("Sample")+
    scale_x_discrete( labels = c("23 BH1", "56 BH1", "68 BH1", "86 BH1", "106 BH1", "166 BH1", "186 BH1", "259 BH1", "295 BH1", "305 BH1",
                                 "8 BH2", "38 BH2", "84 BH2", "146 BH2", "195 BH2", "225 BH2",
                                 "50 BH6", "100 BH6", "162 BH6", "198 BH6", "285 BH6", "345 BH6", "535 BH6", "564 BH6", "655 BH6", "695 BH6", "711 BH6"))+
    theme(axis.title.y = element_text(color="black", size=16, vjust = 3),
          axis.text.y = element_text(color="black", size=14),
          axis.title.x = element_text(color="black", size=16, vjust = 0.1),
          axis.text.x = element_text(color="black", size=14)) +
    scale_fill_manual(values= speciesPalette_bact_class) +
    geom_bar(stat="identity", position="stack", color="black", linewidth=0.02)  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
  barplot_class_bact
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_bact_class_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot_class_bact
  dev.off()
  
  # get legend 
  
  # Define a function to extract the legend
  # here barplot_ord - omit theme(legend.position = "none")
  
  get_legend <- function(barplot_class_bact) {
    # Extract the legend from a plot and return it as a plot object
    tmp <- ggplot_gtable(ggplot_build(barplot_class_bact))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Use the function to extract the legend from your plot
  legend_only <- get_legend(barplot_class_bact)
  
  library(grid)
  
  # Draw the legend on a new page
  grid.newpage()
  grid.draw(legend_only)
  
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_class_yedome_legend.pdf",  width = 10, height = 6, family="ArialMT")
  grid.draw(legend_only)
  dev.off()
  
  
  # now archaea 
  # Calculate the total abundance of each class across all samples
  class_abundances_arch <- psmelt(BP.rel.arch_class) %>%
    group_by(class) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    ungroup()
  
  # Identify low-abundance class
  low_abundance_class <- class_abundances_arch %>%
    filter(TotalAbundance < 0.01) %>%
    pull(class)
  
  # Get the taxonomy table as a data frame
  tax_table_df_arch <- as.data.frame(tax_table(BP.rel.arch))
  
  # Update class names for low-abundance class
  tax_table_df_arch$class <- if_else(tax_table_df_arch$class %in% low_abundance_class, "Others", tax_table_df_arch$class)
  
  # Update the taxonomy table in BP.rel
  tax_table(BP.rel.arch) <- tax_table(as.matrix(tax_table_df_arch))
  
  # Aggregate at class level again to combine "Others"
  BP.rel.arch_class_updated <- tax_glom(BP.rel.arch, "class")
  
  # to get same colors as other graphs -
  
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  
  speciesList_arch_class = unique(tax_table(BP.rel.arch_class_updated)[,"class"])
  speciesPalette_arch_class = getPalette(length(speciesList_arch_class))
  names(speciesPalette_arch_class) = speciesList_arch_class
  
  #all samples, by month
  barplot_class_arch <- plot_bar(BP.rel.arch_class_updated, x= "sample", y = "Abundance", fill="class") +
    ylab("Abundance (proportion)")+
    scale_x_discrete( labels = c("23 BH1", "56 BH1", "68 BH1", "86 BH1", "106 BH1", "166 BH1", "186 BH1", "259 BH1", "295 BH1", "305 BH1",
                                 "8 BH2", "38 BH2", "84 BH2", "146 BH2", "195 BH2", "225 BH2",
                                 "50 BH6", "100 BH6", "162 BH6", "198 BH6", "285 BH6", "345 BH6", "535 BH6", "564 BH6", "655 BH6", "695 BH6", "711 BH6"))+
    theme(axis.title.y = element_text(color="black", size=16, vjust = 3),
          axis.text.y = element_text(color="black", size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color="black", size=14)) +
    scale_fill_manual(values= speciesPalette_arch_class) +
    geom_bar(stat="identity", position="stack", color="black", linewidth=0.02)  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"))
  barplot_class_arch
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_arch_class_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot_class_arch
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    # and at class level
  
  # Aggregate at the class level
  BP.rel_class <- tax_glom(BP.rel, "class")
  
  # Calculate the total abundance of each class across all samples
  class_abundances <- psmelt(BP.rel_class) %>%
    group_by(class) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    ungroup()
  
  # Identify low-abundance phyla
  low_abundance_phyla <- class_abundances %>%
    filter(TotalAbundance < 0.01) %>%
    pull(class)
  
  # Get the taxonomy table as a data frame
  tax_table_df <- as.data.frame(tax_table(BP.rel))
  
  # Update class names for low-abundance phyla
  tax_table_df$class <- if_else(tax_table_df$class %in% low_abundance_phyla, "Others", tax_table_df$class)
  
  # Update the taxonomy table in BP.rel
  tax_table(BP.rel) <- tax_table(as.matrix(tax_table_df))
  
  
  # Aggregate at class level again to combine "Others"
  BP.rel_class_updated <- tax_glom(BP.rel, "class")
  
  
  write.csv(tax_table_BP500.rel, "/home/oded/Documents/post/articles/A_station/abarpot/month/barplot_tax_table.csv", row.names = FALSE)
  
  
  # to get same colors as other graphs -
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  speciesList = unique(tax_table(BP.rel)[,"class"])
  speciesPalette = getPalette(length(speciesList))
  names(speciesPalette) = speciesList
  
  # geom_bar makes the borders within the bars
  
  #all samples, by month
  barplot_class <- plot_bar(BP.rel_class_updated, x= "sample", y = "Abundance", fill="class") +
    ylab("Abundance (proportion)")+
    scale_x_discrete( labels = c("23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1",
                                 "8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2",
                                 "50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6", "345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6"))+
    theme(axis.title.y = element_text(color="black", size=16, vjust = 3),
          axis.text.y = element_text(color="black", size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color="black", size=14)) +
    scale_fill_manual(values= speciesPalette) +
    geom_bar(stat="identity", position="stack", color="black", linewidth=0.02)  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
  
  barplot_class
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_class_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot
  dev.off()
  
  
  
  
  # Define a function to extract the legend
  # here barplot_ord - omit theme(legend.position = "none")
  
  get_legend <- function(barplot_class) {
    # Extract the legend from a plot and return it as a plot object
    tmp <- ggplot_gtable(ggplot_build(barplot_class))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Use the function to extract the legend from your plot
  legend_only <- get_legend(barplot_class)
  
  library(grid)
  
  # Draw the legend on a new page
  grid.newpage()
  grid.draw(legend_only)
  
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_class_yedome_legend.pdf",  width = 15, height = 10, family="ArialMT")
  grid.draw(legend_only)
  dev.off()
  
  
  
  
  
  
  # and at order level
  
  
  # Aggregate at the order level
  BP.rel_order <- tax_glom(BP.rel, "order")
  
  # Calculate the total abundance of each order across all samples
  order_abundances <- psmelt(BP.rel_order) %>%
    group_by(order) %>%
    summarize(TotalAbundance = sum(Abundance)) %>%
    ungroup()
  
  # Identify low-abundance phyla
  low_abundance_phyla <- order_abundances %>%
    filter(TotalAbundance < 0.01) %>%
    pull(order)
  
  # Get the taxonomy table as a data frame
  tax_table_df <- as.data.frame(tax_table(BP.rel))
  
  # Update order names for low-abundance phyla
  tax_table_df$order <- if_else(tax_table_df$order %in% low_abundance_phyla, "Others", tax_table_df$order)
  
  # Update the taxonomy table in BP.rel
  tax_table(BP.rel) <- tax_table(as.matrix(tax_table_df))
  
  
  # Aggregate at order level again to combine "Others"
  BP.rel_order_updated <- tax_glom(BP.rel, "order")
  
  
  write.csv(tax_table_BP500.rel, "/home/oded/Documents/post/articles/A_station/abarpot/month/barplot_tax_table.csv", row.names = FALSE)
  
  
  # to get same colors as other graphs -
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  speciesList = unique(tax_table(BP.rel)[,"order"])
  speciesPalette = getPalette(length(speciesList))
  names(speciesPalette) = speciesList
  
  # geom_bar makes the borders within the bars
  
  #all samples, by month
  barplot_ord <- plot_bar(BP.rel_order_updated, x= "sample", y = "Abundance", fill="order") +
    ylab("Abundance (proportion)")+
    scale_x_discrete( labels = c("23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1",
                                 "8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2",
                                 "50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6", "345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6"))+
    theme(axis.title.y = element_text(color="black", size=16, vjust = 3),
          axis.text.y = element_text(color="black", size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(color="black", size=14)) +
    scale_fill_manual(values= speciesPalette) +
    geom_bar(stat="identity", position="stack", color="black", linewidth=0.02)  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
   
  
  barplot_ord
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_order_yedome.pdf",  width = 15, height = 10, family="ArialMT")
  barplot_ord
  dev.off()
  
  
  
  # Define a function to extract the legend
  # here barplot_ord - omit theme(legend.position = "none")

  get_legend <- function(barplot_ord) {
    # Extract the legend from a plot and return it as a plot object
    tmp <- ggplot_gtable(ggplot_build(barplot_ord))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Use the function to extract the legend from your plot
  legend_only <- get_legend(barplot_ord)
  
  library(grid)
  
  # Draw the legend on a new page
  grid.newpage()
  grid.draw(legend_only)
  
  
  pdf("/home/oded/Documents/orit_sivan/raw_data/170-170/cluster100/exclude/alaska/4runs/final/rename/barplots/24.03.11_barplot_order_yedome_legend.pdf",  width = 15, height = 10, family="ArialMT")
  grid.draw(legend_only)
  dev.off()
