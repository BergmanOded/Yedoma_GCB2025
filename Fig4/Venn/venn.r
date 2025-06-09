#venn diagram spring-clusters

library(tidyverse)
library(VennDiagram)

metadata <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/sample-metadata.tsv", sep = "\t", header = T, 
                       comment.char = "", row.names = 1, check.names = FALSE)
str(metadata)

ASV <- read.table(file = "/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/feature-table.tsv", sep = "\t", header = T, 
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

taxonomy<-read.table("/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/taxonomy.tsv", sep = "\t", header = T, comment.char = "", row.names = 1, check.names = FALSE)

# Filter taxonomy file to include only ASVs present in ASV
filtered_taxonomy <- taxonomy[rownames(taxonomy) %in% rownames(ASV), ]

# Reorder taxonomy to match the order of ASVs in ASV
ordered_taxonomy <- filtered_taxonomy[match(rownames(ASV), rownames(filtered_taxonomy)), ]

# Verify the order (optional)
all(rownames(ordered_taxonomy) == rownames(ASV))  # Should return TRUE

# Save the result (optional)
write.table(ordered_taxonomy, "/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/24.11.18_yedoma_filtered_ordered_rarefied_taxonomy.tsv", sep = "\t", quote = FALSE)

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

write.csv(ASV_taxonomy, "/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/24.11.18_ASV_taxonomy_rarefied.csv", row.names = FALSE)

ASV_taxonomy <- cbind(ASV, ordered_taxonomy_parsed)

ASV_arch<- dplyr::filter(ASV_taxonomy, grepl('Archaea', kingdom))
ASV_bact <- dplyr::filter(ASV_taxonomy, grepl('Bacteria', kingdom))

ASV_arch <- ASV_arch[c(1:27)]
ASV_bact <- ASV_bact[c(1:27)]

#bacteria
# Group samples by environment

feature_table_SSHE <- ASV_bact %>%
  select(contains(c("8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2"))) %>% # 'select' function to choose the columns you want 
  mutate(SSHE = rowSums(.))

feature_table_SSLE <- ASV_bact %>%
  select(contains(c("23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1"))) %>% # 'select' function to choose the columns you want 
  mutate(SSLE = rowSums(.))

feature_table_WSLE <- ASV_bact %>%
  select(contains(c("50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6"))) %>% # 'select' function to choose the columns you want 
  mutate(WSLE = rowSums(.))

feature_table_WDLE <- ASV_bact %>%
  select(contains(c("345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6"))) %>% # 'select' function to choose the columns you want 
  mutate(WDLE = rowSums(.))



feature_table_cbind <- cbind(feature_table_SSHE, feature_table_SSLE, feature_table_WSLE, feature_table_WDLE)
str(feature_table_cbind)

feature_table_cluster <- feature_table_cbind %>%
  select("SSHE", "SSLE", "WSLE", "WDLE")
str(feature_table_cluster)

write.csv(feature_table_cluster,"/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/rarefied_feature_table_emv_bact.csv", row.names = TRUE)


SSHE <- rownames(feature_table_cluster)[feature_table_cluster[,"SSHE"] > 0]
SSLE <- rownames(feature_table_cluster)[feature_table_cluster[,"SSLE"] > 0]
WSLE <- rownames(feature_table_cluster)[feature_table_cluster[,"WSLE"] > 0]
WDLE <- rownames(feature_table_cluster)[feature_table_cluster[,"WDLE"] > 0]

venn_SSHE_SSLE <- venn.diagram(
  x = list('SSHE'=SSHE, 'SSLE'=SSLE),
  fill = c("#af58ba", "#a6ad4b"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
  height = 8000, width = 160000)

grid.draw(venn_SSHE_SSLE)


venn_SSLE_WSLE <- venn.diagram(
  x = list('SSLE'=SSLE, 'WSLE'=WSLE),
  fill = c("#a6ad4b", "#800000"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
  height = 8000, width = 160000)

grid.draw(venn_SSLE_WSLE)


venn_WSLE_WDLE <- venn.diagram(
  x = list('WSLE'=WSLE, 'WDLE'=WDLE),
  fill = c("#800000", "#3c93c2"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
height = 8000, width = 160000)

grid.draw(venn_WSLE_WDLE)

library(cowplot)

alpha_yadoma <- plot_grid(venn_SSHE_SSLE, venn_SSLE_WSLE,venn_WSLE_WDLE , nrow=1, ncol=3, align="h", scale = 0.9)

pdf("/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/24.11.18_rarefied_venn_bact_yedoma.pdf",  width = 12, height = 3, family="ArialMT")
alpha_yadoma
dev.off()

# and archaea

feature_table_SSHE_arch <- ASV_arch %>%
  select(contains(c("8_BH2", "38_BH2", "84_BH2", "146_BH2", "195_BH2", "225_BH2"))) %>% # 'select' function to choose the columns you want 
  mutate(SSHE = rowSums(.))

feature_table_SSLE_arch <- ASV_arch %>%
  select(contains(c("23_BH1", "56_BH1", "68_BH1", "86_BH1", "106_BH1", "166_BH1", "186_BH1", "259_BH1", "295_BH1", "305_BH1"))) %>% # 'select' function to choose the columns you want 
  mutate(SSLE = rowSums(.))

feature_table_WSLE_arch <- ASV_arch %>%
  select(contains(c("50_BH6", "100_BH6", "162_BH6", "198_BH6", "285_BH6"))) %>% # 'select' function to choose the columns you want 
  mutate(WSLE = rowSums(.))

feature_table_WDLE_arch <- ASV_arch %>%
  select(contains(c("345_BH6", "535_BH6", "564_BH6", "655_BH6", "695_BH6", "711_BH6"))) %>% # 'select' function to choose the columns you want 
  mutate(WDLE = rowSums(.))



feature_table_cbind_arch <- cbind(feature_table_SSHE_arch, feature_table_SSLE_arch, feature_table_WSLE_arch, feature_table_WDLE_arch)
str(feature_table_cbind_arch)

feature_table_cluster_arch <- feature_table_cbind_arch %>%
  select("SSHE", "SSLE", "WSLE", "WDLE")
str(feature_table_cluster_arch)

write.csv(feature_table_cluster_arch,"/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/rarefied_feature_table_emv_arch.csv", row.names = TRUE)


SSHE_arch <- rownames(feature_table_cluster_arch)[feature_table_cluster_arch[,"SSHE"] > 0]
SSLE_arch <- rownames(feature_table_cluster_arch)[feature_table_cluster_arch[,"SSLE"] > 0]
WSLE_arch <- rownames(feature_table_cluster_arch)[feature_table_cluster_arch[,"WSLE"] > 0]
WDLE_arch <- rownames(feature_table_cluster_arch)[feature_table_cluster_arch[,"WDLE"] > 0]

venn_SSHE_SSLE_arch <- venn.diagram(
  x = list('SSHE'=SSHE_arch, 'SSLE'=SSLE_arch),
  fill = c("#af58ba", "#a6ad4b"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
  height = 8000, width = 160000)

grid.draw(venn_SSHE_SSLE_arch)


venn_SSLE_WSLE_arch <- venn.diagram(
  x = list('SSLE'=SSLE_arch, 'WSLE'=WSLE_arch),
  fill = c("#a6ad4b", "#800000"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
  height = 8000, width = 160000)

grid.draw(venn_SSLE_WSLE_arch)


venn_WSLE_WDLE_arch <- venn.diagram(
  x = list('WSLE'=WSLE_arch, 'WDLE'=WDLE_arch),
  fill = c("#800000", "#3c93c2"),
  alpha = c(0.7, 0.7), cex= 1.5, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
  height = 8000, width = 160000)

grid.draw(venn_WSLE_WDLE_arch)

library(cowplot)

alpha_yadoma_arch <- plot_grid(venn_SSHE_SSLE_arch, venn_SSLE_WSLE_arch,venn_WSLE_WDLE_arch , nrow=1, ncol=3, align="h", scale = 0.9)

pdf("/home/oded/Documents/orit_sivan/articles/oded/NGS/venn_diagram/24.11.16/24.11.18_rarefied_venn_arch_yedoma.pdf",  width = 12, height = 3, family="ArialMT")
alpha_yadoma_arch
dev.off()
