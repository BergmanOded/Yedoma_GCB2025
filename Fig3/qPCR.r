# depth profile and qPCR Aug17 Achsa data----
# create new object with specific columns


library(dplyr)
library(ggplot2)
library(cowplot)

qpcr <- read.table("/home/oded/Documents/orit_sivan/articles/oded/qpcr/24.11.16/24.11.16_yedoma_qPCR_for_r.tsv", sep="\t", header=T, stringsAsFactors=TRUE, comment.char="") 
str(qpcr)

`qpcr`$depth <- as.numeric(`qpcr`$depth)

qpcr <- qpcr %>% arrange(depth)


# Custom label formatting function - to format the x axis
# Custom label formatting functions
format_label_6 <- function(breaks) {
  sapply(breaks, function(x) {
    if (x == 0) {
      "0"
    } else {
      paste0(x / 1e6, " × 10^6")
    }
  })
}

format_label_8 <- function(breaks) {
  sapply(breaks, function(x) {
    if (x == 0) {
      "0"
    } else {
      paste0(x / 1e8, " × 10^8")
    }
  })
}

format_label_5 <- function(breaks) {
  sapply(breaks, function(x) {
    if (x == 0) {
      "0"
    } else {
      paste0(formatC(x / 1e5, format = "f", digits = 0), " x 10^5")
    }
  })
}

format_label_7 <- function(breaks) {
  sapply(breaks, function(x) {
    if (x == 0) {
      "0"
    } else {
      paste0(x / 1e7, " × 10^7")
    }
  })
}

#mcrA

# Extract the required data dynamically from the dataframe
depth1 <- qpcr %>% filter(sample == "N32") %>% pull(depth)
depth2 <- qpcr %>% filter(sample == "N49") %>% pull(depth)
mcrA_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(mcrA)
mcrA_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(mcrA)


mcrA <- ggplot(qpcr, aes(x=mcrA, y=depth, group=exp, color=exp))+
  geom_point(size=2)+
  geom_errorbar(aes(xmin=mcrA-mcrA_SE, xmax=mcrA+mcrA_SE), stat = "identity", width=10)+
  geom_path(linewidth=1)+
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000"))+
scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  scale_x_continuous(labels = format_label_6,  # Use the custom formatting function
                     breaks = seq(0, 7e6, by = 1e6),  # Set breaks at every 1e6
                     limits = c(0, 7e6)) +
  labs(y = "Depth [cm]", x = "mcrA [Copies/gr soil]", size = 14) +
  theme_classic()+   
  # Add the dashed line between depth1 and depth2
  geom_segment(aes(x = mcrA_depth1, y = depth1, xend = mcrA_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))
mcrA 

# pmoA
pmoA_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(pmoA)
pmoA_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(pmoA)

pmoA <- ggplot(qpcr, aes(x = pmoA, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = pmoA - pmoA_SE, xmax = pmoA + pmoA_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  scale_x_continuous(labels = format_label_8,
                     breaks = seq(0, 2e8, by = 5e7),
                     limits = c(0, 2e8)) +
  labs(y = NULL, x = "pmoA [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = pmoA_depth1, y = depth1, xend = pmoA_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

pmoA

# hzsB
hzsB_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(hzsB)
hzsB_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(hzsB)

hzsB <- ggplot(qpcr, aes(x = hzsB, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = hzsB - hzsB_SE, xmax = hzsB + hzsB_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  scale_x_continuous(labels = format_label_6,
                     breaks = seq(0, 5e6, by = 1e6),
                     limits = c(0, 5e6)) +
  labs(y = NULL, x = "hzsB [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = hzsB_depth1, y = depth1, xend = hzsB_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

hzsB

# NC10
NC10_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(NC10)
NC10_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(NC10)

NC10 <- ggplot(qpcr, aes(x = NC10, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = NC10 - NC10_SE, xmax = NC10 + NC10_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_5,
                     breaks = seq(0, 190000, by = 50000),
                     limits = c(0, 190000)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = NULL, x = "NC10 [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = NC10_depth1, y = depth1, xend = NC10_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

NC10

# A6_acm
A6_acm_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(A6.acm)
A6_acm_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(A6.acm)

A6_acm <- ggplot(qpcr, aes(x = A6.acm, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = A6.acm - A6.acm_SE, xmax = A6.acm + A6.acm_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_7,
                     breaks = seq(0, 5e7, by = 1e7),
                     limits = c(0, 5e7)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = NULL, x = "Acidimicrobium A6 [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = A6_acm_depth1, y = depth1, xend = A6_acm_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))
A6_acm

# nirK_bact
nirK_bact_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(nirK_b)
nirK_bact_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(nirK_b)

nirK_bact <- ggplot(qpcr, aes(x = nirK_b, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = nirK_b - nirK_b_SE, xmax = nirK_b + nirK_b_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_6,
                     breaks = seq(0, 4e6, by = 1e6),
                     limits = c(0, 4e6)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = "Depth [cm]", x = "nirK [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = nirK_bact_depth1, y = depth1, xend = nirK_bact_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

nirK_bact

# narG
narG_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(narG)
narG_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(narG)

narG <- ggplot(qpcr, aes(x = narG, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = narG - narG_SE, xmax = narG + narG_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_7,
                     breaks = seq(0, 5e7, by = 1e7),
                     limits = c(0, 5e7)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = NULL, x = "narG [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = narG_depth1, y = depth1, xend = narG_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

narG

# amoA_arc
amoA_arc_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(amoA_arc)
amoA_arc_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(amoA_arc)

amoA_arc <- ggplot(qpcr, aes(x = amoA_arc, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = amoA_arc - amoA_arc_SE, xmax = amoA_arc + amoA_arc_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_5,
                     breaks = seq(0, 400000, by = 100000),
                     limits = c(0, 400000)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = NULL, x = "amoA Arcaea [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = amoA_arc_depth1, y = depth1, xend = amoA_arc_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

amoA_arc

# amoA_bact
amoA_bact_depth1 <- qpcr %>% filter(sample == "N32") %>% pull(amoA_bact)
amoA_bact_depth2 <- qpcr %>% filter(sample == "N49") %>% pull(amoA_bact)

amoA_bact <- ggplot(qpcr, aes(x = amoA_bact, y = depth, group = exp, color = exp)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = amoA_bact - amoA_bact_SE, xmax = amoA_bact + amoA_bact_SE), stat = "identity", width = 10) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("SSHE" = "#af58ba", "SSLE" = "#a6ad4b", "WDLE" = "#3c93c2", "WSLE" = "#800000")) +
  scale_x_continuous(labels = format_label_6,
                     breaks = seq(0, 3e6, by = 1e6),
                     limits = c(0, 3e6)) +
  scale_y_continuous(breaks = seq(0, -750, by = -50)) +
  labs(y = NULL, x = "amoA Bacteria [Copies/gr soil]") +
  theme_classic() +
  geom_segment(aes(x = amoA_bact_depth1, y = depth1, xend = amoA_bact_depth2, yend = depth2),
               linetype = 3, color = "gray", linewidth = 0.8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5))

amoA_bact

qpcr_yadoma <- plot_grid(mcrA, pmoA, hzsB, NC10, A6_acm, 
                                   nirK_bact, narG, amoA_bact, amoA_arc, nrow=2, ncol=5, align="h", scale = 0.9)

pdf("/home/oded/Documents/orit_sivan/articles/oded/qpcr/24.11.16/24.11.16_qpcr_yadoma.pdf", width = 12, height = 10, family="ArialMT")
qpcr_yadoma
dev.off()
