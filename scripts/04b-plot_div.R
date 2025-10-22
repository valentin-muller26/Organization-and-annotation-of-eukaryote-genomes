library(reshape2)
library(tidyverse)
library(data.table)

# get data from parameter
data <- "HiFiasm_Lu1_primary.fa.mod.out.landscape.Div.Rname.tab"

rep_table <- fread(data, header = FALSE, sep = "\t")
rep_table %>% head()

colnames(rep_table) <- c("Rname", "Rclass", "Rfam", 1:50)
rep_table <- rep_table %>% filter(Rfam != "unknown")
rep_table$fam <- paste(rep_table$Rclass, rep_table$Rfam, sep = "/")

table(rep_table$fam)
# How many elements are there in each Superfamily?

rep_table.m <- melt(rep_table)

rep_table.m <- rep_table.m[-c(which(rep_table.m$variable == 1)), ]

# Arrange the data so that they are in the following order:
# NOTE: Verify that all superfamilies in your dataset are included
rep_table.m$fam <- factor(rep_table.m$fam, levels = c(
  "LTR/Copia", "LTR/Gypsy", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron",
  "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM"
))

rep_table.m$distance <- as.numeric(rep_table.m$variable) / 100

# Calculate TE age using T = K / 2r formula
# r = 8.22 × 10^-9 substitutions per synonymous site per year
rep_table.m$age <- rep_table.m$distance / (2 * 8.22e-9) / 1e6  # Age in million years

# Remove helitrons (EDTA annotation issues)
rep_table.m <- rep_table.m %>% filter(fam != "DNA/Helitron")

# Create output directory
dir.create("Plots", showWarnings = FALSE)

# Plot by distance (divergence)
ggplot(rep_table.m, aes(fill = fam, x = distance, weight = value / 1000000)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Divergence from consensus") +
  ylab("Sequence (Mbp)") +
  ggtitle("TE Landscape by Divergence") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

ggsave(filename = "Plots/TE_landscape_distance.pdf", width = 10, height = 5, useDingbats = F)

# Plot by age (million years)
ggplot(rep_table.m, aes(fill = fam, x = age, weight = value / 1000000)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Age (Million years)") +
  ylab("Sequence (Mbp)") +
  ggtitle("TE Landscape by Age") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

ggsave(filename = "Plots/TE_landscape_age.pdf", width = 10, height = 5, useDingbats = F)

# Summary statistics for Copia vs Gypsy
copia_gypsy <- rep_table.m %>% 
  filter(fam %in% c("LTR/Copia", "LTR/Gypsy")) %>%
  group_by(fam) %>%
  summarise(
    total_mbp = sum(value) / 1e6,
    mean_age = weighted.mean(age, value),
    median_age = median(rep(age, value))
  )

print("Copia vs Gypsy comparison:")
print(copia_gypsy)