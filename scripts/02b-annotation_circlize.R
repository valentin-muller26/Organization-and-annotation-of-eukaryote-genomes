# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)

# Create plot directory if it doesn't exist
plot_dir <- "plots/circlizeTE"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Load the TE annotation GFF3 file
gff_file <- "HiFiasm_Lu1_primary.fa.mod.EDTA.TEanno.gff3"
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Check the superfamilies present in the GFF3 file, and their counts
gff_data$V3 %>% table()

# custom ideogram data
custom_ideogram <- read.table("HiFiasm_Lu1_primary.fa.fai", header = FALSE, stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = T), ]

# Select only the first 15 longest scaffolds
custom_ideogram <- custom_ideogram[1:15, ]

# Function to filter GFF3 data based on Superfamily
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
  filtered_data <- gff_data[gff_data$V3 == superfamily, ] %>%
    as.data.frame() %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# ============================================================================
# CIRCOS PLOT 1: TE Density (Gypsy and Copia)
# ============================================================================

pdf(file.path(plot_dir, "02-TE_density.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

circos.genomicDensity(filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#2E8B57", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#DC143C", track.height = 0.07, window.size = 1e5)

lgd <- Legend(
  title = "Superfamily", 
  at = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon"),
  legend_gp = gpar(fill = c("#2E8B57", "#DC143C"))
)
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("Plot saved to:", file.path(plot_dir, "02-TE_density.pdf"), "\n")

# ============================================================================
# CIRCOS PLOT 2: All Top Superfamilies
# ============================================================================

superfamily_counts <- gff_data$V3 %>% table() %>% sort(decreasing = TRUE)
top_superfamilies <- names(superfamily_counts)[1:6]

superfamily_colors <- c(
  "#2E8B57",  # SeaGreen - Gypsy
  "#DC143C",  # Crimson - Copia
  "#1E90FF",  # DodgerBlue
  "#FF8C00",  # DarkOrange
  "#9370DB",  # MediumPurple
  "#20B2AA"   # LightSeaGreen
)

pdf(file.path(plot_dir, "02-TE_density_all_superfamilies.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

for(i in seq_along(top_superfamilies)) {
  circos.genomicDensity(filter_superfamily(gff_data, top_superfamilies[i], custom_ideogram), 
                        count_by = "number", col = superfamily_colors[i], 
                        track.height = 0.07, window.size = 1e5)
}

lgd_all <- Legend(
  title = "TE Superfamilies", 
  at = top_superfamilies,
  legend_gp = gpar(fill = superfamily_colors[1:length(top_superfamilies)]),
  ncol = 2
)
draw(lgd_all, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("All superfamilies plot saved to:", file.path(plot_dir, "02-TE_density_all_superfamilies.pdf"), "\n")

# ============================================================================
# Load TEsorter classification and merge with GFF3
# ============================================================================

copia_clades <- read_tsv("Copia.cls.tsv", comment = "", show_col_types = FALSE)
gypsy_clades <- read_tsv("Gypsy.cls.tsv", comment = "", show_col_types = FALSE)

cat("\nCopia clades loaded:", nrow(copia_clades), "rows\n")
cat("Gypsy clades loaded:", nrow(gypsy_clades), "rows\n")

all_clades <- bind_rows(copia_clades, gypsy_clades)
all_clades$TE_ID <- str_extract(all_clades[[1]], "TE_[0-9]+")
gff_data$TE_ID <- str_extract(gff_data$V9, "(?<=Name=)TE_[0-9]+")

gff_enriched <- gff_data %>%
  left_join(all_clades %>% select(TE_ID, Clade, Complete), by = "TE_ID")

cat("\nClade distribution:\n")
print(table(gff_enriched$Clade, useNA = "ifany"))

# Function to filter by clade
filter_clade <- function(gff_enriched, clade_name, custom_ideogram) {
  filtered_data <- gff_enriched %>%
    filter(Clade == clade_name) %>%
    filter(!is.na(Clade)) %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V7) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# ============================================================================
# CIRCOS PLOT 3: Athila and CRM (Centromeric clades)
# ============================================================================

athila_count <- sum(gff_enriched$Clade == "Athila", na.rm = TRUE)
crm_count <- sum(gff_enriched$Clade == "CRM", na.rm = TRUE)

cat("\nAthila elements found:", athila_count, "\n")
cat("CRM elements found:", crm_count, "\n")

if(athila_count > 0 || crm_count > 0) {
  pdf(file.path(plot_dir, "03-Athila_and_CRM.pdf"), width = 10, height = 10)
  gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
  circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
  
  circos.genomicInitialize(custom_ideogram)
  
  if(athila_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "Athila", custom_ideogram), 
                          count_by = "number", col = "#FF1493", 
                          track.height = 0.07, window.size = 1e5)
  }
  
  if(crm_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "CRM", custom_ideogram), 
                          count_by = "number", col = "#4B0082", 
                          track.height = 0.07, window.size = 1e5)
  }
  
  present_clades <- c()
  present_colors <- c()
  
  if(athila_count > 0) {
    present_clades <- c(present_clades, "Athila")
    present_colors <- c(present_colors, "#FF1493")
  }
  if(crm_count > 0) {
    present_clades <- c(present_clades, "CRM")
    present_colors <- c(present_colors, "#4B0082")
  }
  
  lgd_centro <- Legend(
    title = "Centromeric TE Clades", 
    at = present_clades,
    legend_gp = gpar(fill = present_colors)
  )
  draw(lgd_centro, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))
  
  circos.clear()
  dev.off()
  
  cat("Centromeric clades plot saved to:", file.path(plot_dir, "03-Athila_and_CRM.pdf"), "\n")
} else {
  cat("\nWarning: No Athila or CRM elements found.\n")
}

# ============================================================================
# Prepare clade summaries for barplots
# ============================================================================

copia_summary <- copia_clades %>%
  group_by(Clade) %>%
  summarise(
    Count = n(),
    Complete_elements = sum(Complete == "yes", na.rm = TRUE),
    Incomplete_elements = sum(Complete == "no", na.rm = TRUE)
  ) %>%
  arrange(desc(Count)) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

gypsy_summary <- gypsy_clades %>%
  group_by(Clade) %>%
  summarise(
    Count = n(),
    Complete_elements = sum(Complete == "yes", na.rm = TRUE),
    Incomplete_elements = sum(Complete == "no", na.rm = TRUE)
  ) %>%
  arrange(desc(Count)) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

# ============================================================================
# BARPLOT 1: Abundance of Copia and Gypsy Clades
# ============================================================================

plot_data_separated <- bind_rows(
  copia_summary %>% mutate(Superfamily_type = "Copia"),
  gypsy_summary %>% mutate(Superfamily_type = "Gypsy")
) %>%
  arrange(Superfamily_type, desc(Count))

p_separated <- ggplot(plot_data_separated, 
                      aes(x = reorder(Clade, -Count), y = Count, fill = Superfamily_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Copia" = "#DC143C", "Gypsy" = "#2E8B57")) +
  facet_wrap(~Superfamily_type, scales = "free_x", ncol = 2) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
  ) +
  labs(
    title = "Abundance of Copia and Gypsy Clades",
    x = "Clade",
    y = "Count"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file.path(plot_dir, "06-Copia_Gypsy_clades_abundance.png"), 
       p_separated, width = 14, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "06-Copia_Gypsy_clades_abundance.pdf"), 
       p_separated, width = 14, height = 6)

cat("\nAbundance plot saved to:", file.path(plot_dir, "06-Copia_Gypsy_clades_abundance.pdf"), "\n")

# ============================================================================
# BARPLOT 2: Complete vs Incomplete Elements - Copia and Gypsy
# ============================================================================

plot_data_separated_long <- bind_rows(
  copia_summary %>% mutate(Superfamily_type = "Copia"),
  gypsy_summary %>% mutate(Superfamily_type = "Gypsy")
) %>%
  select(Clade, Superfamily_type, Count, Complete_elements, Incomplete_elements) %>%
  pivot_longer(cols = c(Complete_elements, Incomplete_elements),
               names_to = "Status",
               values_to = "Element_Count") %>%
  mutate(Status = str_replace(Status, "_elements", "")) %>%
  arrange(Superfamily_type, desc(Count))

plot_data_separated_labels <- bind_rows(
  copia_summary %>% mutate(Superfamily_type = "Copia"),
  gypsy_summary %>% mutate(Superfamily_type = "Gypsy")
) %>%
  mutate(label = paste0("C:", Complete_elements, " I:", Incomplete_elements))

p_separated_stacked <- ggplot(plot_data_separated_long, 
                              aes(x = reorder(Clade, -Count), y = Element_Count, fill = Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(data = plot_data_separated_labels,
            aes(x = reorder(Clade, -Count), y = Count, label = label),
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Complete" = "#228B22", "Incomplete" = "#FFD700"),
                    labels = c("Complete", "Incomplete")) +
  facet_wrap(~Superfamily_type, scales = "free_x", ncol = 2) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  ) +
  labs(
    title = "Complete vs Incomplete Elements: Copia and Gypsy Clades",
    x = "Clade",
    y = "Count",
    fill = "Element Status"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave(file.path(plot_dir, "07-Copia_Gypsy_Complete_Incomplete.png"), 
       p_separated_stacked, width = 16, height = 7, dpi = 300)
ggsave(file.path(plot_dir, "07-Copia_Gypsy_Complete_Incomplete.pdf"), 
       p_separated_stacked, width = 16, height = 7)

cat("Complete vs Incomplete plot saved to:", file.path(plot_dir, "07-Copia_Gypsy_Complete_Incomplete.pdf"), "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n\n")

cat("Generated plots:\n")
cat("  - 02-TE_density.pdf (Circos: Gypsy & Copia density)\n")
cat("  - 02-TE_density_all_superfamilies.pdf (Circos: Top 6 superfamilies)\n")
cat("  - 03-Athila_and_CRM.pdf (Circos: Centromeric clades)\n")
cat("  - 06-Copia_Gypsy_clades_abundance.pdf (Barplot: Clade abundance)\n")
cat("  - 07-Copia_Gypsy_Complete_Incomplete.pdf (Barplot: Complete vs Incomplete)\n")

cat("\nTotal Copia clades:", nrow(copia_summary), "\n")
cat("Total Gypsy clades:", nrow(gypsy_summary), "\n")