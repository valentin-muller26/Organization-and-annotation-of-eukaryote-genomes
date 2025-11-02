# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Create plot directory if it doesn't exist
plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
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
sum(custom_ideogram$end[1:20])

# Select only the first 20 longest scaffolds
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

# Save to plot folder
pdf(file.path(plot_dir, "02-TE_density.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

# Initialize the circos plot with the custom ideogram
circos.genomicInitialize(custom_ideogram)

# Plot te density
circos.genomicDensity(filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#2E8B57", track.height = 0.07, window.size = 1e5) # SeaGreen
circos.genomicDensity(filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#DC143C", track.height = 0.07, window.size = 1e5) # Crimson

lgd <- Legend(
  title = "Superfamily", 
  at = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon"),
  legend_gp = gpar(fill = c("#2E8B57", "#DC143C"))
)
# Center the legend in the middle of the plot
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("Plot saved to:", file.path(plot_dir, "02-TE_density.pdf"), "\n")

# Now plot all your most abundant TE superfamilies in one plot

# Get the most abundant superfamilies
superfamily_counts <- gff_data$V3 %>% table() %>% sort(decreasing = TRUE)
top_superfamilies <- names(superfamily_counts)[1:6]  # Top 6 superfamilies

# Define brighter, more sympathetic colors for each superfamily
superfamily_colors <- c(
  "#2E8B57",  # SeaGreen - Gypsy
  "#DC143C",  # Crimson - Copia
  "#1E90FF",  # DodgerBlue
  "#FF8C00",  # DarkOrange
  "#9370DB",  # MediumPurple
  "#20B2AA"   # LightSeaGreen
)

# Save all superfamilies plot to plot folder
pdf(file.path(plot_dir, "02-TE_density_all_superfamilies.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

# Plot each superfamily
for(i in seq_along(top_superfamilies)) {
  circos.genomicDensity(filter_superfamily(gff_data, top_superfamilies[i], custom_ideogram), 
                        count_by = "number", col = superfamily_colors[i], 
                        track.height = 0.07, window.size = 1e5)
}

# Add centered legend using ComplexHeatmap for better positioning
lgd_all <- Legend(
  title = "TE Superfamilies", 
  at = top_superfamilies,
  legend_gp = gpar(fill = superfamily_colors[1:length(top_superfamilies)]),
  ncol = 2  # Two columns for better layout
)
draw(lgd_all, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("All superfamilies plot saved to:", file.path(plot_dir, "02-TE_density_all_superfamilies.pdf"), "\n")

# Alternative version with base R legend (if you prefer)
pdf(file.path(plot_dir, "02-TE_density_all_superfamilies_alt.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

# Plot each superfamily
for(i in seq_along(top_superfamilies)) {
  circos.genomicDensity(filter_superfamily(gff_data, top_superfamilies[i], custom_ideogram), 
                        count_by = "number", col = superfamily_colors[i], 
                        track.height = 0.07, window.size = 1e5)
}

# Add centered legend using base R
legend(x = 0.5, y = 0.5, 
       legend = top_superfamilies,
       fill = superfamily_colors[1:length(top_superfamilies)],
       title = "TE Superfamilies",
       bg = "white",
       xpd = TRUE,  # Allow plotting outside plot area
       bty = "n",   # No box around legend
       ncol = 2)    # Two columns

circos.clear()
dev.off()

cat("Alternative plot saved to:", file.path(plot_dir, "02-TE_density_all_superfamilies_alt.pdf"), "\n")

# ============================================================================
# Plot the distribution of Athila and CRM clades
# ============================================================================

# Load TEsorter classification files using read_tsv (works better than read.table)
copia_clades <- read_tsv("Copia.cls.tsv", comment = "", show_col_types = FALSE)
gypsy_clades <- read_tsv("Gypsy.cls.tsv", comment = "", show_col_types = FALSE)

cat("\nCopia clades loaded:", nrow(copia_clades), "rows\n")
cat("Gypsy clades loaded:", nrow(gypsy_clades), "rows\n")

# Combine both classifications
all_clades <- bind_rows(copia_clades, gypsy_clades)

# Extract TE_ID from the first column (which contains TE_00000012_INT#LTR/...)
# The column name is `#TE` which R reads as the first column
all_clades$TE_ID <- str_extract(all_clades[[1]], "TE_[0-9]+")

# Extract TE_ID from GFF3 annotation (from Name= field in column V9)
gff_data$TE_ID <- str_extract(gff_data$V9, "(?<=Name=)TE_[0-9]+")

# Verification
cat("\nTE_IDs extracted from clades:", sum(!is.na(all_clades$TE_ID)), "out of", nrow(all_clades), "\n")
cat("TE_IDs extracted from GFF3:", sum(!is.na(gff_data$TE_ID)), "out of", nrow(gff_data), "\n")

# Check for common IDs
common_ids <- intersect(all_clades$TE_ID[!is.na(all_clades$TE_ID)], 
                        gff_data$TE_ID[!is.na(gff_data$TE_ID)])
cat("Number of matching TE_IDs between files:", length(common_ids), "\n")

# Merge GFF3 with clade information
gff_enriched <- gff_data %>%
  left_join(all_clades %>% select(TE_ID, Clade, Complete), by = "TE_ID")

# Check which clades are present and their counts
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

# Check if Athila and CRM are present
athila_count <- sum(gff_enriched$Clade == "Athila", na.rm = TRUE)
crm_count <- sum(gff_enriched$Clade == "CRM", na.rm = TRUE)

cat("\nAthila elements found:", athila_count, "\n")
cat("CRM elements found:", crm_count, "\n")

if(athila_count > 0 || crm_count > 0) {
  pdf(file.path(plot_dir, "03-Athila_and_CRM.pdf"), width = 10, height = 10)
  gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
  circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
  
  circos.genomicInitialize(custom_ideogram)
  
  # Plot Athila (if present)
  if(athila_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "Athila", custom_ideogram), 
                          count_by = "number", col = "#FF1493", 
                          track.height = 0.07, window.size = 1e5) # DeepPink
  }
  
  # Plot CRM (if present)
  if(crm_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "CRM", custom_ideogram), 
                          count_by = "number", col = "#4B0082", 
                          track.height = 0.07, window.size = 1e5) # Indigo
  }
  
  # Create legend only for present clades
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
  cat("Available clades:", paste(unique(gff_enriched$Clade[!is.na(gff_enriched$Clade)]), collapse = ", "), "\n")
}

# ============================================================================
# SECTION: Statistical analysis of TE clades
# ============================================================================

# Create a summary table with counts per superfamily and clade
clade_summary <- gff_enriched %>%
  filter(!is.na(Clade)) %>%
  group_by(Superfamily = V3, Clade) %>%
  summarise(
    Count = n(),
    Complete_elements = sum(Complete == "yes", na.rm = TRUE),
    Incomplete_elements = sum(Complete == "no", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Superfamily, desc(Count))

# Separate statistics for Copia and Gypsy
copia_summary <- clade_summary %>%
  filter(str_detect(Superfamily, "Copia")) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 2))

gypsy_summary <- clade_summary %>%
  filter(str_detect(Superfamily, "Gypsy")) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 2))

# Overall clade summary (all superfamilies)
overall_clade_summary <- gff_enriched %>%
  filter(!is.na(Clade)) %>%
  group_by(Clade) %>%
  summarise(
    Total_Count = n(),
    Complete_elements = sum(Complete == "yes", na.rm = TRUE),
    Incomplete_elements = sum(Complete == "no", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Count)) %>%
  mutate(Percentage = round(Total_Count / sum(Total_Count) * 100, 2))

# Print summaries to console
cat("\n" , rep("=", 80), "\n")
cat("COPIA ELEMENTS BY CLADE\n")
cat(rep("=", 80), "\n\n")
print(copia_summary, n = Inf)
cat("\nTotal Copia elements with clade classification:", sum(copia_summary$Count), "\n")

cat("\n" , rep("=", 80), "\n")
cat("GYPSY ELEMENTS BY CLADE\n")
cat(rep("=", 80), "\n\n")
print(gypsy_summary, n = Inf)
cat("\nTotal Gypsy elements with clade classification:", sum(gypsy_summary$Count), "\n")

cat("\n" , rep("=", 80), "\n")
cat("OVERALL CLADE SUMMARY (ALL SUPERFAMILIES)\n")
cat(rep("=", 80), "\n\n")
print(overall_clade_summary, n = Inf)

cat("\n" , rep("=", 80), "\n")
cat("MOST ABUNDANT CLADES IN THE GENOME:\n")
cat(rep("=", 80), "\n")
top_5_clades <- head(overall_clade_summary, 5)
for(i in 1:nrow(top_5_clades)) {
  cat(sprintf("%d. %s: %d elements (%.2f%%)\n", 
              i, 
              top_5_clades$Clade[i], 
              top_5_clades$Total_Count[i],
              top_5_clades$Percentage[i]))
}

# ============================================================================
# SAVE RESULTS TO FILES
# ============================================================================

# Create results directory if it doesn't exist
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Save Copia summary
write.table(copia_summary, 
            file = file.path(results_dir, "Copia_clades_summary.tsv"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Save Gypsy summary
write.table(gypsy_summary, 
            file = file.path(results_dir, "Gypsy_clades_summary.tsv"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Save overall summary
write.table(overall_clade_summary, 
            file = file.path(results_dir, "Overall_clades_summary.tsv"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Create a comprehensive report
comprehensive_report <- list(
  "Copia_Clades" = copia_summary,
  "Gypsy_Clades" = gypsy_summary,
  "All_Clades" = overall_clade_summary
)

# Save as Excel file (optional, requires openxlsx package)
if(require(openxlsx, quietly = TRUE)) {
  write.xlsx(comprehensive_report, 
             file = file.path(results_dir, "TE_clades_comprehensive_report.xlsx"))
  cat("\nExcel report saved to:", file.path(results_dir, "TE_clades_comprehensive_report.xlsx"), "\n")
} else {
  cat("\nNote: Install 'openxlsx' package to generate Excel report\n")
}

cat("\nTSV files saved to:\n")
cat("  -", file.path(results_dir, "Copia_clades_summary.tsv"), "\n")
cat("  -", file.path(results_dir, "Gypsy_clades_summary.tsv"), "\n")
cat("  -", file.path(results_dir, "Overall_clades_summary.tsv"), "\n")

# ============================================================================
# OPTIONAL: Create a visualization of clade abundance
# ============================================================================

# Barplot of top clades
pdf(file.path(plot_dir, "04-Clade_abundance_barplot.pdf"), width = 12, height = 8)

# Top 15 clades
top_15_clades <- head(overall_clade_summary, 15)

par(mar = c(8, 5, 4, 2))
barplot(top_15_clades$Total_Count,
        names.arg = top_15_clades$Clade,
        las = 2,
        col = rainbow(15),
        main = "Top 15 Most Abundant TE Clades",
        ylab = "Number of Elements",
        cex.names = 0.8)

# Add percentage labels on top of bars
text(x = 1:nrow(top_15_clades) * 1.2 - 0.5,
     y = top_15_clades$Total_Count,
     labels = paste0(top_15_clades$Percentage, "%"),
     pos = 3,
     cex = 0.7)

dev.off()

# Comparison plot: Copia vs Gypsy clades
pdf(file.path(plot_dir, "05-Copia_vs_Gypsy_clades.pdf"), width = 14, height = 6)

par(mfrow = c(1, 2))

# Copia clades
if(nrow(copia_summary) > 0) {
  par(mar = c(8, 5, 4, 2))
  barplot(copia_summary$Count,
          names.arg = copia_summary$Clade,
          las = 2,
          col = "#DC143C",
          main = "Copia Clades Distribution",
          ylab = "Number of Elements",
          cex.names = 0.8)
  text(x = 1:nrow(copia_summary) * 1.2 - 0.5,
       y = copia_summary$Count,
       labels = paste0(copia_summary$Percentage, "%"),
       pos = 3,
       cex = 0.7)
}

# Gypsy clades
if(nrow(gypsy_summary) > 0) {
  par(mar = c(8, 5, 4, 2))
  barplot(gypsy_summary$Count,
          names.arg = gypsy_summary$Clade,
          las = 2,
          col = "#2E8B57",
          main = "Gypsy Clades Distribution",
          ylab = "Number of Elements",
          cex.names = 0.8)
  text(x = 1:nrow(gypsy_summary) * 1.2 - 0.5,
       y = gypsy_summary$Count,
       labels = paste0(gypsy_summary$Percentage, "%"),
       pos = 3,
       cex = 0.7)
}

dev.off()

cat("\nBarplots saved to:\n")
cat("  -", file.path(plot_dir, "04-Clade_abundance_barplot.pdf"), "\n")
cat("  -", file.path(plot_dir, "05-Copia_vs_Gypsy_clades.pdf"), "\n")

# ============================================================================
# FINAL SUMMARY STATISTICS
# ============================================================================

cat("\n" , rep("=", 80), "\n")
cat("FINAL SUMMARY\n")
cat(rep("=", 80), "\n\n")

total_TEs_in_gff <- nrow(gff_data)
total_TEs_with_clade <- sum(!is.na(gff_enriched$Clade))
total_copia <- sum(copia_summary$Count)
total_gypsy <- sum(gypsy_summary$Count)

cat("Total TE annotations in GFF3:", total_TEs_in_gff, "\n")
cat("Total TEs with clade classification:", total_TEs_with_clade, 
    sprintf("(%.2f%%)\n", total_TEs_with_clade/total_TEs_in_gff*100))
cat("Total Copia elements classified:", total_copia, "\n")
cat("Total Gypsy elements classified:", total_gypsy, "\n")
cat("\nNumber of different clades identified:", nrow(overall_clade_summary), "\n")
cat("  - Copia clades:", nrow(copia_summary), "\n")
cat("  - Gypsy clades:", nrow(gypsy_summary), "\n")

cat("\n" , rep("=", 80), "\n\n")