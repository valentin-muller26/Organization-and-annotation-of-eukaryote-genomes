# Load required packages
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Create plot directory if it doesn't exist
plot_dir <- "plots/circlize_with_genes"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ============================================================================
# LOAD DATA FILES
# ============================================================================

# Load the TE annotation GFF3 file
gff_file <- "HiFiasm_Lu1_primary.fa.mod.EDTA.TEanno.gff3"
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Load MAKER gene annotation GFF3 file
gene_gff_file <- "filtered.genes.renamed.gff3"
gene_data <- read.table(gene_gff_file, header = FALSE, sep = "\t", 
                        stringsAsFactors = FALSE, comment.char = "#")

# Filter only gene features from MAKER annotation
gene_features <- gene_data %>%
  filter(V3 == "gene")

cat("\nGene annotations loaded:", nrow(gene_features), "genes\n")
cat("TE annotations loaded:", nrow(gff_data), "TEs\n")

# Custom ideogram data
custom_ideogram <- read.table("HiFiasm_Lu1_primary.fa.fai", header = FALSE, 
                              stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = T), ]

# Select only the first 15 longest scaffolds
custom_ideogram <- custom_ideogram[1:15, ]

# ============================================================================
# CREATE COMBINED TIR CATEGORY
# ============================================================================

# Identify all TIR transposons and combine them
tir_categories <- c("Mutator_TIR_transposon", "hAT_TIR_transposon", 
                    "PIF_Harbinger_TIR_transposon", "CACTA_TIR_transposon", 
                    "Tc1_Mariner_TIR_transposon")

# Create a new column for combined categories
gff_data_combined <- gff_data %>%
  mutate(V3_combined = ifelse(V3 %in% tir_categories, "TIR_transposon_combined", V3))

cat("\nOriginal TIR counts:\n")
for(tir in tir_categories) {
  count <- sum(gff_data$V3 == tir)
  cat(tir, ":", count, "\n")
}

cat("\nCombined TIR count:", sum(gff_data_combined$V3_combined == "TIR_transposon_combined"), "\n")

# ============================================================================
# FILTERING FUNCTIONS
# ============================================================================

# Function to filter GFF3 data based on Superfamily
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
  filtered_data <- gff_data[gff_data$V3 == superfamily, ] %>%
    as.data.frame() %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# Function to filter GFF3 data based on Superfamily (using combined categories)
filter_superfamily_combined <- function(gff_data, superfamily, custom_ideogram) {
  filtered_data <- gff_data[gff_data$V3_combined == superfamily, ] %>%
    as.data.frame() %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# Function to filter genes
filter_genes <- function(gene_data, custom_ideogram) {
  filtered_data <- gene_data %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V7) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

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
# LOAD CLADE DATA FOR ATHILA/CRM PLOT
# ============================================================================

# Load TEsorter classification files
copia_clades <- read_tsv("Copia.cls.tsv", comment = "", show_col_types = FALSE)
gypsy_clades <- read_tsv("Gypsy.cls.tsv", comment = "", show_col_types = FALSE)

cat("\nCopia clades loaded:", nrow(copia_clades), "rows\n")
cat("Gypsy clades loaded:", nrow(gypsy_clades), "rows\n")

# Combine both classifications
all_clades <- bind_rows(copia_clades, gypsy_clades)

# Extract TE_ID from the first column
all_clades$TE_ID <- str_extract(all_clades[[1]], "TE_[0-9]+")

# Extract TE_ID from GFF3 annotation
gff_data$TE_ID <- str_extract(gff_data$V9, "(?<=Name=)TE_[0-9]+")

# Merge GFF3 with clade information
gff_enriched <- gff_data %>%
  left_join(all_clades %>% select(TE_ID, Clade, Complete), by = "TE_ID")

# Check if Athila and CRM are present
athila_count <- sum(gff_enriched$Clade == "Athila", na.rm = TRUE)
crm_count <- sum(gff_enriched$Clade == "CRM", na.rm = TRUE)

cat("\nAthila elements found:", athila_count, "\n")
cat("CRM elements found:", crm_count, "\n")

# ============================================================================
# PLOT 1: ALL TE SUPERFAMILIES (WITH COMBINED TIR) + GENES
# ============================================================================

cat("\n=== Creating Plot 1: All Superfamilies + Genes ===\n")

# Get the most abundant superfamilies (using combined categories)
superfamily_counts <- gff_data_combined$V3_combined %>% table() %>% sort(decreasing = TRUE)

cat("\nSuperfamily counts (with combined TIR):\n")
print(superfamily_counts[1:10])

# Select top superfamilies to plot
priority_superfamilies <- c("Gypsy_LTR_retrotransposon", 
                            "Copia_LTR_retrotransposon",
                            "helitron",
                            "TIR_transposon_combined",
                            "rRNA_gene",
                            "tRNA_SINE_retrotransposon")

# Keep only those that exist in the data
superfamilies_to_plot <- priority_superfamilies[priority_superfamilies %in% names(superfamily_counts)]

# Categories to exclude
exclude_categories <- c("repeat_fragment", "LTR_retrotransposon", "Long_terminal_repeat", 
                        "LTR", "long_terminal_repeat")

cat("\nSuperfamilies to plot:", length(superfamilies_to_plot), "\n")
print(superfamilies_to_plot)

# Define colors for each superfamily + genes
superfamily_colors <- c(
  "Gypsy_LTR_retrotransposon" = "#2E8B57",      # SeaGreen
  "Copia_LTR_retrotransposon" = "#DC143C",      # Crimson
  "helitron" = "brown",                        # DodgerBlue
  "TIR_transposon_combined" = "#FF8C00",        # DarkOrange - COMBINED TIR
  "rRNA_gene" = "violet",                      # LightSeaGreen
  "tRNA_SINE_retrotransposon" = "black",
  "LINE_retrotransposon" = "yellow"
)

# Assign colors to superfamilies
plot_colors <- c()
for(sf in superfamilies_to_plot) {
  if(sf %in% names(superfamily_colors)) {
    plot_colors <- c(plot_colors, superfamily_colors[sf])
  } else {
    plot_colors <- c(plot_colors, rainbow(20)[length(plot_colors) + 1])
  }
}

# Add gene color
plot_colors <- c(plot_colors, "#4169E1")  # RoyalBlue for genes

pdf(file.path(plot_dir, "TE_density_all_superfamilies_combined_TIR_with_genes.pdf"), 
    width = 18, height = 18)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

# Plot each superfamily
plotted_superfamilies <- c()
plotted_colors <- c()

for(i in seq_along(superfamilies_to_plot)) {
  filtered_data <- filter_superfamily_combined(gff_data_combined, superfamilies_to_plot[i], custom_ideogram)
  
  # Only plot if data exists
  if(nrow(filtered_data) > 0) {
    circos.genomicDensity(filtered_data, 
                          count_by = "number", col = plot_colors[i], 
                          track.height = 0.06, window.size = 1e5)
    
    plotted_superfamilies <- c(plotted_superfamilies, superfamilies_to_plot[i])
    plotted_colors <- c(plotted_colors, plot_colors[i])
    
    cat("Plotted", superfamilies_to_plot[i], ":", nrow(filtered_data), "elements\n")
  } else {
    cat("Warning:", superfamilies_to_plot[i], "has no data, skipping...\n")
  }
}

# Add genes track
circos.genomicDensity(filter_genes(gene_features, custom_ideogram), 
                      count_by = "number", col = "#4169E1", 
                      track.height = 0.06, window.size = 1e5)

plotted_superfamilies <- c(plotted_superfamilies, "Genes")
plotted_colors <- c(plotted_colors, "#4169E1")

# Clean up superfamily names for legend
legend_names <- plotted_superfamilies
legend_names <- gsub("_LTR_retrotransposon", " LTR-RT", legend_names)
legend_names <- gsub("_retrotransposon", "", legend_names)
legend_names <- gsub("_gene", "", legend_names)
legend_names <- gsub("TIR_transposon_combined", "TIR transposons (all)", legend_names)
legend_names <- tools::toTitleCase(legend_names)

lgd_all <- Legend(
  title = "Genomic Features", 
  title_gp = gpar(fontsize = 22, fontface = "bold"),
  at = legend_names,
  legend_gp = gpar(fill = plotted_colors),
  labels_gp = gpar(fontsize = 18),
  ncol = 1,
  grid_height = unit(12, "mm"),
  grid_width = unit(12, "mm")
)
draw(lgd_all, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("\nPlot 1 saved to:", file.path(plot_dir, "TE_density_all_superfamilies_combined_TIR_with_genes.pdf"), "\n")

# ============================================================================
# PLOT 2: ATHILA, CRM CLADES AND GENES
# ============================================================================

cat("\n=== Creating Plot 2: Athila, CRM and Genes ===\n")

if(athila_count > 0 || crm_count > 0) {
  pdf(file.path(plot_dir, "Athila_CRM_and_Genes.pdf"), width = 18, height = 18)
  gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
  circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
  
  circos.genomicInitialize(custom_ideogram)
  
  # Plot Athila (if present)
  if(athila_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "Athila", custom_ideogram), 
                          count_by = "number", col = "#FF1493", 
                          track.height = 0.08, window.size = 1e5)
    cat("Plotted Athila:", athila_count, "elements\n")
  }
  
  # Plot CRM (if present)
  if(crm_count > 0) {
    circos.genomicDensity(filter_clade(gff_enriched, "CRM", custom_ideogram), 
                          count_by = "number", col = "#4B0082", 
                          track.height = 0.08, window.size = 1e5)
    cat("Plotted CRM:", crm_count, "elements\n")
  }
  
  # Plot Genes
  circos.genomicDensity(filter_genes(gene_features, custom_ideogram), 
                        count_by = "number", col = "#4169E1", 
                        track.height = 0.08, window.size = 1e5)
  cat("Plotted Genes:", nrow(gene_features), "genes\n")
  
  # Create legend
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
  present_clades <- c(present_clades, "Genes")
  present_colors <- c(present_colors, "#4169E1")
  
  lgd_centro <- Legend(
    title = "Genomic Features", 
    title_gp = gpar(fontsize = 22, fontface = "bold"),
    at = present_clades,
    legend_gp = gpar(fill = present_colors),
    labels_gp = gpar(fontsize = 18),
    ncol = 1,
    grid_height = unit(12, "mm"),
    grid_width = unit(12, "mm")
  )
  draw(lgd_centro, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))
  
  circos.clear()
  dev.off()
  
  cat("Plot 2 saved to:", file.path(plot_dir, "Athila_CRM_and_Genes.pdf"), "\n")
} else {
  cat("\nWarning: No Athila or CRM elements found. Skipping Plot 2.\n")
}

