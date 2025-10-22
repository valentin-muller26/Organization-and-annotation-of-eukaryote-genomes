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

# Plot the distribution of Athila and CRM clades (known centromeric TEs in Brassicaceae).
# You need to run the TEsorter on TElib to get the clades classification from the TE library