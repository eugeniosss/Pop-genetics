# Check if input directory and output PDF file are provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Usage: script.R <input_directory> <output_pdf>")
}

# Input directory and output PDF file
input_dir <- commandArgs(trailingOnly = TRUE)[1]
output_pdf <- commandArgs(trailingOnly = TRUE)[2]

# Load required libraries
library(tidyverse)
library(gridExtra)
library(grid)       

# Define function to read files from the input directory
read_file <- function(filename) {
    file_path <- file.path(input_dir, filename)
    if (!file.exists(file_path)) {
        stop(paste("File not found:", file_path))
    }
    data <- read_delim(file_path, delim = "\t", skip = 1, col_names = FALSE)
    if (ncol(data) == 0) {
        warning(paste("File is empty:", file_path))
    }
    return(data)
}

# Create a PDF device
pdf(output_pdf)

# Variant quality
var_qual <- read_file("site_quality.lqual") %>%
    rename(chr = 1, pos = 2, qual = 3)
ggplot(var_qual, aes(qual)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Variant Quality Distribution",
         x = "Quality Score",
         y = "Density") +
    theme_light()

# Variant depth
var_depth <- read_file("mean_depth_per_site.ldepth.mean") %>%
    rename(chr = 1, pos = 2, mean_depth = 3)
ggplot(var_depth, aes(mean_depth)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Variant Depth Distribution",
         x = "Mean Depth",
         y = "Density") +
    theme_light() +
    xlim(0, 20)

# Variant missingness
var_miss <- read_file("missing_data_per_site.lmiss") %>%
    rename(chr = 1, pos = 2, nchr = 3, nfiltered = 4, nmiss = 5, fmiss = 6)
ggplot(var_miss, aes(fmiss)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Variant Missingness Distribution",
         x = "Proportion of Missing Data",
         y = "Density") +
    theme_light()

# Minor allele frequency
var_freq <- read_file("allele_freq2.frq") %>%
    rename(chr = 1, pos = 2, nalleles = 3, nchr = 4, a1 = 5, a2 = 6) %>%
    mutate(maf = pmin(a1, a2))
ggplot(var_freq, aes(maf)) +
    geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Minor Allele Frequency Distribution",
         x = "Minor Allele Frequency",
         y = "Density") +
    theme_light()

# Mean depth per individual
ind_depth <- read_file("mean_depth_per_individual.idepth") %>%
    rename(ind = 1, nsites = 2, depth = 3)
ggplot(ind_depth, aes(depth)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Mean Depth per Individual Distribution",
         x = "Mean Depth",
         y = "Frequency") +
    theme_light()

# Create a table for mean depth
grid.newpage()
mean_depth_table <- ind_depth %>%
    arrange(desc(depth)) %>%
    mutate(ind = as.character(ind)) %>%
    select(ind, depth) %>%
    rename(Mean_Depth = depth)

# Output the mean depth table to PDF
rows_per_page <- 22  # Adjust this value based on your preference

for (i in seq_len(ceiling(nrow(mean_depth_table) / rows_per_page))) {
    start_row <- (i - 1) * rows_per_page + 1
    end_row <- min(i * rows_per_page, nrow(mean_depth_table))
    grid.table(mean_depth_table[start_row:end_row, ], rows = NULL)
    if (i < ceiling(nrow(mean_depth_table) / rows_per_page)) {
        grid.newpage()
    }
}

# Missing data per individual
ind_miss <- read_file("missing_data_per_indv.imiss") %>%
    rename(ind = 1, ndata = 2, nfiltered = 3, nmiss = 4, fmiss = 5)
ggplot(ind_miss, aes(fmiss)) +
    geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
    labs(title = "Missing Data per Individual Distribution",
         x = "Proportion of Missing Data",
         y = "Frequency") +
    theme_light()

# Create a table for missing data
grid.newpage()
missing_data_table <- ind_miss %>%
    arrange(desc(fmiss)) %>%
    mutate(ind = as.character(ind)) %>%
    select(ind, fmiss) %>%
    rename(Missing_Data = fmiss)

# Output the missing data table to PDF
rows_per_page <- 22  # Adjust this value based on your preference

for (i in seq_len(ceiling(nrow(missing_data_table) / rows_per_page))) {
    start_row <- (i - 1) * rows_per_page + 1
    end_row <- min(i * rows_per_page, nrow(missing_data_table))
    grid.table(missing_data_table[start_row:end_row, ], rows = NULL)
    if (i < ceiling(nrow(missing_data_table) / rows_per_page)) {
        grid.newpage()
    }
}
# Close the PDF device
dev.off()

