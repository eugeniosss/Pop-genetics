# Check if 'ggVennDiagram' is installed
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  # Install remotes package if not already installed
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }
  
  # Install ggVennDiagram from GitHub
  remotes::install_github("gaospecial/ggVennDiagram")
}

# Load required libraries
library(tidyverse)
library(gridExtra)
library(grid)
library(optparse)
library(VennDiagram)
library(ggplot2)
library(ggVennDiagram)
library(UpSetR)

option_list <- list(
  make_option(c("-f", "--folders"), type = "character", action = "store", default = "", 
              help = "One or more folder paths (separated by multiple -f options)"),
  make_option(c("-m", "--matrix"), type = "character", help = "Path to intersects matrix file"),
  make_option(c("-o", "--output"), type = "character", help = "Output file name")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure required arguments are provided
#if (is.null(opt$folders) || is.null(opt$matrix) || is.null(opt$output)) {
#  print_help(opt_parser)
#  stop("All arguments (-f, -m, -o) must be provided.", call. = FALSE)
#}

#input_dirs=opt$folders
#output_pdf=opt$output

## Print parsed inputs for debugging
input_dirs <- unlist(strsplit(opt$folders, ","))

# Define function to read and annotate files
read_file <- function(directory, filename, dataset) {
    file_path <- file.path(directory, filename)
    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }
    data <- read_delim(file_path, delim = "\t", skip = 1, col_names = FALSE) %>%
        mutate(Dataset = dataset)
    return(data)
}

# Initialize empty list
all_data <- list(
    var_qual = tibble(),
    var_depth = tibble(),
    var_miss = tibble(),
    var_freq = tibble(),
    ind_depth = tibble(),
    ind_miss = tibble()
)

dataset_names <- gsub("/stats/", "", input_dirs)

# Read and merge data from all dataset directories
for (i in seq_along(input_dirs)) {

    dir <- input_dirs[i]
    dataset <- dataset_names[i]

    all_data$var_qual <- bind_rows(all_data$var_qual, read_file(dir, "site_quality.lqual", dataset))
    all_data$var_depth <- bind_rows(all_data$var_depth, read_file(dir, "mean_depth_per_site.ldepth.mean", dataset))
    all_data$var_miss <- bind_rows(all_data$var_miss, read_file(dir, "missing_data_per_site.lmiss", dataset))
    all_data$var_freq <- bind_rows(all_data$var_freq, read_file(dir, "allele_freq2.frq", dataset))
    all_data$ind_depth <- bind_rows(all_data$ind_depth, read_file(dir, "mean_depth_per_individual.idepth", dataset))
    all_data$ind_miss <- bind_rows(all_data$ind_miss, read_file(dir, "missing_data_per_indv.imiss", dataset))
}

#all_data

# Create a PDF device
pdf(opt$output)

# Generate plots with NA filtering
plot_list <- list(
    ggplot(all_data$var_qual %>% filter(!is.na(X3)), aes(x = X3, color = Dataset)) +
        stat_density(aes(x=X3, colour=Dataset), geom="line",position="identity") +
        labs(title = "Variant Quality Distribution", x = "Quality Score", y = "Density") +
        theme_light(),

    ggplot(all_data$var_depth %>% filter(!is.na(X3)), aes(x = X3, color = Dataset)) +
        stat_density(aes(x=X3, colour=Dataset), geom="line",position="identity") +
        labs(title = "Variant Depth Distribution", x = "Mean Depth", y = "Density") +
        theme_light(),

    ggplot(all_data$var_miss %>% filter(!is.na(X6)), aes(x = X6, color = Dataset)) +
        stat_density(aes(x=X6, colour=Dataset), geom="line",position="identity") +
        labs(title = "Variant Missingness Distribution", x = "Proportion of Missing Data", y = "Density") +
        theme_light(),

    ggplot(all_data$var_freq %>% mutate(MAF = pmin(X5, X6)) %>% filter(!is.na(MAF)), aes(x = MAF, color = Dataset)) +
        stat_density(aes(x=MAF, colour=Dataset), geom="line",position="identity") +
        labs(title = "Minor Allele Frequency Distribution", x = "Minor Allele Frequency", y = "Density") +
        theme_light(),

    ggplot(all_data$ind_depth %>% filter(!is.na(X3)), aes(x = X3, fill = Dataset)) +
        geom_histogram(position = "identity", alpha = 0.25, bins = 30, na.rm = TRUE) +
        labs(title = "Mean Depth per Individual", x = "Mean Depth", y = "Frequency") +
        theme_light(),

    ggplot(all_data$ind_miss %>% filter(!is.na(X5)), aes(x = X5, fill = Dataset)) +
        geom_histogram(position = "identity", alpha = 0.25, bins = 30, na.rm = TRUE) +
        labs(title = "Missing Data per Individual", x = "Proportion of Missing Data", y = "Frequency") +
        theme_light()
)

# Print each plot to PDF
for (p in plot_list) {
    print(p)
}

grid.newpage()

grid.text("Missing Data per Individual", 
          x = 0.2, y = 0.95)

missing_data_table <- all_data$ind_miss %>%
    arrange(desc(X5)) %>%
    mutate(ind = as.character(X1)) %>%
    select(ind, X5, Dataset) %>%
    pivot_wider(
        names_from = Dataset,      # Make the software names the column names
        values_from = X5,        # Values are from the 'fmiss' (Missing_Data)
        values_fill = list(fmiss = -1)  # Fill NA values with 0 (or another value if you prefer)
    ) 

missing_data_table <- missing_data_table %>%
    select(ind, sort(names(.)[-1])) 

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


grid.newpage()

grid.text("Mean depth per Individual", 
          x = 0.2, y = 0.95)

missing_data_table <- all_data$ind_depth %>%
    arrange(desc(X3)) %>%
    mutate(ind = as.character(X1)) %>%
    select(ind, X3, Dataset) %>%
    pivot_wider(
        names_from = Dataset,      # Make the software names the column names
        values_from = X3,        # Values are from the 'fmiss' (Missing_Data)
        values_fill = list(fmiss = -1)  # Fill NA values with 0 (or another value if you prefer)
    ) 

missing_data_table <- missing_data_table %>%
    select(ind, sort(names(.)[-1])) 

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

# Read SNP presence file (assuming tab-separated)
snp_data <- read.table(opt$matrix, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)

# Convert columns into lists of SNPs per software
snp_sets <- lapply(snp_data, function(col) col[!is.na(col)])

# Name the lists by the column names (software names)
names(snp_sets) <- colnames(snp_data)

snp_sets_cleaned <- lapply(snp_sets, function(x) x[x != ""])

# Now calculate total SNPs per software (length of cleaned lists)
total_snps_cleaned <- sapply(snp_sets_cleaned, length)

# Create a data frame for plotting
snps_df_cleaned <- data.frame(
  Dataset = names(total_snps_cleaned),
  Total_SNPs = total_snps_cleaned
)

ggplot(snps_df_cleaned, aes(x = Dataset, y = Total_SNPs, fill = Dataset)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Total Variants per Dataset", x = "Dataset", y = "Total Variants") +
  theme(axis.text.x = element_blank())

ggVennDiagram(snp_sets, label_alpha=0) + 
  theme_minimal() +
  ggtitle("Variants Overlap Across Datasets") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_distiller(palette = "RdBu")

# Close the PDF **after** all plots are added
dev.off()