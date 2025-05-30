# great_gene_count.R
# Analyze GREAT output: count promoter- and distal-associated genes based on distances from TSS

# Install and load required packages
required_packages <- c("tidyr", "dplyr", "stringr")
installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load GREAT output file
input_file <- "great.txt"
if (!file.exists(input_file)) {
  stop(paste("File not found:", input_file))
}
great <- read.delim(input_file, stringsAsFactors = FALSE)

# Check required column exists
column_name <- "Species.assembly..hg38"
if (!column_name %in% colnames(great)) {
  stop(paste("Column not found in file:", column_name))
}

# Determine how many columns to split based on commas
max_values <- max(stringr::str_count(great[[column_name]], ",") + 1)

# Dynamically create column names
value_cols <- paste0("value", 1:max_values)

# Separate distances and extract numeric values from gene-distance strings
great_split <- great %>%
  tidyr::separate(
    !!sym(column_name),
    into = value_cols,
    sep = ",",
    fill = "right"
  ) %>%
  dplyr::mutate(across(
    starts_with("value"),
    ~ gsub(".*\\(\\+?(-?\\d+)\\).*", "\\1", .)
  )) %>%
  dplyr::mutate(across(
    starts_with("value"),
    as.numeric
  ))

# Add promoter/distal classification columns
great_split <- great_split %>%
  dplyr::mutate(Promoter = NA_real_, Distal_region = NA_real_)

# Apply classification logic using first 3 associations
great_classified <- great_split %>%
  rowwise() %>%
  mutate(
    Promoter = coalesce(
      ifelse(abs(value1) <= 1000, value1, NA),
      ifelse(abs(value2) <= 1000, value2, NA),
      ifelse(abs(value3) <= 1000, value3, NA)
    ),
    Distal_region = coalesce(
      ifelse(abs(value1) > 1000, value1, NA),
      ifelse(abs(value2) > 1000, value2, NA),
      ifelse(abs(value3) > 1000, value3, NA)
    )
  ) %>%
  ungroup()

# Count non-NA values
promoter_count <- sum(!is.na(great_classified$Promoter))
distal_count   <- sum(!is.na(great_classified$Distal_region))

# Print results
cat("Number of genes in Promoter region (±1 kb from TSS):", promoter_count, "\n")
cat("Number of genes in Distal region (>1 kb from TSS):", distal_count, "\n")