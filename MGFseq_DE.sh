# Description: This script performs differential expression analysis by calculating log fold change (logFC)
#              and p-values for multiple group comparisons. It is generalized for use with any dataset.
# Author: Lissa Cruz-Saavedra
# Date: 02-04-2024
# Version: 1.0

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)

# Step 1: Load the data
# Replace "your_data.txt" with the actual file path to your dataset
file_path <- "your_data.txt"  # Example: "data/expression_data.txt"
if (!file.exists(file_path)) {
  stop("Error: File not found. Please check the file path.")
}
data <- read_delim(file_path, delim = "\t")

# Step 2: Rename the gene column (if needed)
# Replace "id" with the actual column name for genes in your dataset
data <- rename(data, Gene = id)

# Step 3: Define groups
# Replace the sample names with those relevant to your dataset
groups <- list(
  Group_A = c("Sample_A1", "Sample_A2", "Sample_A3"),
  Group_B = c("Sample_B1", "Sample_B2", "Sample_B3"),
  Group_C = c("Sample_C1", "Sample_C2", "Sample_C3"),
  Group_D = c("Sample_D1", "Sample_D2", "Sample_D3"),
  Group_E = c("Sample_E1", "Sample_E2", "Sample_E3")
)

# Step 4: Define comparisons
# Replace "Group_X" with your actual group names
comparisons <- list(
  c("Group_A", "Group_B"), c("Group_A", "Group_C"), c("Group_A", "Group_D"), c("Group_A", "Group_E"),
  c("Group_B", "Group_C"), c("Group_B", "Group_D"), c("Group_B", "Group_E"),
  c("Group_C", "Group_D"), c("Group_C", "Group_E"),
  c("Group_D", "Group_E")
)

# Step 5: Function to calculate logFC and p-value
compute_differential_expression <- function(df, group1, group2) {
  df %>%
    rowwise() %>%
    mutate(
      logFC = log2(mean(c_across(all_of(groups[[group2]]))) + 1) - log2(mean(c_across(all_of(groups[[group1]]))) + 1),
      pval = t.test(c_across(all_of(groups[[group2]])), c_across(all_of(groups[[group1]])))$p.value
    ) %>%
    ungroup()
}

# Step 6: Perform differential expression analysis for all comparisons
results <- list()
for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  cat("Performing differential expression analysis for", group1, "vs", group2, "...\n")
  
  # Perform the analysis
  result <- compute_differential_expression(data, group1, group2)
  
  # Add comparison information
  result <- result %>%
    mutate(Comparison = paste(group1, "vs", group2))
  
  # Store the result
  results[[paste(group1, group2, sep = "_vs_")]] <- result
}

# Step 7: Combine all results into a single data frame
final_results <- bind_rows(results)

# Step 8: Save the results to a file
output_file <- "differential_expression_results.txt"
write_delim(final_results, output_file, delim = "\t")
cat("Differential expression results saved to:", output_file, "\n")
