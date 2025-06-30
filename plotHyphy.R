
# R and Unix code to conduct HyPhy/RELAX analysis for NEottieae and Cephalanthera austiniae


### Software and Tools

- **R v4.3.2 (R Core Team 2023)** — used for tree manipulation, alignment processing, and string handling:
  - `ape` (Paradis & Schliep 2019)
  - `Biostrings` (Pages et al. 2023)
  - `stringr` (Wickham 2019)
  - `fs` (Bryan & Wickham 2023)

- **HyPhy v2.5.50 (Kosakovsky Pond et al. 2005)** — installed via Bioconda; used to run the RELAX model (Wertheim et al. 2015) for detecting shifts in selection intensity.

- **HyPhy Vision** (Spielman & Kosakovsky Pond 2018) — browser-based interactive tool for visualizing JSON output from RELAX and other HyPhy models ([https://vision.hyphy.org](https://vision.hyphy.org)).

- **GNU Bash v5.2** — used for automating RELAX runs across multiple loci in Unix environments.


# Remove stops
for aln in "$MACSE_OUT"/*_NT.fasta; do
    base=$(basename "$aln" _NT.fasta)
    
    java -Xmx32g -jar "$MACSE_JAR" \
        -prog exportAlignment \
		-codonForFinalStop --- \
		-codonForInternalStop --- \
        -align "$aln" \
        -out_NT "$MACSE_OUT/${base}_macse_nostop_NT.fasta"
done

# Concatenate
AMAS.py concat \
  -i *nostop_NT.fasta \
  -f fasta \
  -d dna \
  -u fasta \
  -t supermatrix.fasta \
  -p partitions.txt
  

# HyPhy RELAX

  
library(ape)
library(Biostrings)
library(stringr)

# Step 1: Read the full tree
full_tree <- read.tree("supermatrix_rooted.tre")

# Step 2: List all FASTA files
fasta_files <- list.files(pattern = "\\.fasta$|\\.fa$", full.names = TRUE)

# Step 2.5: Make sure hyphy_trees folder exists
if (!dir.exists("hyphy_trees")) {
  dir.create("hyphy_trees")
}

# Step 3: Loop over each alignment
for (fasta_file in fasta_files) {
  
  # Extract base name of the file
  locus_name <- tools::file_path_sans_ext(basename(fasta_file))
  
  # Read alignment
  alignment <- readDNAStringSet(fasta_file)
  alignment_taxa <- names(alignment)
  
  # Clean taxon names (keep everything before first space)
  alignment_taxa <- sapply(strsplit(alignment_taxa, " "), `[`, 1)

  # Prune tree
  taxa_to_drop <- setdiff(full_tree$tip.label, alignment_taxa)
  pruned_tree <- drop.tip(full_tree, taxa_to_drop)
  
  # Save pruned tree
  output_file <- file.path("hyphy_trees", paste0(locus_name, ".tre"))
  write.tree(pruned_tree, file = output_file)
  
  cat("Saved pruned tree for", locus_name, "to", output_file, "\n")
}



# Load libraries
library(ape)
library(Biostrings)
library(stringr)

# List all FASTA files again
fasta_files <- list.files(pattern = "\\.fasta$|\\.fa$", full.names = TRUE)

# Initialize a log to track any issues
verification_log <- list()

# Loop over each locus
for (fasta_file in fasta_files) {
  
  # Extract base name
  locus_name <- tools::file_path_sans_ext(basename(fasta_file))
  
  # Read alignment
  alignment <- readDNAStringSet(fasta_file)
  alignment_taxa <- names(alignment)
  alignment_taxa <- sapply(strsplit(alignment_taxa, " "), `[`, 1)
  
  # Read pruned tree
  tree_file <- file.path("hyphy_trees", paste0(locus_name, ".tre"))
  if (!file.exists(tree_file)) {
    cat("Warning: Tree file missing for", locus_name, "\n")
    next
  }
  pruned_tree <- read.tree(tree_file)
  tree_taxa <- pruned_tree$tip.label
  
  # Compare sets
  taxa_in_alignment_not_tree <- setdiff(alignment_taxa, tree_taxa)
  taxa_in_tree_not_alignment <- setdiff(tree_taxa, alignment_taxa)
  
  if (length(taxa_in_alignment_not_tree) == 0 && length(taxa_in_tree_not_alignment) == 0) {
    cat(locus_name, "- Taxa match perfectly.\n")
  } else {
    cat(locus_name, "- Taxa mismatch detected.\n")
    verification_log[[locus_name]] <- list(
      alignment_not_in_tree = taxa_in_alignment_not_tree,
      tree_not_in_alignment = taxa_in_tree_not_alignment
    )
  }
}

# Optionally: view mismatches
if (length(verification_log) > 0) {
  cat("\nSummary of mismatches:\n")
  print(verification_log)
} else {
  cat("\nAll files matched perfectly!\n")
}


# Create output folder if it doesn't exist
output_dir <- "hyphy_labelled_trees"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define your myco taxa
myco_taxa <- "Stereosandra_javanica_OS2_Ste_osandra"

# List all rooted tree files
tree_files <- list.files("hyphy_trees", pattern = "\\.tre$", full.names = TRUE)

# Loop through each tree
for (tree_file in tree_files) {
  
  # Read tree as text
  tree_text <- readLines(tree_file)
  
  # Combine into a single string (in case split over lines)
  tree_text <- paste(tree_text, collapse = "")
  
  # Modify tree text: add {myco} to the specified taxa
  for (taxon in myco_taxa) {
    tree_text <- str_replace_all(
      tree_text,
      pattern = paste0("(", taxon, ")(:)", collapse = ""),
      replacement = "\\1{myco}\\2"
    )
  }
  
  # Remove any trailing "Root" label
  tree_text <- str_replace(tree_text, "\\)Root;", ");")
  
  # Save cleaned tree
  output_path <- file.path(output_dir, basename(tree_file))
  writeLines(tree_text, output_path)
  
  cat("Processed tree saved:", output_path, "\n")
}



#!/bin/bash

# Activate conda env (if not already activated outside script)
# conda activate /usr/local/src/conda_envs/binf/

# Create output folder
mkdir -p relax_output

# Loop through all fasta alignments
for aln in nostops/*.fasta; do
    # Extract base filename (no path or extension)
    base=$(basename "$aln" .fasta)

    # Define corresponding tree file
    tree="hyphy_labelled_trees/${base}.tre"

    # Define output JSON
    out="relax_output/${base}_RELAX.json"

    # Check if both alignment and tree exist
    if [[ -f "$tree" ]]; then
        echo "Running RELAX for $base"
        hyphy relax \
            --alignment "$aln" \
            --tree "$tree" \
            --output "$out" \
            --branches Test  # Will use {myco} annotation for test branches \
    else
        echo "Tree not found for $base — skipping"
    fi
done

# Summarize the .json output

conda install -c conda-forge jq  # Conda

#!/bin/bash

# Directory with RELAX output
relax_dir="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/CDS_NUC/macse/relax_output"

# Output table
output_file="relax_summary.tsv"
echo -e "Gene\tK-value\tStatus\tLRT\tp-value" > "$output_file"

# Loop through all *_RELAX.json files
for json_file in "$relax_dir"/*_RELAX.json; do
    gene=$(basename "$json_file" _RELAX.json)

    # Use jq to extract values
    k=$(jq -r '.["test results"]["relaxation or intensification parameter"]' "$json_file")
    lrt=$(jq -r '.["test results"]["LRT"]' "$json_file")
    pval=$(jq -r '.["test results"]["p-value"]' "$json_file")

    # Determine intensification or relaxation
    if [[ $(echo "$k > 1" | bc -l) -eq 1 ]]; then
        status="Intensification"
    else
        status="Relaxation"
    fi

    echo -e "${gene}\t${k}\t${status}\t${lrt}\t${pval}" >> "$output_file"
done

echo "✅ Summary written to $output_file"

# Plot in R

library(tidyverse)

# Load data
df <- read.table("relax_summary.tsv", header = TRUE, sep = "\t")

# Add significance, functional category, and direction for coloring
df <- df %>%
  mutate(
    significant = p < 0.05,
    category = case_when(
      str_detect(Gene, "psa|psb|pet|rbc|atp") ~ "Photosynthesis",
      TRUE ~ "Housekeeping"
    ),
    direction = ifelse(K > 1, "Intensification", "Relaxation")
  )

# Custom centered log scale for K
df <- df %>%
  mutate(K_log = case_when(
    K < 1 ~ -log10(1 / K),
    K == 1 ~ 0,
    K > 1 ~ log10(K)
  ))

# Reorder gene labels
ordered_genes <- df %>%
  distinct(Gene, category) %>%
  arrange(factor(category, levels = c("Photosynthesis", "Housekeeping")), Gene) %>%
  pull(Gene)

df$Gene <- factor(df$Gene, levels = rev(ordered_genes))

# Plot
ggplot(df, aes(x = K_log, y = Gene)) +
  geom_segment(aes(x = 0, xend = K_log, yend = Gene), color = "gray70") +
  geom_point(aes(color = direction, shape = significant), size = 3, stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) +
  scale_color_manual(values = c("Relaxation" = "red", "Intensification" = "blue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
    labels = c("0.03", "0.1", "0.3", "1", "3", "10", "30"),
    name = "K value (log-scaled, centered at K = 1)"
  ) +
  labs(
    y = "Gene (locus)",
    shape = "Significant (p < 0.05)",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")



