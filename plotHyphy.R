
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


#----------------------------------------------------------------------
# 1. Drop tips from the full tree if the taxa are not in the alignment, for each file
#----------------------------------------------------------------------


library(ape)
library(Biostrings)
library(stringr)

# Step 1: Read the full tree
full_tree <- read.tree("mixfinder.contree")

# Step 2: List all FASTA files
fasta_files <- list.files("codon_aligned", pattern = "\\.fasta$|\\.fa$", full.names = TRUE)

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


#----------------------------------------------------------------------
# 2. Verify trees and alignments contain identical taxa
#----------------------------------------------------------------------


# Load libraries
library(ape)
library(Biostrings)
library(stringr)

# List all FASTA files again
fasta_files <- list.files("codon_aligned", pattern = "\\.fasta$|\\.fa$", full.names = TRUE)

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
    cat("✅", locus_name, "- Taxa match perfectly.\n")
  } else {
    cat("❌", locus_name, "- Taxa mismatch detected.\n")
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

#----------------------------------------------------------------------
# 3. Set test branches with [myco] to indicate for HyPhy/RELAX
#----------------------------------------------------------------------


library(ape)
library(stringr)

# Define your myco taxa
myco_taxa <- "Stereosandra_javanica"
# Create output folder if it doesn't exist
output_dir <- "hyphy_labelled_trees"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List all rooted tree files
tree_files <- list.files("hyphy_trees_rooted", pattern = "\\.tre$", full.names = TRUE)

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
  
  cat("✅ Processed tree saved:", output_path, "\n")
}

#----------------------------------------------------------------------
# 4. Run HyPhy/RELAX in a loop with 32 threads
#----------------------------------------------------------------------


#!/bin/bash

# Activate conda env (if not already activated outside script)
# conda activate /usr/local/src/conda_envs/binf/

# Create output folder
mkdir -p relax_output

# Loop through all fasta alignments
for aln in codon_aligned/*.fasta; do
    # Extract base filename (no path or extension)
    base=$(basename "$aln" .fasta)

    # Define corresponding tree file
    tree="hyphy_labelled_trees/${base}.tre"

    # Define output JSON
    out="relax_output/${base}_RELAX.json"

    # Check if both alignment and tree exist
    if [[ -f "$tree" ]]; then
        echo "🧪 Running RELAX for $base"
        hyphy relax \
            --alignment "$aln" \
            --tree "$tree" \
            --output "$out" \
            --branches "Test"  # Will use {myco} annotation for test branches
    else
        echo "⚠️  Tree not found for $base — skipping"
    fi
done

#####


#----------------------------------------------------------------------
# 5. Prep trees with ONLY Cephalanthera austiniae as the test branch
#----------------------------------------------------------------------

# remove {myco} flag from all but C austiniae

library(stringr)
library(fs)

# Define folders
input_folder <- "hyphy_labelled_trees"
output_folder <- "hyphy_ceph_labelled"

# Create output folder if it doesn't exist
dir_create(output_folder)

# List all tree files (assumes .tre extension)
tree_files <- dir_ls(input_folder, glob = "*.tre")

# Define the taxon that should keep {myco}
target_taxon <- "Cephalanthera_austiniae_2010_Ce_stiniae"

for (file in tree_files) {
  # Read tree as single string
  tree_str <- readLines(file, warn = FALSE)
  tree_str <- paste(tree_str, collapse = "")  # In case it's multi-line

  # Remove all occurrences of {myco}
  tree_str <- str_replace_all(tree_str, "\\{myco\\}", "")

  # Re-insert {myco} for the target taxon
  tree_str <- str_replace_all(tree_str,
    paste0("(", target_taxon, ")([^\\w\\{])"),
    paste0("\\1{myco}\\2")
  )

  # Write to new file in output folder
  out_file <- path(output_folder, path_file(file))
  writeLines(tree_str, con = out_file)
}

#----------------------------------------------------------------------
# 6. Run hyphy with C austiniae as the test branch, and only run for datasets that include C austiniae (in Unix)
#----------------------------------------------------------------------

#!/bin/bash

mkdir -p relax_ceph_output

for aln in codon_aligned/*.fasta; do
    base=$(basename "$aln" .fasta)

    tree="hyphy_ceph_labelled/${base}.tre"
    out="relax_ceph_output/${base}_RELAX.json"

    if [[ -f "$tree" ]]; then
        echo "🧪 Running RELAX for $base"

        OMP_NUM_THREADS=32 hyphy relax \
            --alignment "$aln" \
            --tree "$tree" \
            --output "$out" \
            --branches myco \
            --test myco

    else
        echo "Tree not found for $base — skipping"
    fi
done

#----------------------------------------------------------------------
# 7. Get list of all fasta files that contain C austiniae (in R)
#----------------------------------------------------------------------

library(Biostrings)

# Define path
fasta_dir <- "codon_aligned"
files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Target taxon
target <- "Cephalanthera_austiniae_2010_Ce_stiniae"

# Create named output vector
annotated_files <- sapply(files, function(f) {
  aln <- readDNAStringSet(f)
  fname <- basename(f)
  if (any(names(aln) == target)) {
    paste0(fname, " *")  # Append asterisk
  } else {
    fname
  }
})

# Print results
cat(annotated_files, sep = "\n")

# Upload all 


#----------------------------------------------------------------------
# 8. Prep data to run hyphy/RELAX with only C austiniae and no other MH taxa (in R)
#----------------------------------------------------------------------

library(Biostrings)
library(ape)
library(stringr)
library(fs)

# Taxa to remove
remove_taxa <- c(
  "Aphyllorchis_montana_KU551262",
  "Cephalanthera_humilis_KU551265",
  "Diplandrorchis_sinica_MZ014629",
  "Diplandrorchis_sinica_OP310037",
  "Neottia_listeroides_KU551272",
  "Neottia_camtschatea_KU551266",
  "Neottia_nidus_avis_JF325876",
  "Neottia_acuminata_KU551268",
  "Limodorum_abortivum_MH590355"
)

# Input/output folders
aln_in <- "codon_aligned"
aln_out <- "codon_aligned_noMH"
tree_in <- "hyphy_ceph_labelled"
tree_out <- "hyphy_ceph_labelled_noMH"

# Create output dirs
dir_create(aln_out)
dir_create(tree_out)

# ---------------------------
# 1. Process FASTA alignments
aln_files <- dir_ls(aln_in, glob = "*.fasta")

for (file in aln_files) {
  aln <- readDNAStringSet(file)
  aln_filtered <- aln[!names(aln) %in% remove_taxa]
  
  # Save filtered alignment
  writeXStringSet(aln_filtered, file = file.path(aln_out, path_file(file)))
}

# ---------------------------
# 2. Process TREE files
tree_files <- dir_ls(tree_in, glob = "*.tre")

for (file in tree_files) {
  tree_str <- readLines(file, warn = FALSE)
  tree_str <- paste(tree_str, collapse = "")
  # Remove each taxon using drop.tip
  tree <- read.tree(text = tree_str)
  keep <- setdiff(tree$tip.label, remove_taxa)
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, keep))
  
  # Preserve {myco} labels if present
  newick_out <- write.tree(tree_pruned)
  
  writeLines(newick_out, file.path(tree_out, path_file(file)))
}

#----------------------------------------------------------------------
# 9. Run hyphy/RELAX
#----------------------------------------------------------------------

#!/bin/bash

mkdir -p relax_ceph_output_noMH

for aln in codon_aligned_noMH/*.fasta; do
    base=$(basename "$aln" .fasta)

    tree="hyphy_ceph_labelled_noMH/${base}.tre"
    out="relax_ceph_output_noMH/${base}_RELAX.json"

    if [[ -f "$tree" ]]; then
        echo "🧪 Running RELAX for $base"

        OMP_NUM_THREADS=32 hyphy relax \
            --alignment "$aln" \
            --tree "$tree" \
            --output "$out" \
            --branches myco \
            --test myco

    else
        echo "⚠️  Tree not found for $base — skipping"
    fi
done


# Re-run ycf1

# Including all MH taxa as test branches
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/codon_aligned/ycf1_mafft_macse_NT.fasta
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/hyphy_labelled_trees/ycf1_mafft_macse_NT.tre

# Cephalanthera austiniae only, includes other MH taxa
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/codon_aligned/ycf1_mafft_macse_NT.fasta
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/hyphy_ceph_labelled/ycf1_mafft_macse_NT.tre

# Cephalanthera austiniae only, no other MH taxa
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/codon_aligned_noMH/ycf1_mafft_macse_NT.fasta
/Data/cbarrett/001_2024_orchid_plastomes/2025_Cephalanthera/hyphy_ceph_labelled_noMH/ycf1_mafft_macse_NT.tre

# Download/install macse_v2
wget https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
sudo cp /usr/local/bin

# Including all MH taxa as test branches
java -jar /usr/local/bin/macse_v2.07.jar \
  -prog alignSequences \
  -seq codon_aligned/ycf1_mafft_macse_NT.fasta \
  -out_NT codon_aligned/ycf1_macse_aligned_NT.fasta \
  -out_AA codon_aligned/ycf1_macse_aligned_AA.fasta

# Cephalanthera austiniae only, no other MH taxa
# Align No-MH dataset
java -jar /usr/local/bin/macse_v2.07.jar \
  -prog alignSequences \
  -seq codon_aligned_noMH/ycf1_mafft_macse_NT.fasta \
  -out_NT codon_aligned_noMH/ycf1_macse_aligned_NT.fasta \
  -out_AA codon_aligned_noMH/ycf1_macse_aligned_AA.fasta



#----------------------------------------------------------------------
# 9. Upload .json files to http://vision.hyphy.org/RELAX
#    Report intensification/relaxation, K, p, and LR
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# 10. #Plotting HyPhy results in R
#----------------------------------------------------------------------

# Table = HypPhy.csv

library(tidyverse)

# Load and reshape
df <- read.csv("HypPhy.csv")
colnames(df) <- c(
  "locus", "taxa", "sites",
  "type_AllMH", "K_AllMH", "p_AllMH", "LR_AllMH",
  "type_Caustiniae", "K_Caustiniae", "p_Caustiniae", "LR_Caustiniae",
  "type_CaustiniaeNoMH", "K_CaustiniaeNoMH", "p_CaustiniaeNoMH", "LR_CaustiniaeNoMH"
)

df_long <- df %>%
  pivot_longer(
    cols = -c(locus, taxa, sites),
    names_to = c(".value", "group"),
    names_sep = "_"
  ) %>%
  mutate(
    significant = p < 0.05,
    # Create functional category
    category = case_when(
      str_starts(locus, "psa|psb|pet|rbc|atp") ~ "Photosynthesis",
      TRUE ~ "Housekeeping"
    )
  )

# Custom transformation of K (centered log scale)
df_long <- df_long %>%
  mutate(K_log = case_when(
    K < 1 ~ -log10(1 / K),
    K == 1 ~ 0,
    K > 1 ~ log10(K)
  ))

# Reorder locus: first by category, then alphabetically
ordered_loci <- df_long %>%
  distinct(locus, category) %>%
  arrange(factor(category, levels = c("Photosynthesis", "Housekeeping")), locus) %>%
  pull(locus)

df_long$locus <- factor(df_long$locus, levels = rev(ordered_loci))  # reverse for top-down plotting

# Set colors
group_colors <- c("AllMH" = "blue", "Caustiniae" = "orange", "CaustiniaeNoMH" = "red")

# Plot
ggplot(df_long, aes(x = K_log, y = locus, color = group)) +
  geom_segment(aes(x = 0, xend = K_log, yend = locus), color = "gray70") +
  geom_point(aes(shape = significant), size = 3, stroke = 0.8, fill = "black") +
  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1, 1.5),
    labels = c("0.1", "0.3", "1", "3", "10", "30"),
    name = "K value (log-scaled, centered at K = 1)"
  ) +
  labs(y = "Gene (locus)", shape = "Significant (p < 0.05)", color = "Group") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
