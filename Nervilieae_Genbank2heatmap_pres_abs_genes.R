
# Code for plotting presence/absence of putatively functional genes across the orchids, focusing on mycoheterotrophic/leafless taxa
# Craig F. Barrett 2025-04-30

# ---- SOFTWARE & RESOURCES USED ----

# tidyverse (v2.0.0)
# Wickham H, Averick M, Bryan J, Chang W, McGowan LDA, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019).
# "Welcome to the tidyverse." *Journal of Open Source Software*, 4(43), 1686. https://doi.org/10.21105/joss.01686

# stringr (v1.5.1)
# Wickham H (2019). *stringr: Simple, Consistent Wrappers for Common String Operations*. R package version 1.5.1. https://CRAN.R-project.org/package=stringr

# dplyr (v1.1.4)
# Wickham H, François R, Henry L, Müller K (2023). *dplyr: A Grammar of Data Manipulation*. https://CRAN.R-project.org/package=dplyr

# tidyr (v1.3.1)
# Wickham H, Girlich M (2023). *tidyr: Tidy Messy Data*. https://CRAN.R-project.org/package=tidyr

# ggplot2 (v3.4.4)
# Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. ISBN 978-3-319-24277-4

# tibble (v3.2.1)
# Müller K, Wickham H (2023). *tibble: Simple Data Frames*. https://CRAN.R-project.org/package=tibble

# readr (v2.1.4)  [Optional if you read CSV/TSV inputs]
# Wickham H, Hester J, Bryan J (2023). *readr: Read Rectangular Text Data*. https://CRAN.R-project.org/package=readr

# External Software and Resources:
# - ChatGPT-4o (OpenAI, 2024). Used to assist in writing and debugging R code and designing visualization workflows.
#   Citation: OpenAI. 2024. *ChatGPT-4o*. https://openai.com/chatgpt

# - NCBI GenBank. Source of plastid genome annotations.
#   Citation: Benson DA, Cavanaugh M, Clark K, et al. GenBank. *Nucleic Acids Res.* 2018;46(D1):D41–D47.
#   https://doi.org/10.1093/nar/gkx1094

# - Geneious v10 (Biomatters Ltd). Used for initial genome annotation, inspection, and taxon thinning.
#   Citation: Kearse M, Moir R, Wilson A, et al. 2012. Geneious Basic: An integrated and extendable desktop software platform for the organization and analysis of sequence data.
#   *Bioinformatics* 28(12):1647–1649. https://doi.org/10.1093/bioinformatics/bts199


# Load libraries and check versions

# Load and print package versions
packages <- c("tidyverse", "stringr", "dplyr", "tidyr", "ggplot2", "tibble")

invisible(lapply(packages, library, character.only = TRUE))
sapply(packages, packageVersion)



# ---- STEP 1: Parse multi-genbank file (// is the delimiter between genomes) ----

# Search NCBI GenBank [Orchidaceae NOT UNVERIFIED ], then download GenBank Full format
# The genbank file is: Orchidaceae_plastomes_thinned.gb. Taxa were thinned to 1 leafy representative per genus in Geneious v. 10 but keeping all leafless species of interest

# Read entire GenBank file as lines
gb_lines <- readLines("Nervilieae_all.gb")

# Find record delimiters (assumes // ends each record)
record_ends <- grep("^//", gb_lines)
record_starts <- c(1, record_ends + 1)
record_ends <- c(record_ends, length(gb_lines))

# Create output directory
dir.create("neerv_orchid_split_gb", showWarnings = FALSE)

# Save each genome to its own file
for (i in seq_along(record_starts)) {
  record <- gb_lines[record_starts[i]:record_ends[i]]
  # Try to extract taxon name from LOCUS or SOURCE or ORGANISM
  locus_line <- grep("^LOCUS", record, value = TRUE)
  taxon_line <- grep("^ +ORGANISM", record, value = TRUE)
  taxon <- if (length(taxon_line) > 0) {
    gsub("^\\s*ORGANISM\\s*", "", taxon_line[1])
  } else if (length(locus_line) > 0) {
    gsub("^LOCUS\\s+|\\s.*$", "", locus_line[1])
  } else {
    paste0("Unknown_", i)
  }
  taxon <- gsub("[^A-Za-z0-9_]", "_", taxon)  # Sanitize
  outfile <- file.path("neerv_orchid_split_gb", paste0(taxon, ".gb"))
  writeLines(record, outfile)
}

### STOP -- need to delete the "UNKNOWN" gb file at the end of list! Then proceed.

# ---- STEP 2: Parse genes and genome lengths ----

gb_files <- list.files("neerv_orchid_split_gb", pattern = "\\.gb$", full.names = TRUE)
target_types <- c("CDS", "tRNA", "rRNA")
gene_records <- list()
genome_lengths <- c()

for (file in gb_files) {
  lines <- readLines(file)
  taxon <- gsub("\\.gb$", "", basename(file))

  if (taxon %in% c("Microtis_unifolia", "Corybas_diemenicus", "Corybas_cheesemanii", "Chiloglottis_cornuta", "Danhatchia_australis")) next

  locus_line <- grep("^LOCUS", lines, value = TRUE)
  genome_length <- as.numeric(str_extract(locus_line, "\\d+"))
  genome_lengths[taxon] <- genome_length

  features_start <- grep("^FEATURES", lines)
  origin_start <- grep("^ORIGIN", lines)
  if (length(features_start) == 0 || length(origin_start) == 0) next
  features <- lines[(features_start + 1):(origin_start - 1)]

  gene_list <- list()
  current_type <- NA
  current_gene <- NA

  for (line in features) {
    if (grepl("^\\s{5}[a-zA-Z_]+", line)) {
      current_type <- trimws(substr(line, 6, 20))
      current_gene <- NA
    }
    if (grepl("/gene=", line)) {
      current_gene <- gsub('.*="/?|"$', '', line)
    }
    if (!is.na(current_gene) && current_type %in% target_types) {
      gene_list[[current_gene]] <- 1
      current_gene <- NA
    }
  }

  if (length(gene_list) > 0) {
    gene_df <- tibble(
      taxon = taxon,
      gene = names(gene_list),
      value = 1
    )
    gene_records[[taxon]] <- gene_df
  }
}

# ---- STEP 3: Create binary matrix ----

all_genes <- bind_rows(gene_records)

mat_bin <- all_genes %>%
  complete(taxon, gene, fill = list(value = 0)) %>%
  pivot_wider(names_from = gene, values_from = value, values_fill = 0)

# Order genes by frequency
gene_totals <- colSums(as.matrix(mat_bin[,-1]))
gene_order <- names(sort(gene_totals, decreasing = TRUE))

# Order genomes by total length
taxa_order <- names(sort(genome_lengths[mat_bin$taxon], decreasing = TRUE))

# Long format for plotting
long_df <- mat_bin %>%
  pivot_longer(-taxon, names_to = "gene", values_to = "value") %>%
  mutate(taxon = factor(taxon, levels = taxa_order),
         gene = factor(gene, levels = gene_order))

# ---- STEP 4: Plot using ggplot2 ----

# ggplot(long_df, aes(x = gene, y = taxon, fill = factor(value))) +
#   geom_tile(color = "black") +
#   scale_fill_manual(values = c("0" = "white", "1" = "black"),
#                     name = "Gene",
#                     labels = c("Absent", "Present")) +
#   theme_minimal(base_size = 10) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     panel.grid = element_blank(),
#     legend.position = "bottom"
#   ) +
#   labs(
#     x = "Gene",
#     y = "Plastid Genome (Orchidaceae)",
#     title = "Plastid Gene Presence/Absence Matrix (Orchidaceae)"
#   )


##################
##################

# Real, cleaned gene names
real_genes <- c(
  "accd","atpa","atpb","atpe","atpf","atph","atpi","ccsa","cema","clpp","infa","matk",
  "ndha","ndhb","ndhc","ndhd","ndhe","ndhf","ndhg","ndhh","ndhi","ndhj","ndhk",
  "peta","petb","petd","petg","petl","petn",
  "psaa","psab","psac","psai","psaj",
  "psba","psbb","psbc","psbd","psbe","psbf","psbh","psbi","psbj","psbk","psbl","psbm","psbn","psbt","psbz",
  "rbcl","rpl14","rpl16","rpl2","rpl20","rpl22","rpl23","rpl32","rpl33","rpl36",
  "rpoa","rpob","rpoc1","rpoc2",
  "rps11","rps12","rps14","rps15","rps16","rps18","rps19","rps2","rps3","rps4","rps7","rps8",
  "rrn16","rrn23","rrn4.5","rrn5",
  "trna-ugc","trnc-gca","trnd-guc","trne-uuc","trnfm-cau","trnf-gaa","trng-gcc","trng-ucc","trnh-gug",
  "trni-cau","trni-gau","trnk-uuu","trnl-caa","trnl-uaa","trnl-uag","trnm-cau","trnn-guu","trnp-ugg",
  "trnq-uug","trnr-acg","trnr-ucu","trns-gcu","trns-gga","trns-uga","trnt-ggu","trnt-ugu","trnv-gac",
  "trnv-uac","trnw-cca","trny-gua",
  "ycf1","ycf2","ycf3","ycf4"
)

# Named vector to map orphan names to real names
orphan_map <- c(
  "clpp1" = "clpp",
  "pbf1" = "psbn",
  "pafii" = "ycf4",
  "pafi" = "ycf3",
  "trnl-gag" = "trnl-uag",
  "trns-cga" = "trns-gga",
  "trna" = "trna-ugc",
  "trni" = "trni-gau",
  "trnd" = "trnd-guc",
  "trne" = "trne-uuc",
  "trnf" = "trnf-gaa",
  "trnf-aaa" = "trnf-gaa",
  "trnk" = "trnk-uuu",
  "trnl-cag" = "trnl-caa",
  "trnn" = "trnn-guu",
  "trnp-ggg" = "trnp-ugg",
  "trnq" = "trnq-uug",
  "trnt-cgu" = "trnt-ggu",
  "trnw" = "trnw-cca",
  "trny" = "trny-gua",
  "lhba" = "psbz",
  "trni-cag" = "trni-cau",
  "trnfm" = "trnfm-cau",
  "rrn16s" = "rrn16",
  "rrn23s" = "rrn23",
  "rrn4.5s" = "rrn4.5",
  "rrn5s" = "rrn5",
  "16srrna" = "rrn16",
  "23srrna" = "rrn23",
  "4.5srrna" = "rrn4.5",
  "5srrna" = "rrn5"
)

# Genes to remove entirely
remove_genes <- c("ycf15", "orf42", "orf56", "ycf68", "nad6", "pebt", "rrn26s", "trnc", "trng", "trnr", "ycf9")

# After cleaning/mapping genes
long_df <- mat_bin %>%
  pivot_longer(-taxon, names_to = "gene", values_to = "value") %>%
  mutate(
    gene = str_to_lower(gene),
    gene = str_replace_all(gene, "_", "-"),
    gene = str_replace_all(gene, "\\(", "-"),
    gene = str_replace_all(gene, "\\)", ""),
    gene = str_replace_all(gene, "\\s", ""),
    gene = recode(gene, !!!orphan_map)
  ) %>%
  filter(!gene %in% remove_genes)

# Collapse to one row per taxon-gene combo
plot_df <- long_df %>%
  group_by(taxon, gene) %>%
  summarise(value = max(value), .groups = "drop")


# Apply full cleanup
plot_df <- plot_df %>%
  mutate(
    gene = str_to_lower(gene),
    gene = str_replace_all(gene, "_", "-"),
    gene = str_replace_all(gene, "\\(", "-"),
    gene = str_replace_all(gene, "\\)", ""),
    gene = str_replace_all(gene, "\\s", ""),
    gene = recode(gene, !!!orphan_map)  # collapse orphans to real names
  ) %>%
  filter(!gene %in% remove_genes)  # remove undesired genes

# Collapse to get binary presence
plot_df_collapsed <- plot_df %>%
  group_by(taxon, gene) %>%
  summarise(value = max(value), .groups = "drop")

# Define highlighted genes
highlight_genes <- c("accd", "clpp", "trne-uuc", "ycf1", "ycf2")

# Add a 'color group' column to encode colors
plot_df_collapsed <- plot_df_collapsed %>%
  mutate(
    color_group = case_when(
      value == 0 ~ "absent",
      gene %in% highlight_genes ~ "highlighted",
      TRUE ~ "present"
    )
  )

# Set custom colors
color_map <- c(
  absent = "white",
  present = "gray45",
  highlighted = "dodgerblue"
)

# Redraw plot
plot_df_collapsed %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    taxon = factor(taxon, levels = taxa_order)
  ) %>%
  ggplot(aes(x = gene, y = taxon, fill = color_group)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = color_map,
                    name = "Gene Presence",
                    labels = c("Absent", "Present", "Core non-bioenergetic")) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Gene (Standardized)",
    y = "Taxon",
    title = "Plastid Gene Presence/Absence (Orchidaceae)",
    subtitle = "Core non-bioenergetic genes highlighted in blue"
  )



# Define highlighted core genes
highlight_genes <- c("accd", "clpp", "trne-uuc", "ycf1", "ycf2")

# Add a 'color group' column to define fill categories
plot_df_colored <- plot_df_collapsed %>%
  mutate(
    color_group = case_when(
      value == 0 ~ "absent",
      gene %in% highlight_genes ~ "highlighted",
      TRUE ~ "present"
    )
  )

# Order genes and taxa
gene_order <- plot_df_colored %>%
  group_by(gene) %>%
  summarise(freq = sum(color_group != "absent")) %>%
  arrange(desc(freq)) %>%
  pull(gene)

taxa_order <- plot_df_colored %>%
  group_by(taxon) %>%
  summarise(n = sum(color_group != "absent")) %>%
  arrange(desc(n)) %>%
  pull(taxon)

# Plot
plot_df_colored %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    taxon = factor(taxon, levels = taxa_order)
  ) %>%
  ggplot(aes(x = gene, y = taxon, fill = color_group)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c(
    "absent" = "white",
    "present" = "gray45",
    "highlighted" = "dodgerblue"
  ),
  name = "Gene Presence",
  labels = c("Absent", "Present", "Core non-bioenergetic")) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Gene (Standardized)",
    y = "Taxon",
    title = "Standardized Gene Presence/Absence Matrix (Orchidaceae)",
    subtitle = "Core non-bioenergetic genes highlighted in blue"
  )
