library(tidyverse)

# Define genes of interest
g2_genes <- c("matk", "clpp", "trna-ugc", "trni-gau", "trnk-uuu", "trnv-uac", "atpf", "rpl2")

# Filter to those genes only, recode presence/absence
df_g2 <- plot_df_collapsed %>%
  filter(gene %in% g2_genes) %>%
  mutate(status = ifelse(value == 1, "present", "absent")) %>%
  select(taxon, gene, status)

# Pivot wider
matrix_df <- df_g2 %>%
  pivot_wider(names_from = gene, values_from = status, values_fill = "absent")

# Filter taxa for either condition
filtered_df <- matrix_df %>%
  filter(
    (matk == "present" & rowSums(select(., -taxon, -matk) == "absent") > 0) |
    (matk == "absent" & rowSums(select(., -taxon, -matk) == "present") > 0)
  )

# Order taxa: matk absent first
filtered_df <- filtered_df %>%
  arrange(matk)

# Pivot longer for plotting
plot_long <- filtered_df %>%
  pivot_longer(cols = -taxon, names_to = "gene", values_to = "status")

# Rank non-matk genes by presence frequency (in filtered data only)
gene_order <- plot_long %>%
  filter(gene != "matk", status == "present") %>%
  count(gene, sort = TRUE) %>%
  pull(gene)

# Create final gene order with matk first
gene_levels <- c("matk", gene_order)

# Ensure gene is a factor with correct levels
plot_long <- plot_long %>%
  filter(gene %in% gene_levels) %>%
  mutate(
    gene = factor(gene, levels = gene_levels),
    taxon = factor(taxon, levels = unique(filtered_df$taxon))
  )

# Assign fill colors
plot_long <- plot_long %>%
  mutate(fill_group = case_when(
    gene == "matk" & status == "present" ~ "matk_present",
    gene == "matk" & status == "absent" ~ "absent",
    gene != "matk" & status == "present" ~ "other_present",
    TRUE ~ "absent"
  ))

# Plot
ggplot(plot_long, aes(x = gene, y = taxon, fill = fill_group)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      matk_present = "#2166ac",     # blue
      other_present = "#1b7837",    # green
      absent = "gray90"
    ),
    name = "Status",
    labels = c("Absent", "matK present", "Other gene present")
  ) +
  labs(
    x = "Gene", y = "Taxon",
    title = "Presence/Absence Heatmap of matK and Group IIA Intron Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
