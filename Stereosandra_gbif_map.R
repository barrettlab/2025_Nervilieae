# ---- 1. Load Required Packages ----
install.packages(c("rgbif", "dplyr", "sf", "ggplot2", 
                   "rnaturalearth", "rnaturalearthdata", "ggrepel"))

library(rgbif)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# ---- 2. Download Occurrence Data ----
spp_name <- "Stereosandra javanica"

occ <- occ_search(scientificName = spp_name, hasCoordinate = TRUE, limit = 1000)
occ_data <- occ$data

# ---- 3. Clean and Convert to sf ----
occ_clean <- occ_data %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>%
  select(species, decimalLatitude, decimalLongitude)

occ_sf <- st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# ---- 4. Get Basemap and Set Extent ----
world <- ne_countries(scale = "medium", returnclass = "sf")

bbox <- st_bbox(occ_sf)
bbox_expanded <- bbox + c(-2, -2, 2, 2)  # Small buffer around occurrences

# ---- 5. Get and Filter Country Labels ----
country_labels <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  filter(name %in% c("Indonesia", "Malaysia", "Thailand", "Vietnam", "Philippines", "Papua New Guinea"))

centroids_coords <- st_coordinates(country_labels)

country_labels_df <- country_labels %>%
  st_drop_geometry() %>%
  bind_cols(as.data.frame(centroids_coords))

# ---- 6. Plot Map ----
ggplot(data = world) +
  geom_sf(fill = "gray95", color = "gray50") +
  geom_sf(data = occ_sf, color = "darkred", size = 2, alpha = 0.8) +
  geom_text_repel(
    data = country_labels_df,
    aes(x = X, y = Y, label = name),
    size = 3, color = "black", fontface = "italic", max.overlaps = 10
  ) +
  coord_sf(xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
           ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"])) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of *Stereosandra javanica*",
       caption = "Occurrence data from GBIF") +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue"),
    plot.title = element_text(face = "bold", size = 16)
  )
