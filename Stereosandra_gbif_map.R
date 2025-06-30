# Package citations:
# rgbif: Chamberlain, S. A., et al. (2024). rgbif: Interface to the Global Biodiversity Information Facility API. R package version 3.7.7. https://CRAN.R-project.org/package=rgbif
# dplyr: Wickham, H., et al. (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.4. https://CRAN.R-project.org/package=dplyr
# sf: Pebesma, E. (2018). Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal*, 10(1), 439–446. https://doi.org/10.32614/RJ-2018-009
# ggplot2: Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org
# rnaturalearth: South, A. (2017). rnaturalearth: World Map Data from Natural Earth. R package version 0.1.0. https://github.com/ropensci/rnaturalearth
# rnaturalearthdata: Natural Earth (Public domain). https://www.naturalearthdata.com/downloads/
# ggrepel: Slowikowski, K. (2024). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.5. https://CRAN.R-project.org/package=ggrepel
# OpenAI. (2024). ChatGPT (April 2024 version) [Large language model]. https://chat.openai.com. Used for code debugging.

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
