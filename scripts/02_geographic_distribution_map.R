# install.packages(c("ggplot2","dplyr","readr","sf","rnaturalearth","rnaturalearthdata","ggrepel"))
library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(rnaturalearth)
library(ggrepel)

# ---- Input ----
infile <- "WNV_Envelope_Mutation_Summary_with_Accession.csv"
df <- read_csv(infile, show_col_types = FALSE)

# ---- Counts per country ----
counts <- df %>%
  filter(!is.na(Country),
         !grepl("^Unknown", Range_Classification)) %>%
  mutate(Range = ifelse(grepl("^Historical", Range_Classification),
                        "Historical", "Novel")) %>%
  count(Country, Range, name = "n")

# ---- Bubble positions (calibrated to the reference figure) ----
pos <- tribble(
  ~Country,                           ~lon,  ~lat,
  "Belgium",                           2.3,  52.1,
  "Estonia",                          23.4,  58.0,
  "France",                           -0.7,  46.2,
  "Switzerland",                       6.0,  45.2,
  "Italy",                            11.0,  40.0,
  "Hungary",                          16.4,  46.8,
  "Slovakia",                         18.7,  49.0,
  "Greece",                           18.9,  38.6,
  "Cyprus",                           30.3,  34.4,
  "Israel",                           32.1,  30.9,
  "Iran",                             50.7,  32.0,
  "Kazakhstan",                       62.6,  47.7,
  "Russia",                           59.3,  55.4,
  "Senegal",                         -14.3,  14.5,
  "Central African Republic",         18.4,   6.5,
  "Uganda",                           30.3,   1.1,
  "Kenya",                            34.9,   0.5,
  "Democratic Republic of the Congo", 20.4,  -4.2,
  "Zambia",                           23.0, -13.4,
  "Madagascar",                       43.7, -19.4,
  "Namibia",                          15.4, -23.5,
  "South Africa",                     24.0, -29.0
)

bub <- counts %>%
  left_join(pos, by = "Country") %>%
  mutate(Range = factor(Range, levels = c("Historical", "Novel"))) %>%
  arrange(desc(n))                       # large bubbles drawn first, small on top

# ---- Base map ----
world <- ne_countries(scale = "medium", returnclass = "sf")

# ---- Plot ----
p <- ggplot() +
  # solid grey land, no country border lines
  geom_sf(data = world, fill = "#DDDDDD", colour = NA) +
  geom_point(data = bub,
             aes(lon, lat, size = n, fill = Range),
             shape = 21, colour = "white", stroke = 1) +
  geom_text(data = bub, aes(lon, lat, label = n),
            colour = "white", fontface = "bold", size = 2.6) +
  geom_text_repel(data = bub, aes(lon, lat, label = Country),
                  size = 2.8, fontface = "bold", colour = "black",
                  bg.color = "white", bg.r = 0.08,
                  box.padding = 0.9, point.padding = 1.2,
                  segment.size = 0.2, segment.color = "grey30",
                  min.segment.length = 0.5, max.overlaps = Inf, seed = 42) +
  scale_fill_manual(values = c("Historical" = "#F28E2B", "Novel" = "#2E86C1"),
                    name = "Range") +
  scale_size_continuous(range = c(4.5, 15), limits = c(1, 85),
                        breaks = c(1, 10, 30, 60, 85),
                        name = "Number of\nsequences") +
  guides(
    fill = guide_legend(order = 1,
                        override.aes = list(size = 4, colour = NA)),
    size = guide_legend(order = 2,
                        override.aes = list(fill = "grey70", colour = "grey35",
                                            stroke = 0.8))
  ) +
  coord_sf(xlim = c(-22, 92), ylim = c(-37, 70), expand = FALSE) +
  labs(title = "Geographic Distribution of WNV Sequences",
       subtitle = "Bubble size and number inside bubble = number of sequences per country") +
  theme_void(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 15,
                                    margin = margin(t = 6)),
    plot.subtitle    = element_text(hjust = 0.5, size = 8, margin = margin(b = 4)),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    legend.key.height = unit(0.9, "cm"),
    legend.position  = "right",
    plot.background  = element_rect(fill = "white", colour = "grey90", linewidth = 3),
    plot.margin      = margin(5, 5, 5, 5)
  )

print(p)
ggsave("WNV_geographic_distribution.png", p, width = 8.5, height = 6, dpi = 300)