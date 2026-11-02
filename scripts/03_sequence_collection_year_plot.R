library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ---- Load data ----
df <- read_csv("C:/Users/Dell/Downloads/Asia WNV project/WNV_Envelope_Mutation_Summary with Accession.csv",
               show_col_types = FALSE)

# ---- Prepare data ----
plot_df <- df %>%
  mutate(
    Year = as.numeric(str_sub(Collection_Date, 1, 4)),   # handles YYYY, YYYY-MM, YYYY-MM-DD
    Category = case_when(
      str_detect(Range_Classification, "^Historical") ~ "Historical",
      str_detect(Range_Classification, "^Novel")      ~ "Novel",
      TRUE ~ NA_character_                              # drops "Unknown" rows
    )
  ) %>%
  filter(!is.na(Year), !is.na(Category)) %>%
  mutate(Category = factor(Category, levels = c("Historical", "Novel"))) %>%
  arrange(Year, Category) %>%                           # historical dots at the bottom of each stack
  group_by(Year) %>%
  mutate(Count = row_number()) %>%                      # one dot per sequence, stacked
  ungroup()

# ---- Plot ----
p <- ggplot(plot_df, aes(x = Year, y = Count, color = Category)) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("Historical" = "#F28E2B", "Novel" = "#3A78B5")) +
  scale_x_continuous(breaks = seq(1955, 2025, 5),
                     limits = c(1953, 2027), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title    = "WNV Sequences by Collection Year",
    subtitle = "Each point represents one sequence",
    x = "Collection year",
    y = "Sequences collected",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 10),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.key = element_blank()
  )

print(p)

# ---- Save ----
ggsave("WNV_sequences_by_year.png", p, width = 7.5, height = 4.5, dpi = 300)