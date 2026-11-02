# ---- Packages ----
pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ---- Load data ----
df <- read_csv("C:/Users/Dell/Downloads/Asia WNV project/WNV_Envelope_Mutation_Summary with Accession.csv",
               show_col_types = FALSE)

# Check: should include Range_Classification and Envelope_Substitutions
print(names(df))
print(table(df$Range_Classification))

# ---- Assign group; drop sequences with unknown range ----
df <- df %>%
  mutate(Group = case_when(
    str_detect(Range_Classification, "^Historical") ~ "Historical",
    str_detect(Range_Classification, "^Novel")      ~ "Novel",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

group_n <- df %>% count(Group, name = "Total")

# ---- One row per sequence x mutation, then count sequences carrying each ----
mut_counts <- df %>%
  filter(!is.na(Envelope_Substitutions)) %>%
  select(Accession, Group, Envelope_Substitutions) %>%
  separate_rows(Envelope_Substitutions, sep = ",\\s*") %>%
  rename(Mutation = Envelope_Substitutions) %>%
  filter(Mutation != "") %>%
  distinct(Accession, Group, Mutation) %>%
  count(Group, Mutation, name = "n")

# ---- Keep only mutations unique to one range ----
unique_muts <- mut_counts %>%
  group_by(Mutation) %>%
  filter(n_distinct(Group) == 1) %>%
  ungroup()

# ---- Top 5 per group ----
top5 <- unique_muts %>%
  left_join(group_n, by = "Group") %>%
  mutate(Pct = 100 * n / Total) %>%
  arrange(Group, desc(n), Mutation) %>%
  group_by(Group) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    Label = paste0(n, " (", round(Pct, 1), "%)"),
    Group = factor(Group, levels = c("Historical", "Novel")),
    Key   = paste(Mutation, Group, sep = "___")
  ) %>%
  arrange(Group, desc(n), Mutation) %>%
  mutate(Key = factor(Key, levels = rev(unique(Key))))

# ---- Plot ----
p <- ggplot(top5, aes(x = n, y = Key, fill = Group)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = Label), hjust = -0.08, fontface = "bold", size = 3.4) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Historical" = "#CD7F3C", "Novel" = "#2B5C8A"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
  labs(title    = "Top Recurrent Mutations by Group",
       subtitle = "Unique to each range",
       x = "Number of sequences",
       y = "Mutations") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 15),
    plot.subtitle      = element_text(size = 12, margin = margin(b = 8)),
    strip.text         = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88"),
    axis.text.y        = element_text(colour = "black", size = 10),
    axis.title.y       = element_text(face = "bold", margin = margin(r = 8)),
    axis.title.x       = element_text(size = 12, margin = margin(t = 6)),
    panel.spacing.y    = unit(1.2, "lines"),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

print(p)
ggsave("WNV_top_recurrent_mutations.png", p, width = 6, height = 5.5, dpi = 300)