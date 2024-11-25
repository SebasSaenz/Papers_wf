# Johan S. Sáenz, University of Hohenheim
# Calculate the genome size of isolates, MAGs, Bacteria and Archaea

# Load libraries ---------------------------------------------------------------
pacman::p_load(
  tidyverse,
  data.table,
  ggtext,
  patchwork,
  viridis
)

colors <- c("#1f78b4", "#1b9e77")
colors_isolate <- c("grey", "black")

#  Load other data frames ------------------------------------------------------
quality_genomes <- read_tsv("defence_rummen/rawdata/quality_report_genomes.tsv") %>%
  rename(genome = Name)

taxonomy <- read_tsv("defence_rummen/rawdata/new_gtdbtk.bac120.tsv") %>%
  rename(genome = user_genome)

hq_genomes <- read_tsv("defence_rummen/rawdata/hq_genomes.txt")

# Clean taxonomy data frame
taxonomy <- taxonomy %>%
  select(genome, classification) %>%
  rename_all(tolower) %>%
  mutate(classification = str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into = c(
             "kingdom", "phylum", "class", "order", "family", "genus",
             "species"
           ),
           sep = ";"
  ) %>%
  inner_join(hq_genomes)

# Genome size ------------------------------------------------------------------
size_kingdom <- quality_genomes %>% 
  right_join(taxonomy, by = "genome") %>%
  select(kingdom, Genome_Size) %>%
  mutate(kbp = Genome_Size/1000000) %>% 
  ggplot(aes(
    x = kingdom,
    y = kbp,
    fill = kingdom
  )) +
  geom_violin(
    aes(fill = kingdom, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  geom_line(data=tibble(x=c(1,2), 
                        y=c(7.5, 7.5)), 
            aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=7.7), 
            aes(x=x, y=y, label="***"), inherit.aes=FALSE) +
  scale_y_continuous(limits = c(0 ,8),
                     breaks = seq(0, 8, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(
    y = "Genome size (Mbp)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )


kruskal.test(kbp ~ kingdom, data = size_kingdom)


size_isolates <- quality_genomes %>% 
  right_join(taxonomy, by = "genome") %>%
  select(genome, kingdom, Genome_Size) %>%
  mutate(kbp = Genome_Size/1000000) %>% 
  mutate(isolate = case_when(grepl("MGY", genome) ~ as.character("MAG"),
                             .default = as.character("Isolate"))) %>% 
  ggplot(aes(
    x = isolate,
    y = kbp,
    fill = isolate
  )) +
  geom_violin(
    aes(fill = isolate, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  geom_line(data=tibble(x=c(1,2), 
                                    y=c(7.5, 7.5)), 
                        aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=7.7), 
            aes(x=x, y=y, label="***"), inherit.aes=FALSE) +
  scale_y_continuous(limits = c(0 ,8),
                     breaks = seq(0, 8, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colors_isolate) +
  labs(
    y = "Genome size (Mbp)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(kbp ~ isolate, data = size_isolates)


# Compose figure

size_plot <- size_kingdom + size_isolates +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))


ggsave(size_plot, file = "defence_rummen/new_plots/genome_size.png",
       width = 6, height = 4, dpi = 350)
