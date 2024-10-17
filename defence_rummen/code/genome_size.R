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