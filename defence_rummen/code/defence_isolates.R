# Johan S. Sáenz, University of Hohenheim
# Calculate the number of defence systems between isolates and MAGs


# Load libraries ---------------------------------------------------------------
pacman::p_load(
  tidyverse,
  data.table,
  pheatmap,
  ggridges,
  ggtext,
  ggdist,
  gghalves,
  patchwork,
  viridis
)

# Concatenate files-----------------------------------------------------
all_paths <- list.files(
  path = "defence_rummen/rawdata/padloc/csv/",
  pattern = "*.csv", # List paths
  full.names = TRUE
)

all_content <- all_paths %>% # Join files
  lapply(read.table,
         header = TRUE,
         sep = ",",
         encoding = "UTF-8"
  )

all_filenames <- all_paths %>% # Manipulate file paths
  basename() %>%
  as.list()

all_lists <- mapply(c,
                    all_content,
                    all_filenames,
                    SIMPLIFY = FALSE
)

all_result <- rbindlist(all_lists, fill = T)
names(all_result)[20] <- "genome" # change column name

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

# Size of genomes

df_isolates <- hq_genomes %>% 
  mutate(isolate = case_when(grepl("MGY", genome) ~ as.character("MAG"),
                             .default = as.character("Isolate")))

number_isolates <- hq_genomes %>% 
  mutate(isolate = case_when(grepl("MGY", genome) ~ as.character("MAG"),
                             .default = as.character("Isolate"))) %>% 
  select(isolate) %>%
  count(isolate, name = "total")


# Read a file for DS families-------------------------------------------------
list_system <- read_tsv("defence_rummen/rawdata/list_system.txt")

# Select colors
colors_isolate <- c("grey", "black")
colors <- c("#1f78b4", "#1b9e77")
# Prevalence plot


# Number of genomes wiht DS
barplots_total_ds_iso <- all_result %>%
  filter(!grepl("_other", system)) %>%
  select(genome) %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  unique() %>%
  inner_join(df_isolates, by = "genome") %>%
  select(isolate) %>%
  count(isolate) %>%
  inner_join(number_isolates, by = "isolate") %>%
  mutate(percent = 100 * (n / total)) %>%
  ggplot(aes(
    x = isolate,
    y = percent,
    fill = isolate
  )) +
  geom_col() +
  scale_fill_manual(values = colors_isolate) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    y = "Genomes with defence system (%)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )


# Number of systems per genome
number_system_iso <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(df_isolates, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  select(isolate, n) %>%
  ggplot(aes(
    x = isolate,
    y = n,
    fill = isolate
  )) +
  geom_violin(
    aes(fill = isolate, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  scale_fill_manual(values = colors_isolate) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 3),
    oob = scales::squish
  ) +
  geom_line(data=tibble(x=c(1,2), 
                        y=c(35, 35)), 
            aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=35.2), 
            aes(x=x, y=y, label="***"), inherit.aes=FALSE) +
  labs(
    y = "Defence system per genome",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(n ~ isolate, data = number_system)


# Number of DS families per genome
family_system_iso <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  inner_join(list_system, by = "system") %>%
  select(genome, system_family) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(df_isolates, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  ggplot(aes(
    x = isolate,
    y = n,
    fill = isolate
  )) +
  geom_violin(
    aes(fill = isolate, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  scale_fill_manual(values = colors_isolate) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 3),
    oob = scales::squish
  ) +
  geom_line(data=tibble(x=c(1,2), 
                        y=c(16, 16)), 
            aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=16.2), 
            aes(x=x, y=y, label="***"), inherit.aes=FALSE) +
  labs(
    y = "Defence system families per genome",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(n ~ isolate, data = family_system)

# kingdom plots ----------------------------------------------------------------

# Prevalence plot
kingdom_number <- taxonomy %>%
  select(kingdom) %>%
  count(kingdom, name = "total")

# Number of genomes wiht DS
barplots_total_ds <- all_result %>%
  filter(!grepl("_other", system)) %>%
  select(genome) %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  unique() %>%
  inner_join(taxonomy, by = "genome") %>%
  select(kingdom) %>%
  count(kingdom) %>%
  inner_join(kingdom_number, by = "kingdom") %>%
  mutate(percent = 100 * (n / total)) %>%
  ggplot(aes(
    x = kingdom,
    y = percent,
    fill = kingdom
  )) +
  geom_col() +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    y = "Genomes with defence system (%)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

# Number of systems per genome
number_system <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  select(kingdom, n) %>%
  ggplot(aes(
    x = kingdom,
    y = n,
    fill = kingdom
  )) +
  geom_violin(
    aes(fill = kingdom, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  
  # geom_point(
  #   shape = 21,
  #   position = position_jitter(
  #     width = 0.2,
  #     seed = 0
  #   ),
  #   size = 3,
  #   alpha = 0.8,
  #   color = "black"
  # ) +
  # stat_summary(
  #   fun = median,
  #   show.legend = FALSE,
  #   geom = "crossbar",
  #   color = "black",
  #   linewidth = 0.5
  # ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 3),
    oob = scales::squish
  ) +
  labs(
    y = "Defence system per genome",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(n ~ kingdom, data = number_system)


# Number of DS families per genome
family_system <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  inner_join(list_system, by = "system") %>%
  select(genome, system_family) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  ggplot(aes(
    x = kingdom,
    y = n,
    fill = kingdom
  )) +
  geom_violin(
    aes(fill = kingdom, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  # geom_point(
  #   shape = 21,
  #   position = position_jitter(
  #     width = 0.2,
  #     seed = 0
  #   ),
  #   size = 3,
  #   alpha = 0.8,
  #   color = "black"
  # ) +
  # stat_summary(
  #   fun = median,
  #   show.legend = FALSE,
  #   geom = "crossbar",
  #   color = "black",
  #   linewidth = 0.5
  # ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 3),
    oob = scales::squish
  ) +
  labs(
    y = "Defence system families per genome",
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(n ~ kingdom, data = family_system)

# System density ---------------------------------------------------------------

number_system_iso_density <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(df_isolates, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(quality_genomes, by = "genome") %>%  
  select(isolate, n, Genome_Size) %>%
  mutate(kbp = Genome_Size/1000,
         density = n/kbp,
         density_1000 = density*1000) %>% 
  ggplot(aes(
    x = isolate,
    y = density_1000,
    fill = isolate
  )) +
  geom_violin(
    aes(fill = isolate, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
  ) +
  geom_boxplot(
    width = .1, size = 0.3, outlier.size = 0.5
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 1) +
  scale_fill_manual(values = colors_isolate) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 13),
    breaks = seq(0, 13, 1),
    oob = scales::squish
  ) +
  geom_line(data=tibble(x=c(1,2), 
                        y=c(12, 12)), 
            aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=12.2), 
            aes(x=x, y=y, label="***"), inherit.aes=FALSE) +
  labs(
    y = bquote("Defence system / kbp"~(x10^-3)),
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(density_1000 ~ isolate, data = number_system_iso_density)

number_system_density <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(quality_genomes, by = "genome") %>%  
  select(kingdom, n, Genome_Size) %>%
  mutate(kbp = Genome_Size/1000,
         density = n/kbp,
         density_1000 = density*1000) %>% 
  ggplot(aes(
    x = kingdom,
    y = density_1000,
    fill = kingdom
  )) +
  geom_violin(
    aes(fill = kingdom, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5
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
                        y=c(12, 12)), 
            aes(x=x, y=y), inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.5, 
                        y=12.2), 
            aes(x=x, y=y, label="**"), inherit.aes=FALSE) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 13),
    breaks = seq(0, 13, 1),
    oob = scales::squish
  ) +
  labs(
    y = bquote("Defence system / kbp"~(x10^-3)),
    x = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

kruskal.test(density_1000 ~ kingdom, data = number_system_density)

# Make compose plot ------------------------------------------------------------
figure1 <- (barplots_total_ds_iso + family_system_iso + number_system_iso)/(barplots_total_ds + family_system + number_system) / 
  (number_system_density + number_system_iso_density + plot_spacer()) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave(figure1, filename = "defence_rummen/new_plots/Figure_1.png", 
       width = 8, height = 9, dpi = 350)
