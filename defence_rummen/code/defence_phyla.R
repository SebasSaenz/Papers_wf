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



# Prevalence plot
phyla_number <- taxonomy %>%
  select(phylum) %>%
  count(phylum, name = "total")

barplots_total_ds <- all_result %>%
  filter(!grepl("_other", system)) %>%
  select(genome) %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  unique() %>%
  inner_join(taxonomy, by = "genome") %>%
  select(phylum) %>%
  count(phylum) %>%
  inner_join(phyla_number, by = "phylum") %>%
  mutate(phylum = str_remove(phylum, "_.*")) %>%
  group_by(phylum) %>% 
  summarise(sum_n = sum(n),
            sum_total = sum(total)) %>% 
  filter(sum_n > 10) %>% 
  mutate(percent = 100 * (sum_n/sum_total)) %>% 
  ggplot(aes(
    x = percent,
    y = fct_reorder(phylum, percent),
  )) +
  geom_col() +
  geom_text(aes(x = percent - 8,
                y = phylum,
                label = sum_n)) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0),
  ) +
  labs(
    x = "Genomes with defence system (%)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

barplots_total_ds %>% 
  summarise(mean = mean(percent))



# Number of systems -----------------------------------------------------------

# Order of phyla
phyla_order <- c("Verrucomicrobiota", "Fibrobacterota", "Methanobacteriota",
                 "Bacteroidota", "Actinomycetota", "Cyanobacteriota",
                 "Bacillota", "Elusimicrobiota", "Spirochaetota",
                 "Pseudomonadota")

number_system <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  inner_join(phyla_number, by = "phylum") %>%
  mutate(n = replace_na(n, 0),
         phylum = str_remove(phylum, "_.*"),
         phylum = factor(phylum, levels = phyla_order)) %>%
  filter(total > 10) %>%
  select(phylum, n) %>%
  ggplot(aes(
    x = n,
    y = fct_rev(phylum) 
  )) +
  geom_boxplot(
    width = .3, size = 0.3, outlier.size = 0.5
  ) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 36),
    breaks = seq(0, 36, 3)
  ) +
  labs(
    x = "Defence systems per genome",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_line(linetype = 2, linewidth = 0.3),
    axis.text.y = element_blank()
  )


# Density plot ----------------------------------------------------------------
number_system_density <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  right_join(taxonomy, by = "genome") %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(quality_genomes, by = "genome") %>%  
  select(phylum, n, Genome_Size) %>%
  filter(phylum %in% phyla_order) %>% 
  mutate(kbp = Genome_Size/1000,
         density = n/kbp,
         density_1000 = density*1000,
         phylum = str_remove(phylum, "_.*"),
         phylum = factor(phylum, levels = phyla_order)) %>% 
  
  ggplot(aes(
    y = fct_rev(phylum),
    x = density_1000)) +
  geom_boxplot(
    width = .3, size = 0.3, outlier.size = 0.5
  ) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 13),
    breaks = seq(0, 13, 1),
    oob = scales::squish
  ) +
  labs(
    x = bquote("Defence systems / kbp"~(x10^-3)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_line(linetype = 2, linewidth = 0.3),
    axis.text.y = element_blank()
  )



# Compose plot

compose_plot <- barplots_total_ds + number_system + number_system_density + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14))

ggsave(compose_plot, file = "defence_rummen/new_plots/phylum_defence.png",
       width = 10, height = 5, dpi = 300)


x <- number_system %>% 
  group_by(phylum) %>% 
  summarise(mean = mean(n), .groups = "drop")
