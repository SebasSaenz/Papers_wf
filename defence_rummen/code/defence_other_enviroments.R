# Libraries---------------------------------------------------------------------
library(tidyverse)
library(patchwork)

# Load data--------------------------------------------------------------------
df <- read_tsv("defence_rummen/rawdata/systems_enviroments.txt")

color_val <- c('grey', '#fc8d59', '#ffffbf', '#91bfdb')

factor_order <- c("Rumen", "Human", "Soil", "Marine")

# Prevalence systems-----------------------------------------------------------
rumen_prevalence <- tibble(Enviroment = "Rumen",
                              percent = 89.071
                              )

prevalence <- df %>% 
  mutate(prevalence = if_else(systems >=1, 1, 0)) %>% 
  group_by(Enviroment) %>% 
  summarise(percent = 100 * (sum(prevalence)/n())) %>%
  rbind(rumen_prevalence) %>% 
  mutate(Enviroment = factor(Enviroment, levels = factor_order)) %>% 
  ggplot(aes(x = Enviroment,
             y = percent,
             fill = Enviroment)) + 
  geom_col() +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     breaks = seq(0, 100, 10),
                     oob = scales::squish) +
  scale_fill_manual(values = color_val) +
  labs(x = NULL,
       y = "Genomes with defence systems (%)") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )
  

# Number of defence systems ---------------------------------------------------


rumen_n_system <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  mutate(n = replace_na(n, 0)) %>% 
  mutate(Enviroment = if_else(is.numeric(n), "Rumen", "")) %>% 
  rename(systems=n) %>% 
  select(systems, Enviroment)


number_ds <- df %>% 
  select(systems, Enviroment) %>% 
  rbind(rumen_n_system) %>% 
  mutate(Enviroment = factor(Enviroment, levels = factor_order)) %>% 
  ggplot(aes( x = Enviroment,
              y = systems,
              fill = Enviroment)) +
  geom_violin(aes(
    fill = Enviroment, fill = after_scale(colorspace::lighten(fill, .5))),
              size = 1.2, color = NA,
              bw = 0.5) +
  geom_boxplot(width = .1, size = 0.3, outlier.size = 0.5) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 36),
                     breaks = seq(0, 36, 3),
                     oob = scales::squish) +
  scale_fill_manual(values = color_val) +
  labs(x = NULL,
       y = "Defence systems per genome") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )


# Defence system density ------------------------------------------------------

rumen_density <- all_result %>% # all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl("_other", system)) %>%
  unique() %>%
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>%
  count(genome) %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(quality_genomes, by = "genome") %>%  
  #select(kingdom, n, Genome_Size) %>%
  mutate(kbp = Genome_Size/1000,
         density = n/kbp) %>% 
  mutate(Enviroment = if_else(is.numeric(density), "Rumen", "")) %>% 
  rename(Density=density) %>% 
  select(Density, Enviroment)

density <- df %>%
  select(Density, Enviroment) %>% 
  rbind(rumen_density) %>% 
  mutate(Density_1000 = Density*1000,
         Enviroment = factor(Enviroment, levels = factor_order)) %>% 
  ggplot(aes( x = Enviroment,
              y = Density_1000,
              fill = Enviroment)) +
  geom_violin(aes(
    fill = Enviroment, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, color = NA,
    bw = 0.5) +
  geom_boxplot(width = .1, size = 0.3, outlier.size = 0.5) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 16),
                     breaks = seq(0, 16, 2),
                     oob = scales::squish) +
  scale_fill_manual(values = color_val) +
  labs(x = NULL,
       y = bquote("Defence system / kbp"~(x10^-3))) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3)
  )

compose_plot <- prevalence + number_ds + density 

ggsave(compose_plot, file = "defence_rummen/new_plots/defence_enviroments.png",
       width = 10, height = 5, dpi = 300)

