# set up working space ---------------------------------------------------------
set.seed(150988)
library(tidyverse)
library(ggrepel)
library(here)
library(UpSetR)
library(remotes)

taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  rename(genome=Genome)

clean_metadata <- readRDS("conservation_protocol/output/clean_meta.rds")

gtdb <- read_tsv("conservation_protocol/rawdata/new_gtdbtk.bac120.tsv") %>% 
  rename(Genome=user_genome) %>%   
  rename_all(tolower) %>%
  mutate(classification=str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into=c("kingdom", "phylum", "class", "order", "family", "genus",
                  "species"),
           sep=";")

# Taxonomy analysis ------------------------------------------------------------
taxonomy %>%
  inner_join(gtdb, by = "genome") %>% 
  select(phylum, starts_with("intensity")) %>% 
  pivot_longer(-phylum, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, phylum) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(rel_abun = 100 * (sum_intensity/sum(sum_intensity))) %>% 
  inner_join(clean_metadata, by = "sample") %>% 
  #filter(phylum == "Spirochaetota") %>% 
  ggplot(aes(x = treatment,
             y = rel_abun)) +
  geom_boxplot(outlier.color = "white", width = 0.5) +
  geom_jitter(aes(color = substrate),width = 0.3, size = 1)+
  facet_wrap(~phylum, scales = "free_y", ncol = 5,)+
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.02)) +
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(filename = "conservation_protocol/plots/all_tax_phylum.png",
       width = 10, height = 8)
# boxplot incubation -----------------------------------------------------------

taxonomy %>%
  inner_join(gtdb, by = "genome") %>% 
  select(phylum, starts_with("intensity")) %>% 
  pivot_longer(-phylum, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, phylum) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(rel_abun = 100 * (sum_intensity/sum(sum_intensity))) %>% 
  inner_join(clean_meta, by = "sample") %>% 
  filter(phylum == "Methanobacteriota") %>% 
  ggplot(aes(x = incubation,
             y = rel_abun)) +
  geom_boxplot(width = 0.5, outlier.color = "white")+
  geom_jitter(aes(color = incubation), width = 0.3) +
  facet_wrap(~treatment, scales = "free_y",) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1 , 0.1)) +
  labs(x = NULL,
       y = "Relative abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
        strip.background = element_blank(),
        axis.line.x = element_line(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.02)) +
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(filename = "conservation_protocol/plots/methano_boxplot_incu.png",
       width = 5, height = 5, dpi = 450)

# correlation ---------------
  
correlation_df <- taxonomy %>%
  inner_join(gtdb, by = "genome") %>% 
  select(phylum, starts_with("intensity")) %>% 
  pivot_longer(-phylum, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, phylum) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop")  %>% 
  inner_join(clean_metadata, by = "sample")

compared_t <- 

  correlation_plot <- function(compared_t, compared_incubation) {
    plot <<- correlation_df %>% 
      filter(treatment == "C" | treatment == {{compared_t}},
             incubation == compared_incubation) %>% 
      group_by(treatment, incubation, phylum) %>% 
      summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
      mutate(log_intensity = log2(mean_intensity + 1)) %>% 
      pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity)  %>% ggplot(aes(x = C,
                                                                                                         y = compared_t))
    return(plot) 
  }
  
correlation_plot("FR", "8")
  
correlation_plot <- function(compared_t, compared_incubation) {
  plot <- correlation_df %>% 
  filter(treatment == "C" | treatment == compared_t,
         incubation == compared_incubation) %>% 
  group_by(treatment, incubation, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = C,
             y = {{compared_t}})) +
  geom_abline(linetype = "dashed", color = "grey") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "W-8h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()
  return(plot)
}

correlation_plot("FD", "8")
plot

C_W_24 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "W",
         incubation == "No incubation" |incubation == "24") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = W,
             y = C)) +
  geom_abline(linetype = "dashed", color = "grey") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "W- 24h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_F_8 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "F",
         incubation == "No incubation" |incubation == "8") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = F,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "F-8h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_F_24 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "F",
         incubation == "No incubation" |incubation == "24") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = F,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "F-24h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_FN_8 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "FN",
         incubation == "No incubation" |incubation == "8") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = FN,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "FN-8h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_FN_24 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "FN",
         incubation == "No incubation" |incubation == "24") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = FN,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "FN-24h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_FD_8 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "FD",
         incubation == "No incubation" |incubation == "8") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = FD,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "FD-8h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

C_FD_24 <- correlation_df %>% 
  filter(treatment == "C" | treatment == "FD",
         incubation == "No incubation" |incubation == "24") %>% 
  group_by(treatment, phylum) %>% 
  summarise(mean_intensity = median(sum_intensity), .groups = "drop") %>%
  mutate(log_intensity = log2(mean_intensity + 1)) %>% 
  pivot_wider(id_cols = phylum, names_from = treatment, values_from = log_intensity) %>% 
  ggplot(aes(x = FD,
             y = C)) +
  geom_abline(linetype = "dashed") +
  geom_point()  +
  geom_text_repel(aes(label=phylum), size = 3,
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 40),
                     breaks = seq(0, 40, 5))+
  labs(x = "FD-24h incubation (Log2 intensity)",
       y = "C (Log2 intensity)") +
  theme_classic()

regression_plot <- (C_W_8 + C_W_24) / (C_F_8 + C_F_24) / (C_FN_8 + C_FN_24) / (C_FD_8 +C_FD_24)




ggsave(regression_plot,
       filename="conservation_protocol/plots/regression.png",
       dpi = 450, width = 10, height = 10)


# Upset plot -------------------------------------------------------------------

number_samples <- tibble(treatment = c("C", "W", "F", "FN", 
                                       "FD", "FD_noN", "FD_noN3Wo"),
                         number = c(8, 21, 19, 21, 21 ,20, 6))

x <- taxonomy %>% 
  select(genome, contains("Intensity")) %>%
  inner_join(gtdb, by ="genome") %>% 
  select(family, contains("Intensity")) %>% 
  pivot_longer(-family, names_to = "sample", values_to = "intensity") %>%
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, family) %>% 
  summarise(sum_inte = sum(intensity), .groups = "drop") %>% 
  mutate(sum_inte = if_else(sum_inte > 0, 1, 0)) %>% 
  right_join(clean_meta, by = "sample") %>% 
  group_by(family, treatment) %>% 
  summarise(sum_inte = sum(sum_inte), .groups = "drop") %>% 
  inner_join(number_samples, by = "treatment") %>% 
  mutate(relative = sum_inte/number) %>% 
  mutate(present = if_else(relative >= 0.7, 1, 0)) %>% 
  pivot_wider(id_cols = family, names_from = treatment, values_from = present) %>% 
  drop_na()


test <- read_tsv("conservation_protocol/rawdata/test.txt")

  
m <- make_comb_mat(x, mode = "distinct")

set_size(m)

# plot the upsetplot -----------------------------------------------------------
png(filename = "conservation_protocol/plots/upset.png",
    width = 15, height = 8, units = "cm", res=1200)

UpSet(m,
              top_annotation = upset_top_annotation(m, add_numbers = TRUE),
              right_annotation = upset_right_annotation(m, add_numbers = TRUE))
dev.off()
