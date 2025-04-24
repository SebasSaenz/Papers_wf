# set up working space ---------------------------------------------------------
library(here)
library(tidyverse)
library(MicrobiomeStat)
library(broom)
library(patchwork)
library(ComplexHeatmap)


base_color <- c("#d8b365", "#67a9cf", "#2166ac", "#91cf60", "#1a9850")



library(devtools)
#install_github("jokergoo/ComplexHeatmap") Need this for Upsetplots
library(ggrepel)
library(UpSetR)
library(remotes)
library(rstatix)
library(stats)
library(DHARMa)
library(MASS)

# Load files -------------------------------------------------------------------
taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  rename(genome=Genome)

clean_metadata <- readRDS("conservation_protocol/output/clean_meta.rds") %>% 
  filter(treatment != "F")

gtdb <- readRDS("conservation_protocol/output/gtdb.rds") %>% 
  filter(grepl("MGYG", genome))

# Heatmap with significance ----------------------------------------------------

taxonomy_df <- taxonomy %>% 
  inner_join(gtdb, by = "genome") %>% 
  dplyr::select(phylum, starts_with("Intensity")) %>%
  pivot_longer(-phylum, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(phylum, sample) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(rel_abun = sum_intensity/sum(sum_intensity)) %>%
  ungroup() %>%
  dplyr::select(phylum, sample, rel_abun) %>%
  pivot_wider(names_from = sample, values_from = rel_abun) %>%
  column_to_rownames(var = "phylum") %>% 
  as.matrix()

linda_meta <- clean_metadata %>% 
  column_to_rownames(var = "sample")

phylum_linda <- linda(feature.dat = taxonomy_df,
      meta.dat = linda_meta,
      formula = "~treatment+incubation",
      feature.dat.type = c('proportion'),
      prev.filter =0.5)


combined_df <- do.call(rbind, phylum_linda[["output"]]) |> 
  rownames_to_column(var = "test") |>
  separate(test,
           into = c("test", "phylum"),
           sep = "\\.") %>%
  mutate(stars = case_when(padj > 0.05 ~ "",
                           padj <= 0.05 & padj > 0.01 ~ "*",
                           padj <= 0.01 & padj > 0.001 ~ "**",
                           padj <= 0.001 ~ "***"),
         phylum = str_replace(phylum, "_", " "))


hm_treatment <- combined_df %>% 
  filter(grepl("treatment", test)) %>% 
  mutate(test = str_remove(test, "treatment"),
          test = factor(test,
                        levels = c("FR", "FRN", "FD","FDN"))) %>% 
  ggplot(aes(y = phylum, x = test, fill = log2FoldChange)) +
  geom_tile(color = "black")+
  coord_equal(ratio = 0.2 ) +
  geom_text(aes(label = stars), color = "black", size = 4,vjust = 0.8) +
  scale_fill_gradient2(low = "#af8dc3",
                                  mid = "#f7f7f7",
                                  high = "#7fbf7b",
                       limits = c(-3.5, 3)) +
  labs(y = NULL,
       x = NULL,
       title = "Control Vs",
       fill = "FC (Log2)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.6,
                                  size = 10),
        panel.grid = element_blank(),
        legend.title = element_text(),
        axis.text.y = element_text(size = 7))

hm_incubation <- combined_df %>% 
  filter(grepl("incubation", test)) %>% 
  mutate(test = str_remove(test, "incubation"),
         test = factor(test,
                       levels = c(8, 24))) %>% 
  ggplot(aes(y = phylum, x = test, fill = log2FoldChange)) +
  geom_tile(color = "black")+
  coord_equal(ratio = 0.2 ) +
  geom_text(aes(label = stars), color = "black", size = 4,vjust = 0.8) +
  scale_fill_gradient2(low = "#af8dd9",
                       mid = "#f7f7f7",
                       high = "#7fbf7b",
                       limits = c(-3.5, 3)) +
  labs(y = NULL,
       x = NULL,
       title = " Pre-incubation Vs") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.3,
                                  size = 10),
        panel.grid = element_blank(),
        axis.text.y = element_blank())

# Boxplot example --------------------------------------------------------------
taxonomy_bp <- taxonomy %>% 
  inner_join(gtdb, by = "genome") %>% 
  dplyr::select(phylum, starts_with("Intensity")) %>%
  pivot_longer(-phylum, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(phylum, sample) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(rel_abun = sum_intensity/sum(sum_intensity)) %>% 
  ungroup() %>% 
  inner_join(clean_metadata, by = "sample")

set.seed(150988)
bp_treatment <- taxonomy_bp |>
  filter(phylum == "Bacteroidota") |> 
  ggplot(aes(x = treatment,
             y = rel_abun * 100)) +
  geom_boxplot(outlier.colour = "white",
               width = 0.4) +
  geom_jitter(aes(color = treatment),
              width = 0.3,
              alpha = 0.8) +
  scale_y_continuous(limits = c(0, 60),
                     expand = c(0, 0)) +
  scale_color_manual(values = base_color) +
  labs(y = "Relative abundance (%)",
       x = NULL,
       title = "Bacteroidota") +
  theme(panel.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.box.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(size = 10,
                                  face = "italic"))


set.seed(150988)
bp_incubation <- taxonomy_bp |>
  filter(phylum == "Bacteroidota") |> 
  ggplot(aes(x = incubation,
             y = rel_abun * 100)) +
  geom_boxplot(outlier.colour = "white",
               width = 0.3) +
  geom_jitter(aes(color = treatment),
              width = 0.2,
              alpha = 0.8) +
  scale_y_continuous(limits = c(0, 60),
                     expand = c(0, 0)) +
  scale_color_manual(values = base_color) +
  labs(y = "Relative abundance (%)",
       x = NULL,
       title = "Bacteroidota") +
  theme(panel.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.box.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(size = 10,
                                  face = "italic"))
  
# Compose plot
figure2 <- (free(hm_treatment) | free(hm_incubation)) / (bp_treatment | bp_incubation) + 
  plot_layout(guides = 'collect',
              axes = "collect",
              heights = unit(c(7, 5), c('cm', 'cm'))) +
  plot_annotation(tag_levels = 'A')

ggsave(figure2, file="conservation_protocol/plots/figure2.tif", width = 8, height = 6,
       dpi = 600)





# Upset plot -------------------------------------------------------------------

number_samples <- tibble(treatment = c("C", "FR", "FRN", 
                                       "FDN", "FD"),
                         number = c(21, 19, 21, 21 ,20))

counts <- taxonomy %>% 
  select(genome, contains("Intensity")) %>%
  inner_join(gtdb, by ="genome") %>% 
  select(genus, contains("Intensity")) %>% 
  pivot_longer(-genus, names_to = "sample", values_to = "intensity") %>%
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, genus) %>% 
  summarise(sum_inte = sum(intensity), .groups = "drop") %>% 
  mutate(sum_inte = if_else(sum_inte > 0, 1, 0)) %>% 
  right_join(clean_metadata, by = "sample") %>% 
  group_by(genus, treatment) %>% 
  summarise(sum_inte = sum(sum_inte), .groups = "drop") %>% 
  inner_join(number_samples, by = "treatment") %>% 
  mutate(relative = sum_inte/number) %>% 
  mutate(present = if_else(relative >= 0.7, 1, 0)) %>% 
  pivot_wider(id_cols = genus, names_from = treatment, values_from = present) %>% 
  drop_na() 

upset_df_phylum <- counts[rowSums(counts[2:6])>0,] 

m_phylum <- make_comb_mat(upset_df_phylum, mode = "distinct")

set_size(m_phylum)



# plot the upsetplot -----------------------------------------------------------
png(filename = "conservation_protocol/plots/genus_upset.png",
    width = 16, height = 8, units = "cm", res=1200)
UpSet((m_phylum),
      lwd = 2,
      left_annotation = upset_left_annotation(m_phylum, add_numbers = TRUE,),
      top_annotation = upset_top_annotation(m_phylum, add_numbers = TRUE))
dev.off()


# Correlation, this is nice but not useful now --------------------------------

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


# Make table for writing paper

table_abundance <- taxonomy_bp %>% 
  mutate(pencentage = 100 *rel_abun) %>% 
  group_by(genus, treatment, incubation, substrate) %>% 
  summarise(mean_abundance = mean(pencentage),
            sd_abundance = sd(pencentage),
            .groups = "drop")

write_tsv(table_abundance,
          file = "conservation_protocol/output/relative_abundance_table_genus.txt")
