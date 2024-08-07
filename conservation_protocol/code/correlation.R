# Load library------------------------------------------------------------------
library(tidyverse)
library(broom)
library(here)
library(gt)

base_color <- c("#d8b365", "#67a9cf", "#2166ac", "#91cf60", "#1a9850")


# Load data---------------------------------------------------------------------
df_gas <- readxl::read_xlsx("conservation_protocol/rawdata/Sample list Correlations_GP-MO.xlsx") %>% 
  rename(sample = Sample, gas_production = `Gas production (ml/ 200 mg DM)`)

fr_df <-  read_csv("conservation_protocol/rawdata/Record_FR_PCN.csv") %>%
  rename(sample = Sample) %>% 
  mutate(sample = str_remove(sample, "Nr0"),
         sample = str_remove(sample, "Nr"),
         sample = as.numeric(sample))

metadata <- readRDS("conservation_protocol/output/clean_meta.rds")

taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  rename(genome=Genome)

gtdb <- readRDS("conservation_protocol/output/gtdb.rds") %>% 
  filter(grepl("MGYG", genome))


# Correlation functional redundancy
fr_df %>%
  inner_join(metadata, by = "sample") %>% 
  inner_join(df_gas, by = "sample") %>% 
  select(sample, treatment, nFR, gas_production, incubation) %>% 
  ggplot(aes(x = nFR,
             y = gas_production))+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  geom_point(aes(color = treatment, shape = incubation)) +
  scale_color_manual(values =  base_color) +
  labs(x = "Fucntional redundancy (nFR)",
       y = "Gas production (ml/200mg DM)") +
  theme_classic()

ggsave(filename = "conservation_protocol/plots/correlation_fr.png",
       width = 6,
       height = 4,
       dpi = 400)

cor.test(cor_df$nFR, cor_df$gas_production, method = "pearson")


# Correlation phylun abundance 
taxonomy_df <- taxonomy %>%  
inner_join(gtdb, by = "genome")  %>%  
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
  mutate(rel_abun = 100 *(sum_intensity/sum(sum_intensity))) %>% 
  ungroup() %>% 
  inner_join(df_gas, by ="sample")




core_phylum <- taxonomy_df %>% 
  nest(data = -phylum) %>% 
  mutate(test = map(.x=data, ~cor.test(.x$rel_abun, .x$gas_production, method = "pearson") %>% tidy())) %>% 
  unnest(test) %>% 
  mutate(padjust = p.adjust(p.value), method = "BH") %>% 
  select(phylum, estimate, conf.low, conf.high , padjust)


# Save correlation table
core_phylum %>% 
  gt() %>% 
  gtsave(filename = "conservation_protocol/output/taxonomy_correlation.docx")

taxonomy_df %>% 
  filter(phylum == "Bacillota_C") %>% 
  ggplot(aes(x = rel_abun,
             y = gas_production))+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  geom_point() +
  scale_color_manual(values =  base_color) +
  theme_classic()
  
