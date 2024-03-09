setwd("/Users/sebastiansaenz/Documents/defence_rummen")

library(tidyverse)

quality_genomes <- read_tsv("rawdata/quality_report_genomes.tsv") %>% 
  rename(genome=Name)

taxonomy <- read_tsv("rawdata/gtdbtk_all_summary.tsv") %>% 
  rename(genome = user_genome)

hq_genomes <- read_tsv("rawdata/hq_genomes.txt")

code_color <- read_tsv("rawdata/phylum_colorcode.txt")

taxonomy <- taxonomy %>%
  select(genome, classification) %>%
  rename_all(tolower) %>%
  mutate(classification=str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into=c("kingdom", "phylum", "class", "order", "family", "genus",
                  "species"),
           sep=";") %>% 
  inner_join(hq_genomes) %>% 
  select(genome, phylum) %>% 
  inner_join(code_color, by = "phylum")

write_tsv(taxonomy, file = "rawdata/tree_phylum_color.txt")

#---- Counting genomes

test <- taxonomy %>% 
  inner_join(hq_genomes, by= "genome") %>% 
  select(genome, classification) %>% 
  separate(classification,
           into = c("kingdom", "phylum", "class", "order", "family", "genus",
                    "species"),
           sep = ";") %>% 
  select(phylum) %>% 
  mutate(phylum = str_remove(phylum, "p__"),
         phylum = str_remove(phylum, "_[A-Z]")
         ) %>% 
  count(phylum) %>% 
  mutate(abundance = 100*(n/sum(n)))






