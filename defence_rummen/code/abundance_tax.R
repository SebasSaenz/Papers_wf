library(tidyverse)

# Load files -------------------------------------------------------------------
taxonomy <- read_tsv("rawdata/new_gtdbtk.bac120.tsv")
hq_genomes <- read_tsv("rawdata/hq_genomes.txt")

abundance <- taxonomy %>% 
  select(user_genome, classification) %>%
  rename(genome=user_genome) %>% 
  mutate(classification=str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into=c("kingdom", "phylum", "class", "order", "family", "genus",
                  "species"),
           sep=";") %>% 
  inner_join(hq_genomes, by = "genome") %>% 
  select(phylum) %>%
  mutate(phylum = str_remove(phylum, "_.*")) %>% 
  count(phylum) %>% 
  mutate(rel_abun = 100*(n/sum(n)))
  
