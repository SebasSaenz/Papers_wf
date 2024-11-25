# Johan S. Sáenz, University of Hohenheim
# Calculate the relative abundance of the high-quality genomes by phylum

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(here)

# Load files -------------------------------------------------------------------
taxonomy <- read_tsv("defence_rummen/rawdata/new_gtdbtk.bac120.tsv")
hq_genomes <- read_tsv("defence_rummen/rawdata/hq_genomes.txt")

# Calculate relative abundance
abundance <- taxonomy %>%
  select(user_genome, classification) %>%
  rename(genome = user_genome) %>%
  mutate(classification = str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
    into = c(
      "kingdom", "phylum", "class", "order", "family", "genus",
      "species"
    ),
    sep = ";"
  ) %>%
  inner_join(hq_genomes, by = "genome") %>%
  select(phylum) %>%
  mutate(phylum = str_remove(phylum, "_.*")) %>%
  count(phylum) %>%
  mutate(rel_abun = 100 * (n / sum(n)))
