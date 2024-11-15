# Johan S. Sáenz, University of Hohenheim
# Wrangle data to fit into Cytoskape and create a network

# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Load files -------------------------------------------------------------------
df_taxo <- read_tsv("defence_rummen/rawdata/gtdbtk_all_summary.tsv") %>%
  select(genome = user_genome, classification) %>%
  separate(classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";"
  ) %>%
  select(genome, phylum) %>%
  mutate(phylum = str_remove(phylum, "p__"))

df_phages <- read_tsv("defence_rummen/rawdata/mapping_rumen_phages.txt") %>%
  select(phage, 
         #family
         )


df_cas <- read_tsv("defence_rummen/rawdata/spacers_match_proviruses.txt") %>%
  select(c(1:2))

colnames(df_cas) <- c("genome", "phage")

# Clean data -------------------------------------------------------------------
df_clean <- df_cas %>%
  mutate(
    genome = str_remove(genome, "_scaffold.*"),
    genome = str_remove(genome, "_scf.*"),
    genome = str_remove(genome, "_[1-9].*")
  ) %>%
  unique() %>%
  inner_join(df_taxo, by = "genome") %>%
  left_join(df_phages, by = "phage")

# Save file to use in Cytoskape
write_tsv(df_clean, file = "defence_rummen/rawdata/clean_network_file2.txt")


# Count phages infecting several species
x <- df_clean %>% 
  select(genome, phage, genus) %>% 
  count(phage, genus) %>% 
  count(phage)

write_tsv(df_clean, file = "defence_rummen/rawdata/spacer_match_tax.txt")
