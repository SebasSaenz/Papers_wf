library(tidyverse)

df_taxo <- read_tsv("rawdata/gtdbtk_all_summary.tsv") %>% 
  select(genome=user_genome, classification) %>% 
  separate(classification,
           into = c("domain", "phylum", "class", "order", "family", "genues", "species"),
           sep = ";") %>% 
  select(genome, phylum) %>% 
  mutate(phylum=str_remove(phylum, "p__"))
  
df_phages <- read_tsv("rawdata/mapping_rumen_phages.txt") %>% 
  select(phage, family)
  
  
df_cas <- read_tsv("rawdata/spacers_match_proviruses.txt") %>% 
  select(c(1:2))

colnames(df_cas) <- c("genome", "phage")

df_clean <- df_cas %>% 
  mutate(genome = str_remove(genome, "_scaffold.*"),
         genome = str_remove(genome, "_scf.*"),
         genome = str_remove(genome, "_[1-9].*")) %>% 
  unique() %>% 
  inner_join(df_taxo, by = "genome") %>%
  left_join(df_phages, by = "phage")


write_tsv(df_clean, file = "rawdata/clean_network_file.txt")
                             