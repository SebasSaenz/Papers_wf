library(tidyverse)

checkm <- read_tsv("defence_rummen/rawdata/quality_report_genomes.tsv") %>% 
  rename(genome=Name, contamination=Contamination, completeness=Completeness)

taxonomy <- read_tsv("defence_rummen/rawdata/new_gtdbtk.bac120.tsv") %>% 
  rename(genome = user_genome)

hq_genomes <- read_tsv("defence_rummen/rawdata/hq_genomes.txt")

taxonomy <- taxonomy %>%
  select(genome, classification) %>%
  rename_all(tolower) %>%
  mutate(classification=str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into=c("kingdom", "phylum", "class", "order", "family", "genus",
                  "species"),
           sep=";") %>% 
  inner_join(hq_genomes)

sysfamily_list <- c("Abi", "AVAST", "CBASS", "CRISPR-Cas", "DRT", "gabija",
                    "Phosphorothioation", "Pycsar", "RM", "wadjet")

#  Limimorpha ------------------------------------------------------------------
Limimorpha_list <- taxonomy %>% 
  select(genome, genus) %>% 
  filter(genus == "Limimorpha") %>% 
  select(genome)

checkm_limimorpha <- inner_join(checkm, Limimorpha_list) %>% 
  select(genome, completeness, contamination) %>% 
  mutate(genome = str_replace(genome, "$", ".fasta"))

write_tsv(Limimorpha_list, file = "rawdata/Limimorpha_list.txt")
write_tsv(checkm_limimorpha, file = "rawdata/checkm_limimorpha.txt")

tree_annotation_limi <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl('_other', system)) %>% 
  unique() %>% 
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>% 
  select(system, genome) %>% 
  inner_join(Limimorpha_list, by = "genome") %>% 
  inner_join(list_system, by = "system") %>% 
  select(system_family, genome) %>% 
  count(system_family, genome) %>%
  complete(system_family, genome, fill = list(n=0)) %>%
  #filter(system_family %in% sysfamily_list) %>%
  mutate(n=if_else(n>=1, 1,0)) %>% 
  pivot_wider(genome, values_from = n, names_from = system_family)

write_tsv(tree_annotation_limi, file = "rawdata/tree_annotation_limi.txt")

limi <- tree_annotation_limi %>% 
  group_by(system_family) %>% 
  summarise(sum = 100 * sum(n)/41)

#  Prevotella ------------------------------------------------------------------
Prevotella_list <- taxonomy %>% 
  select(genome, genus) %>% 
  filter(genus == "Prevotella") %>% 
  select(genome)

checkm_Prevotella <- inner_join(checkm, Prevotella_list) %>% 
  select(genome, completeness, contamination) %>% 
  mutate(genome = str_replace(genome, "$", ".fasta"))

write_tsv(checkm_Prevotella, file = "rawdata/checkm_Prevotella.txt")
write_tsv(Prevotella_list, file = "rawdata/Prevotella_list.txt")

tree_annotation_prevo <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl('_other', system)) %>% 
  unique() %>% 
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>% 
  select(system, genome) %>% 
  inner_join(Prevotella_list, by = "genome") %>% 
  inner_join(list_system, by = "system") %>% 
  select(system_family, genome) %>% 
  count(system_family, genome) %>%
  complete(system_family, genome, fill = list(n=0)) %>%
  #filter(system_family %in% sysfamily_list) %>%
  mutate(n=if_else(n>=1, 1,0)) %>% 
  pivot_wider(genome, values_from = n, names_from = system_family)

write_tsv(tree_annotation_limi, file = "rawdata/tree_annotation_prevo.txt")

prevo <- tree_annotation_prevo %>% 
  group_by(system_family) %>% 
  summarise(sum = 100 * sum(n)/258)
#  Butirivibrio ------------------------------------------------------------------
Butirivibrio_list <- taxonomy %>% 
  select(genome, genus) %>% 
  filter(genus == "Butyrivibrio") %>% 
  select(genome)

checkm_Butyrivibrio <- inner_join(checkm, Butirivibrio_list) %>% 
  select(genome, completeness, contamination) %>% 
  mutate(genome = str_replace(genome, "$", ".fasta"))

write_tsv(Butirivibrio_list, file = "rawdata/Butirivibrio_list.txt")
write_tsv(checkm_Butyrivibrio, file = "rawdata/checkm_Butyrivibrio.txt")


tree_annotation_buty <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl('_other', system)) %>% 
  unique() %>% 
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>% 
  select(system, genome) %>% 
  inner_join(Butirivibrio_list, by = "genome") %>% 
  inner_join(list_system, by = "system") %>% 
  select(system_family, genome) %>% 
  count(system_family, genome) %>%
  complete(system_family, genome, fill = list(n=0)) %>%
  #filter(system_family %in% sysfamily_list) %>%
  mutate(n=if_else(n>=1, 1,0)) %>% 
  pivot_wider(genome, values_from = n, names_from = system_family)

write_tsv(tree_annotation_buty, file = "rawdata/tree_annotation_buty.txt")

buty <- tree_annotation_buty %>% 
  group_by(system_family) %>% 
  summarise(sum = 100 * sum(n)/55)
#  Ruminococcus ------------------------------------------------------------------
Ruminococcus_list <- taxonomy %>% 
  select(genome, genus) %>% 
  filter(grepl("Ruminococcus", genus)) %>% 
  select(genome)

checkm_Ruminococcus <- inner_join(checkm, Ruminococcus_list) %>% 
  select(genome, completeness, contamination) %>% 
  mutate(genome = str_replace(genome, "$", ".fasta"))

write_tsv(Ruminococcus_list, file = "rawdata/Ruminococcus_list.txt")
write_tsv(checkm_Ruminococcus, file = "rawdata/checkm_Ruminococcus.txt")


Ruminococcus_list <- all_result %>%
  select(genome, system.number, system) %>%
  filter(!grepl('_other', system)) %>% 
  unique() %>% 
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>% 
  select(system, genome) %>% 
  inner_join(Ruminococcus_list, by = "genome") %>% 
  inner_join(list_system, by = "system") %>% 
  select(system_family, genome) %>% 
  count(system_family, genome) %>%
  complete(system_family, genome, fill = list(n=0)) %>%
  #filter(system_family %in% sysfamily_list) %>%
  mutate(n=if_else(n>=1, 1,0)) %>% 
  pivot_wider(id_cols = genome, values_from = n, names_from = system_family)

write_tsv(Ruminococcus_list, file = "rawdata/tree_annotation_rummi.txt")

rummi <- Ruminococcus_list %>% 
  group_by(system_family) %>% 
  summarise(sum = 100 * sum(n)/138)
