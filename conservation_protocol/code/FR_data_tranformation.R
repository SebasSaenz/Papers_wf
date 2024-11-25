# set envitoment ---------------------------------------------------------------
#transform data from metalab-Mag to calculate Functional redundancy
# based on https://github.com/yvonnelee1988/Metaproteome_FRp

# Set working directory --------------------------------------------------------


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(here)

# Load data --------------------------------------------------------------------
functions <- read_tsv("conservation_protocol/rawdata/functions_report.tsv")

taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  select(Genome, GTDB_Name, GTDB_Genus) %>% 
  mutate(Genus = if_else(is.na(GTDB_Genus), GTDB_Name, GTDB_Genus)) %>% 
  select(Genome, Genus)

# Reshape data
x <- functions %>%
  filter(Protein_ID == 1) %>% 
  separate(Name,
           into = "Genome",
           sep = "_",
           remove = FALSE) %>% 
  separate(KEGG_ko,
           into = "KEGG",
           sep = ",") %>% 
  separate(`COG accession`,
           into = "COG",
           sep = ",") %>%
  inner_join(taxonomy, by = "Genome") %>% 
  mutate(KEGG_COG = if_else(KEGG=="-", COG, KEGG),
         KEGG_COG = str_remove(KEGG_COG, "ko:")) %>% 
  select(Name, Genus, KEGG_COG, contains("Intensity ")) %>% 
  rename('Protein IDs'=Name)

names(x) <- gsub(pattern = "Intensity", replacement = " ", x = names(x))

write_tsv(x, "conservation_protocol/rawdata/Full_table_pro_tax_fun_KEGG_COG_for_PCN.tsv")  


taxonomy_fr <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  select(Genome, GTDB_Name, GTDB_Genus, contains("Intensity")) %>% 
  mutate(Genus = if_else(is.na(GTDB_Genus), GTDB_Name, GTDB_Genus)) %>% 
  select(Name=Genus, contains("Intensity")) %>% 
  pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>% 
  group_by(Name, sample) %>% 
  summarise(intensity =sum(intensity), .groups = "drop") %>%
  mutate(sample = str_remove(sample, "Intensity ")) %>% 
  pivot_wider(names_from = sample, values_from = intensity)


write_tsv(taxonomy_fr, "conservation_protocol/rawdata/Taxonomy_table_genus_level_metaproteomics.tsv")
