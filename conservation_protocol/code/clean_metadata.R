set.seed(150988)

# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Load data --------------------------------------------------------------------
metadata <- read_tsv("conservation_protocol/rawdata/metadata.txt")

# Clean metadata ---------------------------------------------------------------
clean_meta <- metadata %>% 
  mutate(treatment = str_replace(treatment, "W-0.", "FF-0a")) %>% 
  separate(treatment,
           into = c("sampling", "treatment", "incubation", "substrate"),
           sep = "-") %>% 
  mutate(treatment = case_when(
    treatment == "W" ~ "C",
    treatment == "FD" ~ "FDN",
    treatment == "FD_noN" ~ "FD",
    treatment == "F" ~ "FR",
    treatment == "FN" ~ "FRN",
    treatment == "FF" ~ "F",
    .default = as.character(treatment)),
    incubation = case_when(incubation == "1" ~ "Pre-incubation",
                           incubation == "0a" ~ "No incubation",
                           .default = as.character(incubation)),
    substrate = case_when(substrate == "079" ~ "Wheat grain",
                          substrate == "057" ~ "Rapeseed meal",
                          substrate == "053" ~ "Maize grain",
                          substrate == "Heu" ~ "Hay",
                          substrate == "KF" ~ "Concentrate",
                          is.na(substrate) ~ "None"),
    incubation = factor(incubation,
                        levels = c("No incubation", "Pre-incubation", "8", "24", "72")),
    treatment = factor(treatment,
                       levels = c("F", "C", "FR", "FRN", "FD", "FDN", "FD_noN3Wo")),
    substrate = factor(substrate,
                       levels = c("None", "Hay", "Rapeseed meal", "Wheat grain",
                                  "Maize grain", "Concentrate"))) %>% 
  filter(sample != "49" & sample != "78" & sample != "87",
         treatment != "FD_noN3Wo"
         ) #Remove samples with 72h

# Save clean metadata ----------------------------------------------------------
saveRDS(clean_meta, file = "conservation_protocol/output/clean_meta.rds")



# Clean Peptides file
peptides <- read_tsv("conservation_protocol/rawdata/final_peptides.tsv") %>% 
  select(Proteins, contains("Intensity ")) %>% 
  filter(!grepl("SHEEP", Proteins))

# Save clean metadata ----------------------------------------------------------
saveRDS(peptides, file = "conservation_protocol/output/peptides.rds")
