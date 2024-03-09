# set up working enviroment ----------------------------------------------------
pacman::p_load(tidyverse,
               here,
               patchwork)

metadata <- read_tsv("conservation_protocol/rawdata/metadata.txt")

proteins <- read_tsv("conservation_protocol/rawdata/functions_report.tsv")

# clean metadata ---------------------------------------------------------------
clean_meta <- metadata %>%
  separate(treatment,
           into = c("sampling", "treatment", "incubation", "substrate"),
           sep = "-") %>% 
  mutate(incubation = str_replace(incubation, "1", "Pre-treated"),
         incubation = str_replace(incubation, "0[a-b]", "No incubation"),
         incubation = factor(incubation,
                             levels = c("No incubation", "Pre-treated", "8", "24", "72")),
         substrate = if_else(is.na(substrate), "None", substrate),
         treatment = factor(treatment,
                            levels = c("C", "W", "F", "FN", "FD", "FD_noN", "FD_noN3Wo")))

# Filter proteins --------------------------------------------------------------

ko_names <- c(
  "K00200" = "fwdA",
  "K00202" = "fwdC",
  "K00203" = "fwdD",
  "K00319" = "mtd",
  "K00320" = "mer",
  "K00399" = "mcrA",
  "K00401" = "mcrB",
  "K00402" = "mcrG",
  "K00577" = "mtrA",
  "K01499" = "mch",
  "K03388" = "hdrA2"
)

proteins %>% 
  filter(Protein_ID == 1,
         Name %in% protein_names) %>% 
  filter(grepl("M00567", KEGG_Module)) %>% 
  select(Name, KEGG_ko, starts_with("Intensity ")) %>% 
  pivot_longer(-c(Name, KEGG_ko), names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample),
         KEGG_ko = str_remove(KEGG_ko, "ko:")) %>% 
  group_by(sample, KEGG_ko) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  mutate(log_intensity = log2(sum_intensity + 1)) %>% 
  inner_join(clean_meta, by = "sample") %>% 
  #filter(K00399) %>% 
  ggplot(aes(x = treatment,
             y = log_intensity)) +
  geom_boxplot(outlier.colour = "white") +
  geom_jitter(aes(color = incubation,
                  shape = substrate)) +
  facet_wrap(~KEGG_ko, scales = "free_y",
             labeller = labeller(KEGG_ko = ko_names)) +
  labs(x = NULL,
       y = "Abundance (Log2 Intensity") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,angle = 45, hjust = 1))

ggsave(filename = "conservation_protocol/plots/methane_proteins.png",
       width = 10, height = 6, dpi = 450)
# Cazyme -----------------------------------------------------------------------



proteins %>% 
  filter(Protein_ID == 1,
         Name %in% protein_names,
         CAZy != "-") %>% 
  select(Name, CAZy, starts_with("Intensity ")) %>% 
  pivot_longer(-c(Name, CAZy), names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  group_by(sample, CAZy) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  mutate(log_intensity = log2(sum_intensity + 1)) %>% 
  inner_join(clean_meta, by = "sample") %>% 
  filter(CAZy == "GH57") %>% 
  ggplot(aes(x = treatment,
             y = log_intensity)) +
  geom_boxplot(outlier.colour = "white") +
  geom_jitter(aes(color = incubation,
                  shape = substrate)) +
  facet_wrap(~CAZy, scales = "free_y") +
  labs(x = NULL,
       y = "Abundance (Log2 Intensity") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,angle = 45, hjust = 1))
