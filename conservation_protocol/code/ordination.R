################################################################################
#                           Code by Johan S. Sáenz                             #
#                           University of Hohenheim                            #
#                                                                              #
################################################################################

# Load data and libraries ------------------------------------------------------
set.seed(150988)

library(tidyverse)
library(vegan)
library(here)
library(patchwork)

taxonomy <- read_tsv("conservation_protocol/rawdata/Genome.tsv") %>% 
  rename(genome=Genome)

metadata <- read_tsv("conservation_protocol/rawdata/metadata.txt")

proteins <- read_tsv("conservation_protocol/rawdata/functions_report.tsv")


F='#8c510a'
C='#d8b365'
FDN='#2166ac'
FD = '#67a9cf'
FRN = '#1a9850'
FR = '#91cf60'


base_color <- c('#8c510a', '#d8b365','#67a9cf', '#2166ac','#91cf60','#1a9850')
# clean metadata ---------------------------------------------------------------
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
                            levels = c("F", "C", "FR", "FRN", "FD", "FDN", "FD_noN3Wo"))) %>% 
  filter(sample != "49" & sample != "78" & sample != "87")


#-----------
df_taxonomy <- taxonomy %>% 
  select(genome, starts_with("Intensity ")) %>% 
  pivot_longer(-genome, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST")) %>% 
  mutate(rel_abun = 100*(intensity/sum(intensity)),
         .by = sample) %>% 
  pivot_wider(id_cols = genome, names_from = sample, values_from = rel_abun)

# Create matrix-----------------------------------------------------------------
matrix_taxonomy <- df_taxonomy[c(2:ncol(df_taxonomy))] %>% 
  t()

dist <- vegdist(matrix_taxonomy,method = "bray")

adonis2(as.dist(dist)~clean_meta$treatment * clean_meta$substrate, permutations = 1000)


# Calculate nMDS
nmds1_taxonomy <- metaMDS(matrix_taxonomy,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)

nmds_best_taxonomy <- metaMDS(matrix_taxonomy, 
                              distance = "bray", #find best nmds
                              try = 20, 
                     trymax = 100, 
                     maxit = 1000, 
                     k = 3,
                     previous.best = nmds1_taxonomy)

nmds_best_taxonomy$stress          

# Extract ndms data points -----------------------------------------------------
data_scores_tax <- as.data.frame(scores(nmds_best_taxonomy, display=c("sites")))

data_scores_tax$sample <- as.numeric(row.names(data_scores_tax))

data_nmds_tax <- data_scores_tax  %>% 
  left_join(clean_meta, by = "sample") %>% 
  filter(treatment != "FD_noN3Wo")



tax_plot <- data_nmds_tax %>% 
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             color = treatment,
             shape = incubation)) +
  geom_point(alpha = 0.9,
             size = 3
             ) +
  scale_color_manual(values = base_color) +
  labs(title = "Taxonomy") +
  theme_classic()+
  theme(legend.position = "none")

tax_plot_feed <- data_nmds_tax %>% 
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             color = treatment,
             shape = substrate)) +
  geom_point(alpha = 0.9,
             size = 3
  ) +
  scale_color_manual(values = base_color) +
  labs(title = "Taxonomy") +
  theme_classic()+
  theme(legend.position = "none")


#------------------------------------------------------------------------------

df_proteins <- proteins %>% 
  filter(Protein_ID == 1) %>%
  select(Name, starts_with("Intensity ")) %>% 
  pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>% 
    mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
           sample = str_remove(sample, "_WDH"),
           sample = str_remove(sample, "_STneu"),
           sample = str_remove(sample, "_wdh"),
           sample = str_remove(sample, "_neuST")) %>% 
    mutate(rel_abun = 100*(intensity/sum(intensity)),
           .by = sample) %>%
    pivot_wider(id_cols = Name, names_from = sample, values_from = rel_abun)

matrix_proteins <- df_proteins[c(2:ncol(df_proteins))] %>% 
  t()

# Calculate nMDS
nmds1_proteins <- metaMDS(matrix_proteins,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)

nmds_best_proteins <- metaMDS(matrix_proteins, distance = "bray", #find best nmds
                     try = 20, 
                     trymax = 100, 
                     maxit = 1000, 
                     k = 3,
                     previous.best = nmds1_proteins)

nmds_best_proteins$stress          

# Extract ndms data points -----------------------------------------------------
data_scores_proteins <- as.data.frame(scores(nmds_best_proteins, display=c("sites")))

data_scores_proteins$sample <- as.numeric(row.names(data_scores_proteins))

data_nmds_proteins <- data_scores_proteins  %>% 
  left_join(clean_meta, by = "sample") %>% 
  filter(treatment != "FD_noN3Wo")



proteins_plot <- data_nmds_proteins %>% 
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             color = treatment,
             shape = incubation)) +
  geom_point(alpha = 0.7,
             size = 3) +
  scale_color_manual(values = base_color) +
  labs(title = "Functions") +
  theme_classic()+
  theme(legend.title = element_blank()) 

proteins_plot_feed <- data_nmds_proteins %>% 
  ggplot(aes(x = NMDS1,
             y = NMDS2,
             color = treatment,
             shape = substrate)) +
  geom_point(alpha = 0.7,
             size = 3) +
  scale_color_manual(values = base_color) +
  labs(title = "Functions") +
  theme_classic()+
  theme(legend.title = element_blank())

ordination <- (tax_plot +proteins_plot) /(tax_plot_feed + proteins_plot_feed)

 ggsave(ordination, file = "conservation_protocol/plots/ordination.png", 
       width = 12, height = 8)


# Ordination core proteins -----------------------------------------------------

protein_names <- proteins %>% 
   filter(Protein_ID == 1) %>%
   select(Name, starts_with("Intensity ")) %>% 
   pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>% 
   mutate(intensity_presence = if_else(intensity > 0, 1, 0),
          sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
          sample = str_remove(sample, "_WDH"),
          sample = str_remove(sample, "_STneu"),
          sample = str_remove(sample, "_wdh"),
          sample = str_remove(sample, "_neuST"),
          sample = as.numeric(sample)) %>% 
   group_by(Name) %>% 
   summarise(sum_presence = sum(intensity_presence)) %>% 
   filter(sum_presence >80) %>% 
   pull(Name)

 
 protein_filter <- proteins %>% 
   filter(Protein_ID == 1) %>%
   select(Name, starts_with("Intensity ")) %>% 
   pivot_longer(-Name, names_to = "sample", values_to = "intensity") %>% 
   filter(Name %in% protein_names) %>% 
   mutate(sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
                                              sample = str_remove(sample, "_WDH"),
                                              sample = str_remove(sample, "_STneu"),
                                              sample = str_remove(sample, "_wdh"),
                                              sample = str_remove(sample, "_neuST")) %>% 
   mutate(rel_abun = 100*(intensity/sum(intensity)),
          .by = sample) %>%
   pivot_wider(id_cols = Name, names_from = sample, values_from = rel_abun)
 
 matrix_proteins_filtered <- protein_filter[c(2:ncol(protein_filter))] %>% 
   t()
 
 # Calculate nMDS
 nmds1_proteins_filtered <- metaMDS(matrix_proteins_filtered,  #perform nmds
                           distance = "bray", 
                           try = 20, 
                           trymax = 100, 
                           maxit = 1000, 
                           k = 3)
 
 nmds_best_proteins_filtered <- metaMDS(matrix_proteins_filtered, distance = "bray", #find best nmds
                               try = 20, 
                               trymax = 100, 
                               maxit = 1000, 
                               k = 3,
                               previous.best = nmds1_proteins_filtered)
 
 nmds_best_proteins_filtered$stress          
 
 # Extract ndms data points -----------------------------------------------------
 data_scores_proteins_filtered <- as.data.frame(scores(nmds_best_proteins_filtered, display=c("sites")))
 
 data_scores_proteins_filtered$sample <- as.numeric(row.names(data_scores_proteins_filtered))
 
 data_nmds_proteins_filtered <- data_scores_proteins_filtered  %>% 
   left_join(clean_meta, by = "sample") %>% 
   filter(treatment != "FD_noN3Wo")
 
 
centroid <-  data_nmds_proteins_filtered %>% 
   group_by(treatment) %>% 
   summarise(NMDS1 = mean(NMDS1),
             NMDS2 = mean(NMDS2))
 
 prot_70p_incubation <- data_nmds_proteins_filtered %>%
   ggplot(aes(x = NMDS1,
              y = NMDS2,
              color = treatment,
              shape = incubation)) +
   geom_point(alpha = 0.8,
              size = 3) +
   geom_point(data = centroid, size = 3,
              shape = 21,
              color = "black",
              aes(fill=treatment)) +
   scale_color_manual(values = base_color,
                      guide=guide_legend(override.aes = list(color = base_color))) +
   scale_fill_manual(values = base_color) +
   labs(title = "Incubation") +
   theme_classic()+
   theme(legend.title = element_blank())
 
 prot_70p_substrate <- data_nmds_proteins_filtered %>% 
   ggplot(aes(x = NMDS1,
              y = NMDS2,
              color = treatment,
              shape = substrate)) +
   geom_point(alpha = 0.8,
              size = 3) +
   geom_point(data = centroid, size = 3,
              shape = 21,
              color = "black",
              aes(fill=treatment)) +
   scale_color_manual(values = base_color,
                      guide=guide_legend(override.aes = list(color = base_color))) +
   scale_fill_manual(values = base_color) +
   labs(title = "Substrate") +
   theme_classic()+
   theme(legend.title = element_blank(),
         legend.background = element_rect(color = NA))
 
# save plot --------------------------------------------------------------------
 prot_70p <- prot_70p_incubation + prot_70p_substrate

 ggsave(prot_70p, filename = "conservation_protocol/plots/ordination_prot70p.png",
        width = 10, height = 4, dpi = 450)

 dist <- vegdist(matrix_proteins_filtered,method = "bray")
 
 adonis2(as.dist(dist)~clean_meta$treatment*clean_meta$incubation*clean_meta$substrate, permutations = 1000)
 