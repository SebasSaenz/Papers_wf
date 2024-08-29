

library(tidyverse)

annotation <- read_csv("microplastics_rumen/rawdata/proteins_annotations.csv") 

proteins <- read_tsv("microplastics_rumen/rawdata/protein_groups.txt")

metadata <- read_tsv("microplastics_rumen/rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid))

pazy <- read_tsv("rawdata/plastic_enzymes.txt") %>%
  pull(EC)

cluster1_barley <- read.csv("rawdata/Cluster_info_1.csv") %>% 
  select(Protein.IDs)
cluster6_barley <- read.csv("rawdata/Cluster_info_6.csv") %>% 
  select(Protein.IDs)
# --------------------
cazy <- annotation %>% 
  select(`Protein IDs`, CAZy) %>% 
  filter(CAZy != "-") %>% 
  separate(CAZy,
           into = "CAZy",
           sep = ",") %>% 
  mutate(classes = case_when(grepl("GH", CAZy) ~ "Glycoside Hydrolases (GHs)",
                             grepl("GT", CAZy) ~ "Glycosyl Transferases (GTs)",
         grepl("PL", CAZy) ~ "Polysaccharide Lyases (PLs)",
         grepl("CE", CAZy) ~ "Carbohydrate Esterases (CEs)",
         grepl("AA", CAZy) ~ "Auxiliary Activities (AAs)",
         grepl("CBM", CAZy) ~ "Carbohydrate-Binding Modules (CBMs)"))


df <-  proteins %>% 
  right_join(cazy, by = "Protein IDs") %>% 
  pivot_longer(cols = -c(`Protein IDs`, CAZy, classes), names_to = "sampleid", values_to = "intensity") %>%
  group_by(sampleid, classes) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  mutate(log_intensity = log2(sum_intensity + 1)) %>% 
  inner_join(metadata, by = "sampleid")

write.table(df, file = "rawdata/Cazy_families.txt", sep = "\t")

df %>%
  filter(matrix == "barley" | matrix == "control") %>%
  filter(group != "RumenJ+Barley",
         group != "uHDPE x uPVC") %>%
  ggplot(aes(x = factor(
    group,
    levels = c(
      "Barley", "RumenJ", "RumenJ+Barley", "nPLA",
      "nPHB", "nHDPE", "nPVC", "uPLA",
      "uPHB", "uHDPE", "uPVC", "uPP"
    )
  ),
  y = log_intensity)) +
  geom_boxplot(outlier.fill = NULL) +
  geom_jitter() +
  facet_wrap( ~ classes, scales = "free_y") + 
  labs(y = "Intensity (Log2)",
       x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8))

ggsave(filename = "figures/CAZy_families_barley.png", width = 8, height = 5, dpi = 450)

df %>%
  filter(matrix == "hay" | matrix == "control") %>%
  filter(group != "RumenJ+Hay",
         group != "uHDPE x uPVC") %>%
  ggplot(aes(x = factor(
    group,
    levels = c(
      "Hay", "HayJ", "RumenJ+Hay", "nPLA",
      "nPHB", "nHDPE", "nPVC", "uPLA",
      "uPHB", "uHDPE", "uPVC", "uPP"
    )
  ),
  y = log_intensity)) +
  geom_boxplot(outlier.fill = NULL) +
  geom_jitter(width = 0.2) +
  facet_wrap( ~ classes, scales = "free_y") + 
  labs(y = "Intensity (Log2)",
       x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8))

ggsave(filename = "figures/CAZy_families_hay.png", width = 8, height = 5, dpi = 450)

######################

x <- annotation %>% 
  select(`Protein IDs`, EC) %>% 
  filter(EC != "-") %>%
  separate(EC,
           into = "EC",
           sep = ",") %>% 
  filter(EC %in% pazy)

##############

cl_1_b <- annotation %>% 
  rename(Protein.IDs=`Protein IDs`) %>% 
  right_join(cluster1_barley, by = "Protein.IDs") %>% 
  select(COG_category) %>% 
  mutate(COG_category = case_when(COG_category == "CO"~"Other",
                                  COG_category == "CQ"~"Other",
                                  COG_category == "FG"~"Other",
                                  COG_category == "IQ"~"Other",
                                  COG_category == "KT"~"Other",
                                  COG_category == "NT"~"Other",
                                  COG_category == "NU"~"Other",
                                  COG_category == "OU"~"Other",
                                  COG_category == "-"~"Other",
                                  is.na(COG_category)~"Other",
                                  .default = as.character(COG_category))
                                  ) %>% 
  count(COG_category) %>% 
  ggplot(aes(x = fct_reorder(COG_category, n, .desc = TRUE),
             y = n)) +
  geom_col() +
  labs(x = "COG category",
       y = "Number of proteins") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 400)) +
  theme_classic() +
  geom_text(aes(label = n), vjust = -1, size = 3)

ggsave(cl_1_b, file = "figures/cl_1_b.png", width = 6, height = 4, dpi = 450) 




cl_6_b <- annotation %>% 
  rename(Protein.IDs=`Protein IDs`) %>% 
  right_join(cluster6_barley, by = "Protein.IDs") %>% 
  select(COG_category) %>% 
  mutate(COG_category =str_replace(COG_category, "..", "Other"),
         COG_category =str_replace(COG_category, "Other.", "Other"),
         COG_category = case_when(COG_category == "-"~"Other",
                                  is.na(COG_category)~"Other",
                                  .default = as.character(COG_category))
  ) %>% 
  count(COG_category) %>% 
  ggplot(aes(x = fct_reorder(COG_category, n, .desc = TRUE),
             y = n)) +
  geom_col() +
  labs(x = "COG category",
       y = "Number of proteins") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 400)) +
  theme_classic() +
  geom_text(aes(label = n), vjust = -1, size = 3)

ggsave(cl_6_b, file = "figures/cl_6_b.png", width = 6, height = 4, dpi = 450)
