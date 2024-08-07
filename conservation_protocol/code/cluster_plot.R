library(tidyverse)

cluster_1 <- read_csv(file = "conservation_protocol/rawdata/Cluster_info_1.csv")
cluster_2 <- read_csv(file = "conservation_protocol/rawdata/Cluster_info_2.csv")
cluster_3 <- read_csv(file = "conservation_protocol/rawdata/Cluster_info_3.csv")

functions  <- read_tsv(file = "conservation_protocol/rawdata/functions_report.tsv") 

# Clean df ----------------------

clean_cluster1 <- cluster_1 |>
  select(`Protein IDs`) |>
  rename(Name = `Protein IDs`)

clean_cluster2 <- cluster_2 |>
  select(`Protein IDs`) |>
  rename(Name = `Protein IDs`)

clean_cluster3 <- cluster_3 |>
  select(`Protein IDs`) |>
  rename(Name = `Protein IDs`)


df_c1 <- functions %>% 
  select(Name, `COG category`) %>%
  mutate(Cluster = "Cluster 1") %>% 
  inner_join(clean_cluster1, by = "Name")

df_c2 <- functions %>% 
  select(Name, `COG category`) %>%
  mutate(Cluster = "Cluster 2") %>% 
  inner_join(clean_cluster2, by = "Name")

df_c3 <- functions %>% 
  select(Name, `COG category`) %>%  
  mutate(Cluster = "Cluster 3") %>% 
  inner_join(clean_cluster3, by = "Name")

rbind(df_c1, df_c2, df_c3) %>%
  rename(cog = `COG category`) %>%
  mutate(cog = case_when(cog == "CHR"~ "Other",
                         cog == "CR"~ "Other",
                         cog == "EF"~ "Other",
                         cog == "EG"~ "Other",
                         cog == "EH"~ "Other",
                         cog == "EM"~ "Other",
                         cog == "EP"~ "Other",
                         cog == "ER"~ "Other",
                         cog == "ET"~ "Other",
                         cog == "FH"~ "Other",
                         cog == "FV"~ "Other",
                         cog == "GK"~ "Other",
                         cog == "GT"~ "Other",
                         cog == "HC"~ "Other",
                         cog == "HE"~ "Other",
                         cog == "HP"~ "Other",
                         cog == "IQ"~ "Other",
                         cog == "JR"~ "Other",
                         cog == "KT"~ "Other",
                         cog == "MO"~ "Other",
                         cog == "PV"~ "Other",
                         cog == "CHR"~ "Other",
                         cog == "TG"~ "Other",
                         cog == "TK"~ "Other",
                         cog == "UW"~ "Other",
                         is.na(cog)~ "Other",
                         .default = cog)) %>%
  count(cog, Cluster) %>% 
  # mutate(abundance = 100*(n/sum(n))) %>% 
  complete(cog, Cluster, fill = list(n = 0)) %>% 
  ggplot(aes(x = cog,
             y = n,
             fill = Cluster)) +
  geom_col(position = "dodge") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 15),
                     breaks = seq(0, 15, 1)) +
  scale_fill_manual(values = c("#e66101", "#fdb863", "#5e3c99")) +
  labs(x = "COG Categories",
       y = "Protein groups") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.6))

ggsave(filename = "conservation_protocol/plots/barplot_cog.png",
       width = 8, height = 4,
       dpi = 400)

test <- functions %>% 
  inner_join(clean_clsuter3, by = "Name") %>% 
  filter(`COG category` == "J")
