# set up working enviroment ----------------------------------------------------
pacman::p_load(tidyverse,
               here,
               patchwork,
               car)

clean_metadata <- readRDS("conservation_protocol/output/clean_meta.rds")

proteins <- read_tsv("conservation_protocol/rawdata/functions_report.tsv")

base_color <- c("#d8b365", "#67a9cf", "#2166ac", "#91cf60", "#1a9850")


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



cazy <- proteins %>% 
  filter(Protein_ID == 1,
               Name %in% protein_names,
               CAZy != "-") %>% 
  select(Name, CAZy) %>% 
  separate(CAZy,
           into = "CAZy",
           sep = ",") %>% 
  mutate(classes = case_when(grepl("GH", CAZy) ~ "Glycoside Hydrolases (GHs)",
                             grepl("GT", CAZy) ~ "Glycosyl Transferases (GTs)",
                             grepl("PL", CAZy) ~ "Polysaccharide Lyases (PLs)",
                             grepl("CE", CAZy) ~ "Carbohydrate Esterases (CEs)",
                             grepl("AA", CAZy) ~ "Auxiliary Activities (AAs)",
                             grepl("CBM", CAZy) ~ "Carbohydrate-Binding Modules (CBMs)"))


cazy_df <- proteins %>%  
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
   separate(CAZy,
            into = "CAZy",
            sep = ",") %>% 
   mutate(classes = case_when(grepl("GH", CAZy) ~ "Glycoside Hydrolases (GHs)",
                              grepl("GT", CAZy) ~ "Glycosyl Transferases (GTs)",
                              grepl("PL", CAZy) ~ "Polysaccharide Lyases (PLs)",
                              grepl("CE", CAZy) ~ "Carbohydrate Esterases (CEs)",
                              grepl("AA", CAZy) ~ "Auxiliary Activities (AAs)",
                              grepl("CBM", CAZy) ~ "Carbohydrate-Binding Modules (CBMs)")) %>% 
  group_by(sample, classes) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  mutate(log_intensity = log2(sum_intensity + 1)) %>% 
  inner_join(clean_metadata, by = "sample") %>% 
   filter(treatment != "F")
 
 
# Make Cazy tables

write_tsv(cazy_df, file = "conservation_protocol/output/Cazy_abundance.txt")
  
pl_df <- cazy_df %>% 
   filter(classes == "Polysaccharide Lyases (PLs)")
  
 
res.aov <- aov(log_intensity ~ treatment, data = cbm_df)
summary(res.aov)
 

leveneTest(log_intensity ~ treatment, data = gh_df)
  
TukeyHSD(res.aov) 


cbm_plot <- cbm_df %>% 
ggplot(aes(x = treatment,
             y = log_intensity)) +
  geom_boxplot(outlier.colour = "white",
               width = 0.4) +
  geom_jitter(aes(color = treatment),
              width = 0.2,
              alpha = 0.7) +
  # geom_line(data=tibble(x=c(1, 4), 
  #                       y=c(30.7, 30.7)), 
  #           aes(x=x, y=y), inherit.aes=FALSE, linewidth = 0.3) +
  # geom_text(data=tibble(x=2.5, 
  #                       y=30.8), 
  #           aes(x=x, y=y, label="**"), inherit.aes=FALSE) +
  # geom_line(data=tibble(x=c(1, 5), 
  #                       y=c(31.5, 31.5)), 
  #           aes(x=x, y=y), inherit.aes=FALSE, linewidth = 0.3) +
  # geom_text(data=tibble(x=3, 
  #                       y=31.6), 
  #           aes(x=x, y=y, label="**"), inherit.aes=FALSE) +
  facet_wrap(~classes, scales = "free_y") +
  labs(x = NULL,
       y = "Abundance (Log2 Intensity)") +
   scale_y_continuous(expan = c(0,0),
                      limits = c(22, 32),
                      breaks = seq(22, 32, 1)) +
  scale_color_manual(values = base_color) +
  theme_minimal() +
  theme(panel.background = element_blank(),
        axis.ticks = element_line(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.box.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


compose_plot <- ce_plot + gh_plot + gt_plot + cbm_plot + plot_layout(guides = "collect",
                                                     axes = "collect")

ggsave(compose_plot, file = "conservation_protocol/plots/cazy_plot.png",
       width = 8,
       height = 6,
       dpi = 400)
