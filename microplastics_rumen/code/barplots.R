# Load libraries ---------------------------------------------------------------
pacman::p_load(tidyverse,
               patchwork,
               showtext)

# Load data, set colors and font -----------------------------------------------
metadata <- read_tsv("microplastics_rumen/rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid))

taxonomy <- read_csv("microplastics_rumen/rawdata/BuiltIn.taxa.refine.csv") 

col_val <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
             '#fdb462','#b3de69','#fccde5','#d9d9d9', 'black', 'violet')

font_add_google("Open Sans", "opensans") #You need these for the font
showtext_auto()

# Wrangle data -----------------------------------------------------------------
df_tax <- taxonomy %>% 
  select(Phylum, contains("LFQ")) %>% 
  pivot_longer(-Phylum, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Phylum=if_else(is.na(Phylum), "Unclassified", Phylum)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "hay" | material == "Control") %>%
  filter(group != "RumenJ+Barley",
         group != "Barley",
         group != "uHDPE x uPVC") %>% 
  group_by(Phylum, sampleid, new_group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  ungroup() %>% 
  group_by(Phylum, new_group) %>% 
  summarise(mean_rel_abun = mean(rel_abun),
            sd = sd(rel_abun),
            .groups = "drop")

write.csv(df_tax, file = "microplastics_rumen/hay_phylum_summary.csv")

pool_tax <- df_tax %>% 
  group_by(Phylum) %>% 
  summarise(pool = max(mean_rel_abun) < 0.5,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax$new_group <- factor(df_tax$new_group,
                       levels = c("Hay+", "RF", "RF+Hay", "sPLA", "sPHB", 
                                  "sHDPE","sPVC", "lPLA", "lPHB", "lHDPE", 
                                  "lPVC", "lPP"))

hay_plot_phylum <- inner_join(pool_tax, df_tax) %>%
  mutate(Phylum = if_else(pool, "Other", Phylum)) %>% 
  group_by(new_group, Phylum) %>% 
  summarise(total_mean_rel = sum(mean_rel_abun),
            .groups = "drop"
            ) %>% 
  #filter(Phylum != "Unclassified") %>% 
  ggplot(aes(x = new_group,
             y = total_mean_rel,
             fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,100, 10)
  ) +
  scale_fill_manual(values = col_val) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        #legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Barley
df_tax_barley <- taxonomy %>% 
  select(Phylum, contains("LFQ")) %>% 
  pivot_longer(-Phylum, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Phylum=if_else(is.na(Phylum), "Unclassified", Phylum)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "barley" | material == "Control") %>%
  filter(group != "RumenJ+Hay",
         group != "Hay",
         group != "uHDPE x uPVC") %>% 
  group_by(Phylum, sampleid, new_group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  ungroup() %>% 
  group_by(Phylum, new_group) %>% 
  summarise(mean_rel_abun = mean(rel_abun),
            sd = sd(rel_abun),
            .groups = "drop")

write.csv(df_tax_barley, file = "microplastics_rumen/barley_phylum_summary.csv")

pool_tax <- df_tax_barley %>% 
  group_by(Phylum) %>% 
  summarise(pool = max(mean_rel_abun) < 0.2,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax_barley$new_group <- factor(df_tax_barley$new_group,
                       levels = c("Barley+", "RF", "RF+Barley", "sPLA", "sPHB", 
                                  "sHDPE","sPVC", "lPLA", "lPHB", "lHDPE", 
                                  "lPVC", "lPP"))

barley_plot_phylum <- inner_join(pool_tax, df_tax_barley) %>%
  mutate(Phylum = if_else(pool, "Other", Phylum)) %>% 
  group_by(new_group, Phylum) %>% 
  summarise(mean_rel = sum(mean_rel_abun),
            .groups = "drop") %>%
  #filter(Family != "Unclassified") %>%
  ggplot(aes(x = new_group,
             y = mean_rel,
             fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = col_val) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        #legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Plots at family level --------------------------------------------------------

df_tax_family_hay <- taxonomy %>% 
  select(Family, contains("LFQ")) %>% 
  pivot_longer(-Family, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Family=if_else(is.na(Family), "Unclassified", Family)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "hay" | material == "Control") %>%
  filter(group != "RumenJ+Barley",
         group != "Barley",
         group != "uHDPE x uPVC") %>% 
  group_by(Family, sampleid, new_group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  ungroup() %>% 
  group_by(Family, new_group) %>% 
  summarise(mean_rel_abun = mean(rel_abun),
            sd = sd(rel_abun),
            .groups = "drop")

write.csv(df_tax_family_hay, file = "microplastics_rumen/hay_family_summary.csv")

pool_tax <- df_tax_family_hay %>% 
  group_by(Family) %>% 
  summarise(pool = max(mean_rel_abun) < 0.6,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax_family_hay$new_group <- factor(df_tax_family_hay$new_group,
                           levels = c("Hay+", "RF", "RF+Hay", "sPLA", "sPHB", 
                                      "sHDPE","sPVC", "lPLA", "lPHB", "lHDPE", 
                                      "lPVC", "lPP"))

hay_plot_family <- inner_join(pool_tax, df_tax_family_hay) %>%
  mutate(Family = if_else(pool, "Other", Family)) %>% 
  group_by(new_group, Family) %>% 
  summarise(total_mean_rel = sum(mean_rel_abun),
            .groups = "drop") %>% 
  filter(Family != "Unclassified") %>% 
  ggplot(aes(x = new_group,
             y = total_mean_rel,
             fill = Family)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 30, 5)
  ) +
  scale_fill_manual(values = col_val) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        #legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Barley
df_tax_family_barley <- taxonomy %>% 
  select(Family, contains("LFQ")) %>% 
  pivot_longer(-Family, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Family=if_else(is.na(Family), "Unclassified", Family)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "barley" | material == "Control") %>%
  filter(group != "RumenJ+Hay",
         group != "Hay",
         group != "uHDPE x uPVC") %>% 
  group_by(Family, sampleid, new_group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  ungroup() %>% 
  group_by(Family, new_group) %>% 
  summarise(mean_rel_abun = mean(rel_abun),
            sd = sd(rel_abun),
            .groups = "drop")

write.csv(df_tax_family_barley, file = "microplastics_rumen/barley_family_summary.csv")

pool_tax <- df_tax_family_barley %>% 
  group_by(Family) %>% 
  summarise(pool = max(mean_rel_abun) < 0.5,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax_family_barley$new_group <- factor(df_tax_family_barley$new_group,
                                  levels = c("Barley+", "RF", "RF+Barley", "sPLA", "sPHB", 
                                             "sHDPE","sPVC", "lPLA", "lPHB", "lHDPE", 
                                             "lPVC", "lPP"))


col_val_barley <- c('#ffffb3','#bebada','#fb8072','#80b1d3',
             '#fdb462','#b3de69','#fccde5','#d9d9d9', 'black', 'brown', 'violet')


barley_plot_family <- inner_join(pool_tax, df_tax_family_barley) %>%
  mutate(Family = if_else(pool, "Other", Family)) %>% 
  group_by(new_group, Family) %>% 
  summarise(mean_rel = sum(mean_rel_abun),
            .groups = "drop") %>%
  filter(Family != "Unclassified") %>%
  ggplot(aes(x = new_group,
             y = mean_rel,
             fill = Family)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values = col_val_barley) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        #legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



hay <- hay_plot_phylum / hay_plot_family + 
  plot_annotation(tag_levels = 'A') &
  theme(#legend.position='bottom',
    text = element_text(size = 40),
    legend.text = element_text(size = 35))

ggsave(filename = "microplastics_rumen/plots/hay_tax_plot.png", width = 8, height = 6) 

barley <- barley_plot_phylum /   barley_plot_family +
  plot_annotation(tag_levels = 'A') &
  theme(#legend.position='bottom',
    text = element_text(size = 40),
    legend.text = element_text(size = 35))

ggsave(filename = "microplastics_rumen/plots/barley_tax_plot.png", width = 8, height = 6) 


