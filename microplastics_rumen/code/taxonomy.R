# Load libraries ---------------------------------------------------------------
pacman::p_load(tidyverse,
               vegan,
               here,
               patchwork,
               showtext)

# Load data and set colors  ----------------------------------------------------
metadata <- read_tsv("microplastics_rumen/rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid))

taxonomy <- read_csv("microplastics_rumen/rawdata/BuiltIn.taxa.refine.csv") 

hay_color <- c('#1b9e77','#66a61e','#e6ab02','#d95f02', '#a6761d','#666666')

barley_color <- c('#7570b3','#66a61e','#e6ab02','#e7298a','#a6761d','#666666')

col_val <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
             '#fdb462','#b3de69','#fccde5','#d9d9d9', 'black', 'violet')

font_add_google("Open Sans", "opensans") #You need these for the font
showtext_auto()

# NMDS using taxonomy ----------------------------------------------------------

# Hay
taxonomy$row_num <- seq.int(nrow(taxonomy)) #make each row unique

tax_pivot <- taxonomy %>%  # make a long data frame for NMDS
  pivot_longer(cols = c(11:94), 
               names_to = "sampleid", 
               values_to = "lfq") %>%
  mutate(sampleid=str_replace(sampleid,"LFQ intensity ", "")) %>% 
  select(row_num, everything()) %>%
  inner_join(metadata, by = "sampleid")

nmds_df <- tax_pivot %>%
  filter(matrix == "hay" | matrix == "RumenJ" |matrix == "RumenJ+Hay",
         group != "uHDPE x uPVC") %>%
  pivot_wider(id_cols = c(1:11), names_from = "sampleid", values_from = "lfq")

nmds_df <- nmds_df[c(12:54)] %>% #calculate relative abundance and create a matrix 
  apply(.,2, function(x){(x/sum(x))*100}) %>%
  t() 

nmds1 <- metaMDS(nmds_df,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)


nmds_best <- metaMDS(nmds_df, distance = "bray", #find best nmds
                     try = 20, 
                     trymax = 100, 
                     maxit = 1000, 
                     k = 3,
                     previous.best = nmds1)

nmds_best$stress  
stressplot(nmds_best)
plot(nmds_best)


#plotting with ggplot
data.scores <- as.data.frame(scores(nmds_best, display=c("sites")))

#Addd metadata to dataframe
data.scores$sampleid <- as.character(row.names(data.scores)) 

#joing metadata nmds scores
data_nmds <- inner_join(data.scores, metadata, by = "sampleid")

#make plot hay
hay_nmds_tax <- data_nmds %>%
  ggplot() +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =new_size,
                 shape= material),
             size = 3,
             alpha = 0.8) +
  scale_color_manual(values = hay_color) +
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8)) +
  #scale_color_manual(values = c('#66C2A5','#FC8D62')) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) #remove legend title

# Barley
nmds_df <- tax_pivot %>%
  filter(matrix == "barley" | matrix == "RumenJ" |matrix == "RumenJ+Barley",
         group != "uHDPE x uPVC") %>%
  pivot_wider(id_cols = c(1:11), names_from = "sampleid", values_from = "lfq")

nmds_df <- nmds_df[c(12:54)] %>% #calculate relative abundance and create a matrix 
  apply(.,2, function(x){(x/sum(x))*100}) %>%
  t() 

nmds1 <- metaMDS(nmds_df,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)


nmds_best <- metaMDS(nmds_df, distance = "bray", #find best nmds
                     try = 20, 
                     trymax = 100, 
                     maxit = 1000, 
                     k = 3,
                     previous.best = nmds1)

nmds_best$stress  
stressplot(nmds_best)
plot(nmds_best)


#plotting with ggplot
data.scores <- as.data.frame(scores(nmds_best, display=c("sites")))

#Addd metadata to dataframe
data.scores$sampleid <- as.character(row.names(data.scores)) 

#joing metadata nmds scores
data_nmds <- inner_join(data.scores, metadata, by = "sampleid")


#make plot barley
barley_nmds_tax <- data_nmds %>%
  ggplot() +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =new_size,
                 shape= material),
             size = 3,
             alpha = 0.8) +
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8)) +
  scale_color_manual(values = barley_color) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Compose plot

compose_nmds <- (hay_nmds_tax | barley_nmds_tax) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom',
        text = element_text(size = 50),
        legend.text = element_text(size = 45))



ggsave("microplastics_rumen/plots/nmds.png", dpi = 450, height = 5, width = 11)


# Barplot
df_tax <- taxonomy %>% 
  select(Family, contains("LFQ")) %>% 
  pivot_longer(-Family, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Family=if_else(is.na(Family), "Unclassified", Family)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "hay" | matrix == "control") %>%
  filter(group != "RumenJ+Barley",
         group != "uHDPE x uPVC") %>% 
  group_by(Family, sampleid, group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  group_by(Family, group) %>% 
  summarise(mean_rel_abun = mean(rel_abun), .groups = "drop")

pool_tax <- df_tax %>% 
  group_by(Family) %>% 
  summarise(pool = max(mean_rel_abun) < 0.5,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax$group <- factor(df_tax$group,
                       levels = c("Hay", "RumenJ", "RumenJ+Hay", "nPLA", "nPHB", 
                                  "nHDPE","nPVC", "uPLA", "uPHB", "uHDPE", 
                                  "uPVC", "uPP"))

inner_join(pool_tax, df_tax) %>%
  mutate(Family = if_else(pool, "Other", Family)) %>% 
  group_by(group, Family) %>% 
  summarise(total_mean_rel = sum(mean_rel_abun),
            .groups = "drop") %>% 
  filter(Family != "Unclassified") %>% 
  ggplot(aes(x = group,
             y = total_mean_rel,
             fill = Family)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,30, 5)
                     ) +
  scale_fill_manual(values = col_val) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("figures/hay_family.png", width = 6, height = 4)

############# Barley ###########



########## Barplots

df_tax <- taxonomy %>% 
  select(Family, contains("LFQ")) %>% 
  pivot_longer(-Family, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Family=if_else(is.na(Family), "Unclassified", Family)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "barley" | matrix == "control") %>%
  filter(group != "RumenJ+Hay",
         group != "uHDPE x uPVC") %>% 
  group_by(Family, sampleid, group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>% 
  group_by(Family, group) %>% 
  summarise(mean_rel_abun = mean(rel_abun), .groups = "drop")


pool_tax <- df_tax %>% 
  group_by(Family) %>% 
  summarise(pool = max(mean_rel_abun) < 1,
            mean = mean(mean_rel_abun),
            .groups = "drop")

df_tax$group <- factor(df_tax$group,
                       levels = c("Barley", "RumenJ", "RumenJ+Barley", "nPLA", "nPHB", 
                                  "nHDPE","nPVC", "uPLA", "uPHB", "uHDPE", 
                                  "uPVC", "uPP"))

inner_join(pool_tax, df_tax) %>%
  mutate(Family = if_else(pool, "Other", Family)) %>% 
  group_by(group, Family) %>% 
  summarise(mean_rel = sum(mean_rel_abun),
            .groups = "drop") %>%
  filter(Family != "Unclassified") %>%
  ggplot(aes(x = group,
             y = mean_rel,
             fill = Family)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,30,5)) +
  scale_fill_manual(values = col_val) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("figures/barley_family.png", width = 6, height = 4)

taxonomy %>% 
  select(Phylum, contains("LFQ")) %>% 
  pivot_longer(-Phylum, 
               names_to = "sampleid",
               values_to = "LFQ") %>% 
  mutate(sampleid = str_remove(sampleid, "LFQ intensity "),
         Phylum=if_else(is.na(Phylum), "Unclassified", Phylum)) %>% 
  inner_join(metadata, by = "sampleid") %>% 
  filter(matrix == "barley" | matrix == "control") %>%
  filter(group != "RumenJ+Hay",
         group != "uHDPE x uPVC") %>% 
  group_by(Phylum, sampleid, group) %>% 
  summarise(sum_lfq = sum(LFQ), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abun = 100*(sum_lfq/sum(sum_lfq))) %>%
  filter(Phylum=="Proteobacteria") %>% 
  ggplot(aes(x = factor(group,
                        levels = c("Barley", "RumenJ", "RumenJ+Barley", "nPLA", "nPHB", 
                                   "nHDPE","nPVC", "uPLA", "uPHB", "uHDPE", 
                                   "uPVC", "uPP")),
             y = rel_abun)) +
  geom_boxplot() +
  labs(x=NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



(hay_nmds_tax | barley_nmds_tax) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')
