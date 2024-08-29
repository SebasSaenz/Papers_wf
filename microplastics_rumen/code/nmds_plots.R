# NMDS PLOT
#  By Johan S. Sáenz

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

proteins <- read_tsv("microplastics_rumen/rawdata/protein_groups.txt") %>%
  rename(proteinid='Protein IDs')

hay_color <- c('#1b9e77','#66a61e','#e6ab02','#d95f02', '#a6761d','#666666')

barley_color <- c('#7570b3','#66a61e','#e6ab02','#e7298a','#a6761d','#666666')

font_add_google("Open Sans", "opensans") # You need these for the font
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
  stat_ellipse(linetype = 2,
                         aes(x = NMDS1, 
                             y = NMDS2,
                             colour = new_size),
                         level = 0.95) +
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
  stat_ellipse(linetype = 2,
               aes(x = NMDS1, 
                   y = NMDS2,
                   colour = new_size),
               level = 0.95) +
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

## NMDS using proteins ##

protein_pivot <- proteins %>%  # make a long data frame for NMDS
  pivot_longer(cols = c(2:85), 
               names_to = "sampleid", 
               values_to = "lfq") %>%
  inner_join(metadata, by = "sampleid") %>% 
  group_by(sampleid) %>% 
  mutate(rel_abun = lfq/sum(lfq)) %>% 
  ungroup()


nmds_df <- protein_pivot %>%
  filter(matrix == "hay" | matrix == "RumenJ" |matrix == "RumenJ+Hay",
         group != "uHDPE x uPVC") %>%
  pivot_wider(id_cols = c(1), names_from = "sampleid", values_from = "rel_abun")

nmds_df <- nmds_df[2:44] %>% 
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

# fast ploting and qc checking
nmds_best$stress  
stressplot(nmds_best)
plot(nmds_best)

#plotting with ggplot
data.scores <- as.data.frame(scores(nmds_best, display=c("sites")))

#Addd metadata to dataframe
data.scores$sampleid <- as.character(row.names(data.scores)) 


#joing metadata nmds scores
data_nmds <- left_join(data.scores, metadata, by = "sampleid")


#make plot 
hay_nmds_proteins <- data_nmds %>%
  ggplot() +
  stat_ellipse(linetype = 2,
               aes(x = NMDS1, 
                   y = NMDS2,
                   colour = new_size),
               level = 0.95) +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =new_size,
                 shape= material),
             size = 3) +
  scale_color_manual(values = hay_color) +
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8)) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "none") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Barley
nmds_df <- protein_pivot %>%
  filter(matrix == "barley" | matrix == "RumenJ" |matrix == "RumenJ+Barley",
         group != "uHDPE x uPVC") %>%
  pivot_wider(id_cols = c(1), names_from = "sampleid", values_from = "rel_abun")


nmds_df <- nmds_df[2:44] %>% 
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

# fast ploting and qc checking
nmds_best$stress  
stressplot(nmds_best)
plot(nmds_best)




#plotting with ggplot
data.scores <- as.data.frame(scores(nmds_best, display=c("sites")))

#Addd metadata to dataframe
data.scores$sampleid <- as.character(row.names(data.scores)) 


#joing metadata nmds scores
data_nmds <- left_join(data.scores, metadata, by = "sampleid")


#make plot 
barley_nmds_proteins <- data_nmds %>%
  ggplot() +
  stat_ellipse(linetype = 2,
               aes(x = NMDS1, 
                   y = NMDS2,
                   colour = new_size),
               level = 0.95) +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =new_size,
                 shape= material),
             size = 3) +
  scale_color_manual(values = barley_color) +
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8)) +
  theme_bw() +
  theme(text = element_text(family = "opensans"),
        panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = "none") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


## Compose plot ##

compose_nmds <- (hay_nmds_tax | hay_nmds_proteins) / (barley_nmds_tax | barley_nmds_proteins) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(#legend.position='bottom',
        text = element_text(size = 60),
        legend.text = element_text(size = 60))


ggsave("microplastics_rumen/plots/nmds_compose.png", dpi = 450, height = 8, width = 12)
