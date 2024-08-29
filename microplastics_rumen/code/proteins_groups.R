pacman::p_load(tidyverse,
               vegan, 
               glue,
               here,
               patchwork,
               showtext)


metadata <- read_tsv("microplastics_rumen/rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid))

proteins <- read_tsv("microplastics_rumen/rawdata/protein_groups.txt") %>%
  rename(proteinid='Protein IDs')

hay_color <- c('#1b9e77','#66a61e','#e6ab02','#d95f02', '#a6761d','#666666')

barley_color <- c('#7570b3','#66a61e','#e6ab02','#e7298a','#a6761d','#666666')


font_add_google("Open Sans", "opensans") #You need these for the font
showtext_auto()




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

bray <- vegdist(nmds_df, method = "bray") 

pcoa <- cmdscale(bray, k=2, eig=TRUE, add=TRUE) 
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")

percent_explained <- 100 *pcoa$eig / sum(pcoa$eig)


positions %>% 
  as_tibble(rownames = "sampleid") %>%
  inner_join(metadata, by ="sampleid") %>% 
  ggplot() +
  geom_point(aes(x = pcoa1,
                 y = pcoa2,
                 color = size)) +
  labs(x="PCo 1 (28.3.7%)",
       y="PCo 2 (12.7%)")


metadata_Hay <- read_tsv("rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid)) %>%
  filter(matrix == "hay",
         group != "uHDPE x uPVC",
         size != "Hay")

adonis2(bray ~ size, 
        data = metadata_Hay, 
        permutations = 999)


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
hay_nmds <- data_nmds %>%
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
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))



######### Barley
nmds_df <- protein_pivot %>%
  filter(matrix == "barley" | matrix == "RumenJ" |matrix == "RumenJ+Barley",
         group != "uHDPE x uPVC") %>%
  pivot_wider(id_cols = c(1), names_from = "sampleid", values_from = "rel_abun")


nmds_df <- nmds_df[2:44] %>% 
  t()

bray <- vegdist(nmds_df, method = "bray") 

metadata_barley <- read_tsv("rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid)) %>%
  filter(matrix == "barley" | matrix == "RumenJ" |matrix == "RumenJ+Barley",
         group != "uHDPE x uPVC")

adonis2(bray ~ group, 
        data = metadata_barley, 
        permutations = 999)



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
barley_nmds <- data_nmds %>%
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
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


# Compose plot
compose_nmds <- (hay_nmds | barley_nmds) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom',
        text = element_text(size = 50),
        legend.text = element_text(size = 45))



ggsave("microplastics_rumen/plots/nmds_proteins.png", dpi = 450, height = 5, width = 11)







nmds_df <- protein_pivot %>%
  filter(matrix == "barley" & material == "HDPE") %>%
  pivot_wider(id_cols = c(1), names_from = "sampleid", values_from = "rel_abun")

nmds_df <- nmds_df[2:9] %>% 
  t()

bray <- vegdist(nmds_df, method = "bray")

metadata_barley <- read_tsv("rawdata/metadata.txt") %>% 
  mutate(sampleid = as.character(sampleid)) %>%
  filter(matrix == "barley" & material == "HDPE")

adonis_result <- adonis2(bray~size, 
        data = metadata_barley, 
        permutations = 999)

adonis_result$`Pr(>F)`

p.adjust(c(0.017, 0.026, 0.062, 0.032, 0.038, 0.028, 0.032, 0.03), 
         method = "BH")
