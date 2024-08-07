pacman::p_load(tidyverse,
               viridis,
               Hmisc,
               ggpubr)

setwd("~/Documents/Lab_stuff/microplastics")

metabolites <- read_tsv("rawdata/metabolites.txt")

lfq_barley <- read_csv("rawdata/Results_lfq_barley.csv") %>% 
  rename(protein_id="Protein IDs")

lfq_hay <- read_csv("rawdata/Results_lfq_hay.csv") %>% 
  rename(protein_id="Protein IDs")

proteins_annotation <- read_csv("rawdata/proteins_annotations.csv") %>% 
  rename(protein_id=X.query)

protein_groups <- read_tsv("rawdata/protein_groups.txt") %>%
  rename(protein_id="Protein IDs")

metadata <- read_tsv("rawdata/metadata.txt")
metadata$sampleid <- as.character(metadata$sampleid)

colors_cog <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99', '#8dd3c7',
                '#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69',
                '#fccde5','#d9d9d9','#bc80bd','#ccebc5')

lfq_barley %>% 
  filter(significant == "TRUE") %>% 
  select(protein_id) %>% 
  inner_join(protein_groups, by="protein_id") %>%
  pivot_longer(-protein_id,
               names_to = "sampleid",
               values_to = "lfq") %>% 
  inner_join(proteins_annotation, by="protein_id") %>% 
  select(protein_id, sampleid, lfq, COG_category) %>% 
  inner_join(metadata, by="sampleid") %>% 
  filter(matrix=="barley" | matrix =="RumenJ",
         #group != "Barley",
         group != "uHDPE x uPVC") %>%
  select(sampleid, group, COG_category, lfq) %>% 
  mutate(COG_category=str_replace(COG_category, "..", "Other"),
         COG_category=str_replace(COG_category, "Other.", "Other"),
         COG_category=str_replace(COG_category, "-", "Other"),
         COG_category=if_else(is.na(COG_category), "Other", COG_category),
         group=factor(group, 
                      levels = c("Barley", "RumenJ", "nHDPE","nPHB",
                                                 "nPLA", "nPVC", "uHDPE", "uPHB", "uPLA",
                                                 "uPP", "uPVC")
                                                 )) %>%
  group_by(sampleid, group, COG_category) %>%
  summarise(sum_lfq = sum(lfq), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abundance = 100*(sum_lfq/sum(sum_lfq))) %>% 
  group_by(group, COG_category) %>%
  summarise(mean_rel_abundance =mean(rel_abundance), .groups = "drop") %>% 
  ggplot(aes(x=group,
             y=mean_rel_abundance,
             fill = COG_category)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors_cog) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
ggsave(filename = "figures/cog_barplot_sig_barley.png", width = 7, height = 5)
  
  

lfq_hay %>% 
  filter(significant == "TRUE") %>% 
  select(protein_id) %>% 
  inner_join(protein_groups, by="protein_id") %>%
  pivot_longer(-protein_id,
               names_to = "sampleid",
               values_to = "lfq") %>% 
  inner_join(proteins_annotation, by="protein_id") %>% 
  select(protein_id, sampleid, lfq, COG_category) %>% 
  inner_join(metadata, by="sampleid") %>% 
  filter(matrix=="hay" | matrix == "RumenJ",
         #group != "Hay",
         group != "uHDPE x uPVC") %>%
  select(sampleid, group, COG_category, lfq) %>% 
  mutate(COG_category=str_replace(COG_category, "..", "Other"),
         COG_category=str_replace(COG_category, "Other.", "Other"),
         COG_category=str_replace(COG_category, "-", "Other"),
         COG_category=if_else(is.na(COG_category), "Other", COG_category),
         group=factor(group, 
                      levels = c("Hay", "RumenJ", "nHDPE","nPHB",
                                 "nPLA", "nPVC", "uHDPE", "uPHB", "uPLA",
                                 "uPP", "uPVC")
         )) %>% 
  group_by(sampleid, group, COG_category) %>%
  summarise(sum_lfq = sum(lfq), .groups = "drop") %>%
  group_by(sampleid) %>% 
  mutate(rel_abundance = 100*(sum_lfq/sum(sum_lfq))) %>% 
  group_by(group, COG_category) %>%
  summarise(mean_rel_abundance =mean(rel_abundance), .groups = "drop") %>%
  filter(#group != "uPLA",
         #group != "uHDPE",
         #group != "uPHB"
    ) %>% 
  ggplot(aes(x=group,
             y=mean_rel_abundance,
             fill = COG_category)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors_cog) +
  labs(y = "Relative abundance (%)",
       x = NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(filename = "figures/cog_barplot_sig_hay.png", width = 7, height = 5)


lfq_barley %>% 
  filter(significant == "TRUE") %>% 
  select(protein_id) %>% 
  inner_join(protein_groups, by="protein_id") %>%
  pivot_longer(-protein_id,
               names_to = "sampleid",
               values_to = "lfq") %>% 
  inner_join(proteins_annotation, by="protein_id") %>% 
  select(protein_id, sampleid, lfq, CAZy) %>% 
  filter(CAZy != "-") %>% 
  group_by(sampleid, CAZy) %>% 
  summarise(sum_lfq = sum(lfq), .groups = "drop") %>% 
  inner_join(metadata, by = "sampleid") %>%
  filter(matrix=="barley",
         group != "Barley",
         group != "uHDPE x uPVC") %>% 
  group_by(group, CAZy) %>% 
  summarise(mean_lfq = mean(sum_lfq), .groups = "drop") %>% 
  ggplot(aes(x=group,
             y=mean_lfq)) +
  geom_col(fill = "blue") + 
  facet_wrap(~CAZy, nrow = 4) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "LFQ",
       x = NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(filename = "figures/CAZy_barplots.png", width = 5, height = 5)  


########### match proteins  

barely_results <- proteins_annotation %>% 
                  inner_join(lfq_hay, by = "protein_id") %>% 
  filter(significant == "TRUE")

write_tsv(barely_results, file = "hay_results.txt") 

###### correlation #####

 wide_barley <- lfq_barley %>% 
  filter(significant == "TRUE") %>% 
  select(protein_id) %>% 
  inner_join(protein_groups, by="protein_id") %>%
  pivot_longer(-protein_id,
               names_to = "sampleid",
               values_to = "lfq") %>% 
  inner_join(proteins_annotation, by="protein_id") %>% 
  select(protein_id, sampleid, lfq) %>% 
  inner_join(metadata, by="sampleid") %>% 
  filter(matrix=="barley",
         group != "Barley",
         group != "uHDPE x uPVC",
         sampleid != "21") %>%
  select(sampleid, protein_id, lfq) %>% 
  pivot_wider(sampleid, 
              names_from = protein_id,
              values_from =  lfq)

wide_barley$sampleid <- as.character(wide_barley$sampleid)
metabolites$sampleid <- as.character(metabolites$sampleid)



full_dataset <- inner_join(wide_barley, metabolites,
           by = "sampleid")

res2 <- rcorr(as.matrix(full_dataset[2:2468]), type = "pearson")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

list_metabolites <- metabolites %>% 
  pivot_longer(-sampleid,
               names_to = "meta",
               values_to = "values") %>% 
  select(meta) %>% 
  unique()

correlation <- flattenCorrMatrix(res2$r, res2$P) %>%
  mutate(p_adjusted=p.adjust(p,method = "fdr")) %>% 
  filter(p_adjusted < 0.001) %>% 
  filter(!grepl("MGY", row) | !grepl("MGY", column)) %>% 
  filter(!grepl("\\.", row) | !grepl("\\.", column)) %>%
  rename(protein_id=row) %>% 
  inner_join(proteins_annotation, by="protein_id") 


write.table(correlation, 
            file = "rawdata/pearson_result.txt", 
            sep = "\t",
            row.names = FALSE)

protein_groups %>% 
  filter(protein_id == "MGYG000290862_01036") %>% 
  pivot_longer(-protein_id,
               names_to = "sampleid",
               values_to = "lfq") %>% 
  inner_join(metadata, by= "sampleid") %>% 
  inner_join(metabolites, by= "sampleid") %>% 
  filter(matrix=="barley",
         group != "Barley",
         group != "uHDPE x uPVC",
         sampleid != "21"
         ) %>% 
  rename(b7_46='7.46') %>% 
  ggscatter(x = "lfq", 
            y = "Dimethylamine",
            alpha = 0,
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson") +
  geom_point(aes(color = size,
                 shape = material),
             size=3) +
  scale_shape_manual(values = c(15, 17 ,18,19, 9)) +
  scale_color_manual(values = c("red", "blue")) +
  labs(x="Pyruvate formate lyase-like \n(MGYG000290862_01036)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin())



ggsave(filename = "figures/correlation_4.png", width = 6, height = 6)
