# setting ----------------------------------------------------------------------
setwd("/Users/sebastiansaenz/Documents/defence_rummen")

library(tidyverse)
library(data.table)
library(pheatmap)
library(ggridges)
library(ggtext)
library(ggdist)
library(gghalves)
library(patchwork)
library(viridis)

colors <- c('#1f78b4', '#1b9e77')

# concatenate files ------------------------------------------------------------
all_paths <- list.files(path = "rawdata/padloc/csv/",
                        pattern = "*.csv",#List paths
                        full.names = TRUE)

all_content <- all_paths %>% #Join files
  lapply(read.table,
         header = TRUE,
         sep =",",
         encoding ="UTF-8")

all_filenames <- all_paths %>% #Manipulate file paths
  basename() %>%
  as.list()

all_lists <- mapply(c, 
                    all_content, 
                    all_filenames, 
                    SIMPLIFY = FALSE)

all_result <- rbindlist(all_lists, fill = T)
names(all_result)[20] <- "genome" #change column name

# xxx -------------------------------------------------------------------------
quality_genomes <- read_tsv("rawdata/quality_report_genomes.tsv") %>% 
  rename(genome=Name)

taxonomy <- read_tsv("rawdata/gtdbtk_all_summary.tsv") %>% 
  rename(genome = user_genome)

hq_genomes <- read_tsv("rawdata/hq_genomes.txt")

taxonomy <- taxonomy %>%
  select(genome, classification) %>%
  rename_all(tolower) %>%
  mutate(classification=str_replace_all(classification, "\\w+__", "")) %>%
  separate(classification,
           into=c("kingdom", "phylum", "class", "order", "family", "genus",
                  "species"),
           sep=";") %>% 
  inner_join(hq_genomes)


# number-system ----------------------------------------------------------------
number_system <- all_result %>% #all genomes with DS
  select(genome, system.number, system) %>%
  filter(!grepl('_other', system)) %>% 
  unique() %>% 
  mutate(genome = str_remove(genome, ".fasta_padloc.csv")) %>% 
  count(genome) %>% 
  right_join(taxonomy, by = "genome") %>% 
  mutate(n = replace_na(n, 0)) %>% 
  select(genome, kingdom, n) 

# proviruses -------------------------------------------------------------------

viral_qc <- read_tsv("rawdata/viral_quality_summary.tsv")

provirus <- viral_qc %>% 
  filter(provirus == "Yes",
         viral_genes >=1,
         proviral_length > 10000) %>% 
  select(contig_id) %>% 
  mutate(contig_id=str_replace_all(contig_id, "_sc.*", ""),
         contig_id=str_replace_all(contig_id, "_\\d+.*", ""),
         contig_id=str_replace_all(contig_id, "_deg.*", "")) %>% 
  group_by(contig_id) %>% 
  summarise(proviruses=n()) %>% 
  rename(genome=contig_id)

provirus_tree <- provirus %>%
  mutate(proviruses = if_else(proviruses >0, 1,0))
  
write_tsv(provirus_tree, file = "rawdata/provirus_counts.txt")

ds_provirus <- left_join(taxonomy, provirus, by = "genome") %>%
  select(genome, proviruses) %>% 
  left_join(number_system, provirus, by = "genome") %>% 
  mutate(proviruses=if_else(is.na(proviruses), 0, proviruses),
         proviruses = factor(proviruses, levels = c(0, 1, 2, 3, 4, 5))) %>% 
  ggplot(aes(x = proviruses,
             y = n,
             fill = kingdom)) +
  geom_point(shape = 21,
             position = position_jitter(width = 0.2,
                                        seed = 0),
             size = 3,
             alpha = 0.8,
             color = "black") +
  stat_summary(fun=median, 
               show.legend=FALSE, 
               geom = "crossbar",
               color="black",
               linewidth = 0.5) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 35),
                     breaks = seq(0, 35, 5),
                     oob = scales::squish) +
  labs(y = "Defence system per genome",
       x = "Number of proviruses") +
  facet_grid(cols = vars(kingdom)) +
  theme_classic() +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        panel.spacing = unit(2, "lines")) 

ds_provirus  


