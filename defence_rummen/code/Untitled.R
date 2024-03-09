#set working directory
setwd("/Users/sebastiansaenz/Documents/defence_rummen")


library(tidyverse)

taxonomy <- read_tsv("rawdata/gtdbtk_all_summary.tsv")

 x <- taxonomy %>% 
  select(user_genome, classification) %>% 
  separate(classification,
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>% 
  select("genus") %>% 
  mutate(genus = str_remove(genus, ".__")) %>%
   filter(!is.na(genus)) %>% 
  unique() %>% 
   count()


 
 tibble(tax_level = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
        count = c(2, 27, 34, 73, 136, 578))
 