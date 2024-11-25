# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Load data --------------------------------------------------------------------
ph_df <- read_tsv("silage_temperature/rawdata/ph.txt")
metadata <- read_tsv("silage_temperature/rawdata/metadata.txt")

base_color <- c('#a6611a','#dfc27d','#80cdc1','#018571')

# Wrangle data -----------------------------------------------------------------

ph_df %>% 
  drop_na() %>% 
  inner_join(metadata,
            by = c("sampleid", "sample_number")) %>% 
  group_by(day, temperature) %>% 
  summarise(mean_ph = mean(ph),
            sd_ph = sd(ph),
            .groups = "drop") %>% 
  complete(day,
          temperature,
          fill=list(mean_ph=0)) %>% 
  mutate(day = factor(
         day,
         levels = c(0 ,2, 15, 45)),
         temperature = factor(temperature,
                              levels = c("Fresh", 25, 30, 35))) %>% 
  ggplot(aes(x = day,
             y = mean_ph,
             fill  = temperature)) + 
  geom_bar(stat="identity", 
           position=position_dodge2(), 
           width = 0.8) + 
        geom_errorbar(aes(ymin=mean_ph - sd_ph, ymax=mean_ph + sd_ph),
                        width=.2,
                        position=position_dodge(0.80)) +
  scale_y_continuous(limits = c(0, 7),
                     expand = c(0,0),
                     breaks = seq(0, 7, 1)) +
  scale_fill_manual(values = base_color) +
  labs(x = "Day",
       y = "pH") +
  theme_bw() +
  theme(
     legend.title = element_blank(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank())


