# Load libraries ----------------------------------------------------------------
library(tidyverse)

# Load files -------------------------------------------------------------------
df_quality <- read_tsv("Brazil_samples/final_summary.tsv")

df_origin <- read_tsv("Brazil_samples/final_proteins.tsv")

# Wrangle data -----------------------------------------------------------------

# Quality check

df_quality %>% 
  filter(`Raw file` != "Total") %>% 
  ggplot(aes(x = `MS/MS Identified [%]`)) +
  geom_histogram(binwidth = 1, aes(y=after_stat(density))) +
  geom_density(color = "red")+ 
  # scale_y_continuous(limits = c(0, 8),
  #                    breaks = seq(1, 8, 1),
  #                    expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  theme_bw() +
  theme(
        panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())


df_quality %>% 
  filter(`Raw file` != "Total") %>% 
  ggplot(aes(x = `Peptide Sequences Identified`)) +
  geom_histogram(binwidth = 100, aes(y=after_stat(density))) +
  geom_density(color = "red")+ 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 8500, 1000)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank())


# Origin of proteins
protein_origin <- df %>% 
  select(`Majority protein IDs`) %>% 
  separate(`Majority protein IDs`,
           into =c("protein"),
           sep = ";") %>% 
  mutate(origin = if_else(grepl("MGY", protein), "Microbial", "Host")) %>% 
  count(origin) 

protein_origin %>% 
  ggplot(aes(x=origin, y=n)) +
  geom_col()+
  scale_y_continuous(limits = c(0,1300),
                     breaks = seq(0, 1300, 200),
                     expand = c(0,0)) +
  labs(x = NULL,
       y = "Number of proteins") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
  )



df_origin %>% 
  select(`Majority protein IDs`, contains("Intensity")) %>% 
  separate(`Majority protein IDs`,
           into =c("protein"),
           sep = ";") %>% 
  mutate(origin = if_else(grepl("MGY", protein), "Microbial", "Host")) %>% 
  select(-c(Intensity, protein)) %>% 
  pivot_longer(-origin, names_to = "sample", values_to = "intensity") %>% 
  mutate(sample = str_remove(sample, "Intensity 240410_")) %>% 
  group_by(origin, sample) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(percent = sum_intensity/sum(sum_intensity)) %>% 
  ggplot(aes(x = sample,
             y = percent,
             fill = origin)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("grey", "red")) +
  labs(x = NULL,
       y = "Relative abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
  
