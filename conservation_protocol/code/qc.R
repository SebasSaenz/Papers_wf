# Johan S. Sáenz, Hohenheim University

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(here)
library(gt)

# Load dataset -----------------------------------------------------------------

quality <- read_tsv("conservation_protocol/rawdata/final_summary.tsv")

proteins <- read_tsv("conservation_protocol/rawdata/final_proteins.tsv")


# Make a summary table ---------------------------------------------------------

table <- quality %>% 
  select(-Experiment) %>% 
  filter(`Raw file` != "Total") %>% 
  pivot_longer(cols = -1, names_to = "Variable", values_to = "value") %>% 
  group_by(Variable) %>% 
  summarise(Max = round(max(value), 1),
            Min = round(min(value), 1),
            Mean = round(mean(value), 1),
            Sd = round(sd(value), 1))

table %>% 
  gt() %>% 
  tab_stubhead(label = "landmass") %>% 
  cols_label() %>% 
  gtsave("tab_1.docx")

# MS Percentage histogram
percentage <- quality %>%
  mutate(`Raw file` =str_remove(`Raw file`, "Seifert_Hammel_Nr")) %>% 
  ggplot(aes(x = `MS/MS Identified [%]`)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept= 22.7, linetype = "dashed", 
             color = "red", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 50, by = 2)) +
  scale_y_continuous(breaks = seq(0, 20, by = 2),
                     expand = c(0, 0)) +
  labs(y = "Number of samples") +
  theme_classic()

mean(quality$`MS/MS Identified [%]`)

# Peptides identified histogram
peptides <- quality %>%
  filter(`Raw file` != "Total") %>% 
  mutate(`Raw file` =str_remove(`Raw file`, "Seifert_Hammel_Nr")) %>% 
  ggplot(aes(x = `Peptide Sequences Identified`)) +
  geom_histogram() +
  geom_vline(xintercept = 5911, linetype = "dashed", 
             color = "red", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 12000, by = 1000)) +
  scale_y_continuous(breaks = seq(0, 20, by = 2),
                     expand = c(0, 0)) +
  labs(y = "Number of samples") +
  theme_classic()


mean(quality$`Peptide Sequences Identified`)


# Protein origin----------------------------------------------------------------
 filter_proteins <- proteins %>% 
  filter(grepl("SHEEP", `Majority protein IDs`))


number_proteins <- tibble(source = c("Microbiome", "Host"),
       number_proteins = c(24068, 586)) %>% 
  mutate(percentage = 100 * (number_proteins/sum(number_proteins))) %>% 
  ggplot(aes(x = source,
             y = percentage)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 100),
                     breaks = seq(0, 100, 10)) +
  labs(x = NULL,
       y = "Annotated proteins (%)") +
  geom_text(label ="24068", x = 2, y =90)+
  geom_text(label ="586", x = 1, y =5)+
  theme_classic()

ggsave(filename = "conservation_protocol/plots/percentage_proteins.png",
       width = 3, height = 3, dpi = 400)

# Compose plot -----------------------------------------------------------------
plot <- ((peptides / percentage) | number_proteins)  + 
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = 'A')

ggsave(plot, filename = "conservation_protocol/plots/qc_plot.png",
       width = 8, height = 4, dpi = 450)
