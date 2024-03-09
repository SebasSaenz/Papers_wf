metadata <- read_tsv("conservation_protocol/rawdata/metadata.txt")

clean_meta <- metadata %>%
  separate(treatment,
           into = c("sampling", "treatment", "incubation", "substrate"),
           sep = "-") %>% 
  mutate(incubation = str_replace(incubation, "1", "Pre-treated"),
         incubation = str_replace(incubation, "0[a-b]", "No incubation"),
         incubation = factor(incubation,
                             levels = c("No incubation", "Pre-treated", "8", "24", "72")),
         substrate = if_else(is.na(substrate), "None", substrate),
         treatment = factor(treatment,
                            levels = c("C", "W", "F", "FN", "FD", "FD_noN", "FD_noN3Wo")))

peptides <- read_tsv("conservation_protocol/rawdata/final_peptides.tsv") %>% 
  select(Proteins, contains("Intensity "))

df <- peptides %>% 
  pivot_longer(-Proteins, names_to = "sample", values_to = "Intensity") %>% 
  group_by(sample) %>% 
  summarise(sum_intensity = sum(Intensity)) %>% 
  mutate(log2 = log10(sum_intensity),
         sample = str_remove(sample, "Intensity Seifert_Hammel_Nr"),
         sample = str_remove(sample, "_WDH"),
         sample = str_remove(sample, "_STneu"),
         sample = str_remove(sample, "_wdh"),
         sample = str_remove(sample, "_neuST"),
         sample = as.numeric(sample)) %>% 
  inner_join(clean_meta, by = "sample") %>% 
  filter(treatment != "FD_noN3Wo")


model <- aov(log2 ~ treatment * incubation * substrate, data=df)

summary(model)


df %>% 
  ggplot(aes(x = incubation,
             y = log2)) +
  geom_boxplot()+
  geom_jitter(aes(colour = substrate),
              width = 0.2) +
  scale_y_continuous(limits = c(11,11.30),
                     breaks = seq(11, 11.30, 0.05)) +
  labs(x = NULL,
       y = "Peptides intensity (Log2)") +
  theme_bw()
