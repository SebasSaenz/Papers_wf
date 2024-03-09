library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)

metadata <- read_tsv("conservation_protocol/rawdata/metadata.txt")

fr <- read_csv("conservation_protocol/rawdata/Record_FR_PCN.csv")

#------------------------------------------------------------------------------
clean_meta <- metadata %>%
  separate(
    treatment,
    into = c("sampling", "treatment", "incubation", "substrate"),
    sep = "-"
  ) %>%
  mutate(
    incubation = str_replace(incubation, "1", "Pre-treated"),
    incubation = str_replace(incubation, "0[a-b]", "No incubation"),
    incubation = factor(
      incubation,
      levels = c("No incubation", "Pre-treated", "8", "24", "72")
    ),
    substrate = if_else(is.na(substrate), "None", substrate),
    treatment = factor(
      treatment,
      levels = c("C", "W", "F", "FN", "FD", "FD_noN", "FD_noN3Wo")
    )
  )

fr_df <- fr %>%
  rename(sample = Sample) %>%
  mutate(
    sample = str_remove(sample, "Nr[0*]"),
    sample = str_remove(sample, "Nr"),
    sample = as.numeric(sample)
  ) %>%
  inner_join(clean_meta, by = "sample") %>%
  filter(treatment != "FD_noN3Wo",
         treatment != "F" | incubation != "Pre-treated",
         treatment != "F" | incubation != "8" | substrate != "079",
         treatment != "FD_noN" | incubation != "8" | substrate != "057",
         treatment != "FD" | incubation != "8" | substrate != "057") %>% 
  mutate(incubation = factor(
    incubation,
    levels = c("No incubation", "Pre-treated", "8", "24", "72")
  ))



fr_df %>%
  group_by(treatment, incubation, substrate) %>%
  identify_outliers(nFR)


model <- aov(nFR ~ treatment * incubation * substrate, data=fr_df)
summary(model)

# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


fr_df %>%
  group_by(treatment, incubation, substrate) %>%
  shapiro_test(nFR)



ggqqplot(fr_df, "nFR", ggtheme = theme_bw()) +
  facet_grid(treatment + incubation ~ substrate, labeller = "label_both")




model  <- lm(nFR ~ incubation*substrate*treatment, data = fr_df)
summary(model)

treatment.effect <- fr_df %>%
  group_by(incubation, substrate) %>%
  anova_test(nFR ~ treatment, error = model)




fr_df %>%
  group_by(treatment) %>%
  anova_test(nFR ~ incubation*substrate, error = model)

fr_df %>% 
  ggplot(aes(x = treatment,
             y = nFR)) +
  geom_boxplot(outlier.colour = "white", width = 0.4) +
  geom_jitter(aes(color = incubation,
                  shape = substrate), width = 0.3) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 0.1),
                     breaks = seq(0, 0.1, 0.01)) +
  labs(x = NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 7))

ggsave(filename = "conservation_protocol/plots/FR.png", width = 6, height = 3)
