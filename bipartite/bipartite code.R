library(broom)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
TCGA_LUAD_EXPRESSION_METHYLATION_data <- read_csv("TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
methylation_data <- TCGA_LUAD_EXPRESSION_METHYLATION_data[c(245:720)]
for (i in methylation_data) {
   lm_fit <- lm(TUBB2B ~ i, data = TCGA_LUAD_EXPRESSION_METHYLATION_data)
}
tidy(lm_fit)
#new code
lm_fit <- lm(TUBB2B ~ cg02331561, data = TCGA_LUAD_EXPRESSION_METHYLATION_data)
summary(lm_fit)
tidy(lm_fit)
#creating a for loop
methylation_data <- TCGA_LUAD_EXPRESSION_METHYLATION_data[c(245:720)]
for (i in colnames(methylation_data)) {
  x <- i
  }
TCGA_LUAD_EXPRESSION_METHYLATION_data %>%
  nest(data = -submitter_id) %>%
  mutate(
    fit = map(data, ~ lm(TUBB2B ~ cg02001279, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)

lm_fit <- lm(TUBB2B ~ cg02331561, data = TCGA_LUAD_EXPRESSION_METHYLATION_data)
summary(lm_fit)


#grap
