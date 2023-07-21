library(broom)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
TCGA_LUAD_EXPRESSION_METHYLATION_data <- read_csv("TCGA_LUAD_EXPRESSION_METHYLATION_data.csv")
lm_fit <- lm(TUBB2B ~ cg02331561, data = TCGA_LUAD_EXPRESSION_METHYLATION_data)
summary(lm_fit)
tidy(lm_fit)
lm_fit <- TCGA_LUAD_EXPRESSION_METHYLATION_data %>%
  nest(data = -submitter_id) %>%
  mutate(
    fit = map(data, ~ lm(TUBB2B ~ cg04658841, data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)


#creating a for loop
for (i in TCGA_LUAD_EXPRESSION_METHYLATION_data[c(245:2040)]){
  list <- lm(paste(TUBB2B,  '~', methylation_data[i]),
                 data = TCGA_LUAD_EXPRESSION_METHYLATION_data)
}

