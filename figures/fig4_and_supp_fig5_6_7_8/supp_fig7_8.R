library(tidyverse)

df_exp_pred_inact <- 
  readr::read_csv(
    file = "exp_data.csv"
  ) %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct < 6 
  )

df_exp_pred_inact %>% 
  ggplot2::ggplot(
    aes(x = activity_value)
  ) +
  geom_histogram(
    fill ="#ECAA84"
  ) +
  geom_vline(
    xintercept = 25,
    color = "#E85D75",
    linetype = "dashed"
  ) +
  xlab("%Inhibition @1000 nM") +
  ylab("Number of compound-kinase pairs") +
  theme_bw()

ggsave(
  "supp_fig7.pdf",
  width = 4.2,
  height = 4
)


df_exp_pred_inact %>% 
  ggplot2::ggplot(
    aes(x = activity_value)
  ) +
  geom_histogram(
    fill = "#ECAA84"
  ) +
  facet_wrap(
    vars(hgnc_symbol),
    ncol = 4
  ) +
  geom_vline(
    xintercept = 25,
    color = "#E85D75",
    linetype = "dashed"
  ) +
  xlab("%Inhibition @1000 nM") +
  ylab("Number of compound-kinase pairs") +
  theme_bw()

ggsave(
  "supp_fig8.pdf",
  width = 7.21,
  height = 4.5
)
