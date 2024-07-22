library(tidyverse)

df_exp_pred_act <- 
  readr::read_csv(
    file = "exp_data.csv"
  ) %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct >= 6 
  )

pl_exp <- list()

pl_exp[[1]] <- 
  df_exp_pred_act %>% 
  dplyr::filter(
    .data$activity_type == "kd"
  ) %>% 
  dplyr::mutate(
    pkd = -log10(.data$activity_value * 10^-9)
  ) %>% 
  ggplot2::ggplot(
    aes(x = pkd)
  ) +
  geom_histogram(
    fill = "#25719E"
  ) +
  geom_vline(
    xintercept = 6,
    color = "#E85D75",
    linetype = "dashed"
  ) +
  xlab(expression(pK[d])) +
  ylab("Number of compound-kinase pairs") +
  theme_bw()

pl_exp[[2]] <- 
  df_exp_pred_act %>% 
  dplyr::filter(
    .data$activity_type == "inhibition"
  ) %>% 
  ggplot2::ggplot(
    aes(x = activity_value)
  ) +
  geom_histogram(
    fill = "#56A9A1"
  ) +
  geom_vline(
    xintercept = 75,
    color = "#E85D75",
    linetype = "dashed"
  ) +
  xlab("%Inhibition @1000 nM") +
  ylab("Number of compound-kinase pairs") +
  theme_bw()


pl_exp %>% 
  patchwork::wrap_plots(
    ncol = 2,
    nrow = 1
  ) + 
  patchwork::plot_annotation(
    tag_levels = "A"
  )

ggsave(
  "supp_fig5.pdf",
  width = 7.5,
  height = 4
)


df_exp_pred_act %>% 
  dplyr::filter(
    .data$activity_type == "inhibition"
  ) %>% 
  ggplot2::ggplot(
    aes(x = activity_value)
  ) +
  geom_histogram(
    fill = "#56A9A1"
  ) +
  facet_wrap(
    vars(hgnc_symbol)
  ) +
  geom_vline(
    xintercept = 75,
    color = "#E85D75",
    linetype = "dashed"
  ) +
  xlab("%Inhibition @1000 nM") +
  ylab("Number of compound-kinase pairs") +
  theme_bw()

ggsave(
  "supp_fig6.pdf"
)
