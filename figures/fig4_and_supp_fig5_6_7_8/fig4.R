library(tidyverse)

# Read data -----------------------------------------------------------------------------------------------------
df_exp_res <- 
  readr::read_csv(
    file = "exp_data.csv"
  ) %>% 
  dplyr::mutate(
    exp_result_kd = dplyr::case_when(
      .data$activity_type == "kd" & .data$activity_value <= 1000 & .data$activity_relation != ">" ~ "active",
      .data$activity_type == "kd" & .data$activity_value > 1000 ~ "inactive",
      .data$activity_type == "inhibition" ~ NA_character_
    ),
    exp_result = dplyr::case_when(
      .data$activity_type == "kd" ~ exp_result_kd,
      .data$activity_type == "inhibition" & .data$activity_value >= 75 ~ "active",
      .data$activity_type == "inhibition" & .data$activity_value < 75 ~ "inactive"
    )
  )

df_pred_active <- 
  df_exp_res %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct >= 6
  )

df_pred_inactive <- 
  df_exp_res %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct < 6
  )


# Calculate hit rates -------------------------------------------------------------------------------------------
inh_thr <- seq(from = 50, to = 95, by = 5)

df_hit_rate <- 
  purrr::map_dfr(
    inh_thr,
    function(thr){
      
      df_count <-
        df_pred_active %>%
        dplyr::mutate(
          true_class = dplyr::case_when(
            .data$activity_type == "kd" ~ exp_result_kd,
            .data$activity_type == "inhibition" & .data$activity_value >= thr ~ "active",
            .data$activity_type == "inhibition" & .data$activity_value < thr ~ "inactive"
          )
        ) %>% 
        dplyr::count(.data$true_class)
      
      hit_rate_perc <- 
        (
          ((df_count %>% 
              dplyr::filter(.data$true_class == "active") %>% 
              dplyr::pull(.data$n)) / nrow(df_pred_active) )* 100
        ) %>% 
        round(digits = 2) 
      
      
      df_res <- 
        tibble::tribble(
          ~inh_thr, ~hit_rate,
          thr, hit_rate_perc,
        )
      
      return(df_res)
    }
  )


df_pred_active__new_cmp <- 
  df_pred_active %>% 
  dplyr::filter(
    .data$nearest_train_tc < 1
  )

df_hit_rate__new_cmp <- 
  purrr::map_dfr(
    inh_thr,
    function(thr){
      
      df_count <-
        df_pred_active__new_cmp %>%
        dplyr::mutate(
          true_class = dplyr::case_when(
            .data$activity_type == "kd" ~ exp_result_kd,
            .data$activity_type == "inhibition" & .data$activity_value >= thr ~ "active",
            .data$activity_type == "inhibition" & .data$activity_value < thr ~ "inactive"
          )
        ) %>% 
        dplyr::count(.data$true_class)
      
      hit_rate_perc <- 
        (
          ((df_count %>% 
              dplyr::filter(.data$true_class == "active") %>% 
              dplyr::pull(.data$n)) / nrow(df_pred_active__new_cmp) )* 100
        ) %>% 
        round(digits = 2) 
      
      df_res <- 
        tibble::tribble(
          ~inh_thr, ~hit_rate,
          thr, hit_rate_perc
        )
      
      return(df_res)
    }
  )



# Calculate negative predictive value ----------------------------------------------------------------------------
inh_thr_tn <- seq(from = 5, to = 50, by = 5)

df_npv <- 
  purrr::map_dfr(
    inh_thr_tn,
    function(thr){
      
      df_count <-
        df_pred_inactive %>%
        dplyr::mutate(
          true_class = dplyr::case_when(
            .data$activity_type == "kd" ~ exp_result_kd,
            .data$activity_type == "inhibition" & .data$activity_value >= thr ~ "active",
            .data$activity_type == "inhibition" & .data$activity_value < thr ~ "inactive"
          )
        ) %>% 
        dplyr::count(.data$true_class)
      
      hit_rate_perc <- 
        (
          ((df_count %>% 
              dplyr::filter(.data$true_class == "inactive") %>% 
              dplyr::pull(.data$n)) / nrow(df_pred_inactive) )* 100
        ) %>% 
        round(digits = 2) 
      
      df_res <- 
        tibble::tribble(
          ~inh_thr, ~hit_rate,
          thr, hit_rate_perc
        )
      
      return(df_res)
    }
  )


df_pred_inactive__new_cmp <- 
  df_pred_inactive %>% 
  dplyr::filter(
    .data$nearest_train_tc < 1
  )

df_npv__new_cmp <- 
  purrr::map_dfr(
    inh_thr_tn,
    function(thr){
      
      df_count <-
        df_pred_inactive__new_cmp %>%
        dplyr::mutate(
          true_class = dplyr::case_when(
            .data$activity_type == "kd" ~ exp_result_kd,
            .data$activity_type == "inhibition" & .data$activity_value >= thr ~ "active",
            .data$activity_type == "inhibition" & .data$activity_value < thr ~ "inactive"
          )
        ) %>% 
        dplyr::count(.data$true_class)
      
      hit_rate_perc <- 
        (
          ((df_count %>% 
              dplyr::filter(.data$true_class == "inactive") %>% 
              dplyr::pull(.data$n)) / nrow(df_pred_inactive__new_cmp) )* 100
        ) %>% 
        round(digits = 2) 
      
      df_res <- 
        tibble::tribble(
          ~inh_thr, ~hit_rate,
          thr, hit_rate_perc
        )
      
      return(df_res)
    }
  )


pred_active_n <- nrow(df_pred_active)
pred_active_new_cmp_n <- nrow(df_pred_active__new_cmp)
pred_inactive_n <- nrow(df_pred_inactive)
pred_inactive_new_cmp_n <- nrow(df_pred_inactive__new_cmp)

df_hit_rate_final <- 
  dplyr::bind_rows(
    df_hit_rate__new_cmp %>% 
      dplyr::mutate(
        type = paste0("New compounds\n(", pred_active_new_cmp_n, " CK pairs)")
      ),
    df_hit_rate %>% 
      dplyr::mutate(
        type = paste0("All (", pred_active_n, " CK pairs)")
      )
  ) %>% 
  dplyr::mutate(
    hit_rate = .data$hit_rate %>% round(digits = 1)
  )

df_npv_final <- 
  dplyr::bind_rows(
    df_npv__new_cmp %>% 
      dplyr::mutate(
        type = paste0("New compounds\n(", pred_inactive_new_cmp_n, " CK pairs)")
      ),
    df_npv %>% 
      dplyr::mutate(
        type = paste0("All (", pred_inactive_n, " CK pairs)")
      )
  ) %>% 
  dplyr::mutate(
    hit_rate = .data$hit_rate %>% round(digits = 1)
  )


# Plot figures --------------------------------------------------------------------------------------------------
axis_txt_size = 11
axis_title_size = 12
plot_title_size = 12
legend_size = 8

pl_hr <- 
  ggplot2::ggplot() +
  geom_area(
    data = df_hit_rate, 
    aes(x = inh_thr, y = hit_rate), 
    fill = "#C3E0F0",
    alpha = 0.45
  ) +
  geom_area(
    data = df_hit_rate__new_cmp, 
    aes(x = inh_thr, y = hit_rate), 
    fill = "#7AADCA", 
    alpha = 0.35
  ) +
  geom_point(
    data = df_hit_rate_final,
    aes(x = inh_thr, y = hit_rate, color = type), 
    size = 2.2
  ) +
  geom_line(
    data = df_hit_rate_final, 
    aes(x = inh_thr, y = hit_rate, color = type, linetype = type), 
    size = 0.7
  ) + 
  scale_color_manual(
    name = "Data",
    values = c(
      "All (297 CK pairs)" = "#6EAACD", 
      "New compounds\n(142 CK pairs)" = "#25719E"
      )
    ) +
  scale_linetype_manual(
    name = "Data",
    values = c(
      "All (297 CK pairs)" = "solid", 
      "New compounds\n(142 CK pairs)" = "twodash"
    )
  ) +
  geom_text(
    data = df_hit_rate_final, 
    aes(x = inh_thr, y = hit_rate, label = paste0(round(hit_rate, digits = 0), "%")), 
    color  = "#13567D",
    vjust = -1.1, 
    hjust = 0.3, 
    size = 3.2
  ) +
  geom_segment(
    aes(x = 75, y = -Inf, xend = 75, 
        yend = df_hit_rate %>% dplyr::filter(.data$inh_thr == 75) %>% dplyr::pull(.data$hit_rate)
        ),
    color = "#949494",
    linetype = "dashed",
    size = 0.65
  ) +
  labs(
    color = "Data",
    x = "%Inhibition threshold @1000 nM",
    y = "Hit rate, %"
  ) +
  ylim(
    c(0, 50)
  ) +
  scale_x_continuous(
    breaks = df_hit_rate$inh_thr
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size),
    axis.title = element_text(size = axis_title_size),
    legend.position = c(0.85, 0.895),
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    legend.box.background = element_rect(color = "gray", linewidth = 0.5),
  )
pl_hr


pl_npv <- 
  ggplot2::ggplot() +
  geom_area(
    data = df_npv, 
    aes(x = inh_thr, y = hit_rate), 
    fill = "#E6CDD2",
    alpha = 0.35
  ) +
  geom_area(
    data = df_npv__new_cmp, 
    aes(x = inh_thr, y = hit_rate), 
    fill = "#DBD6D8",
    alpha = 0.45
  ) +
  geom_point(
    data = df_npv_final,
    aes(x = inh_thr, y = hit_rate, color = type), 
    size = 2.2
  ) +
  geom_line(
    data = df_npv_final, 
    aes(x = inh_thr, y = hit_rate, color = type, linetype = type), 
    size = 0.7
  ) + 
  scale_color_manual(
    name = "Data",
    values = c(
      "All (50 CK pairs)" = "#A96E7A", 
      "New compounds\n(32 CK pairs)" = "#9C818A"
    )
  ) +
  scale_linetype_manual(
    name = "Data",
    values = c(
      "All (50 CK pairs)" = "solid", 
      "New compounds\n(32 CK pairs)" = "twodash"
    )
  ) +
  geom_text(
    data = df_npv_final, 
    aes(x = inh_thr, y = hit_rate, label = paste0(round(hit_rate, digits = 0), "%")), 
    color  = "#362025",
    vjust = -0.9, 
    hjust = 0.3, 
    size = 3.2
  ) +
  geom_segment(
    aes(x = 25, y = -Inf, xend = 25, 
        yend = df_npv__new_cmp %>% dplyr::filter(.data$inh_thr == 25) %>% dplyr::pull(.data$hit_rate)
    ),
    color = "#949494",
    linetype = "dashed",
    size = 0.65
  ) +
  labs(
    color = "Data",
    x = "%Inhibition threshold @1000 nM",
    y = "Negative predictive value, %"
  ) +
  scale_x_continuous(
    breaks = df_npv$inh_thr
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size),
    axis.title = element_text(size = axis_title_size),
    legend.position = c(0.145, 0.895),
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    legend.box.background = element_rect(color = "gray", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  )
pl_npv


pl_tc <- 
  df_pred_active %>%
  # Tanimoto distance to the nearest training compound
  dplyr::mutate(
    nearest_train_tanimoto_dist = 1 - .data$nearest_train_tc
  ) %>% 
  ggplot2::ggplot(
    aes(y = nearest_train_tanimoto_dist, x = pwkrr_two_stage_uncertainty)
  ) +
  geom_point(
    color = "#DC8768"
  ) +
  labs(
    y = "Tanimoto distance to the\nnearest training compound", 
    x = "Model uncertainty estimate"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size),
    axis.title.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size-1)
  )
pl_tc


pl_pwkrr_all <- 
  df_pred_active %>% 
  dplyr::mutate(
    exp_result = dplyr::recode(
      .data$exp_result,
      "active" = "Active", 
      "inactive" = "Inactive"
    )
  ) %>% 
  ggplot(
    aes(
      x = pwkrr_two_stage_pred_pAct, 
      y = pwkrr_two_stage_uncertainty, 
      color = exp_result,
      shape = exp_result
    )
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    x = "Predicted pActivity",
    y = "Model uncertainty estimate",
    title = paste0(
      "All (297 CK pairs), ",
      df_hit_rate %>% dplyr::filter(.data$inh_thr == 75) %>% dplyr::pull(.data$hit_rate) %>% round(digits = 0), 
      "% hit rate (inhibition >= 75%)"
    )
  ) +
  scale_color_manual(
    name = "", 
    values = c("#258A25", "#E85D75")
  ) +
  scale_shape_manual(
    values = c("Active" = 17, "Inactive" = 16),
    name = ""
  ) +
  geom_segment(
    aes(x = 6.5, y = -Inf, xend = 6.5, yend = 0.6),
    color = "#3277BC",
    linetype = "dashed",
    size = 0.65
  ) +
  geom_segment(
    aes(x = 6.5, y = 0.6, xend = Inf, yend = 0.6),
    color = "#3277BC",
    linetype = "dashed",
    size = 0.65
  ) +
  guides(
    color = guide_legend(title = "Experimental\nresult"),
    shape = guide_legend(title = "Experimental\nresult")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size),
    axis.title = element_text(size = axis_title_size),
    plot.title = element_text(hjust = 0.5, size = plot_title_size),
    legend.position = c(0.9, 0.88),
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    legend.box.background = element_rect(color = "gray", linewidth = 0.5),
    plot.margin = margin(0.6, 0, 0.6, 0, "cm")
  )
pl_pwkrr_all


pl_pwkrr_new_cmp <- 
  df_pred_active__new_cmp %>% 
  dplyr::mutate(
    exp_result = dplyr::recode(
      .data$exp_result,
      "active" = "Active", 
      "inactive" = "Inactive"
    )
  ) %>% 
  ggplot(
    aes(
      x = pwkrr_two_stage_pred_pAct, 
      y = pwkrr_two_stage_uncertainty, 
      color = exp_result,
      shape = exp_result
    )
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    x = "Predicted pActivity",
    y = "Model uncertainty estimate",
    title = paste0(
      "New compounds (142 CK pairs), ",
      df_hit_rate__new_cmp %>% dplyr::filter(.data$inh_thr == 75) %>% dplyr::pull(.data$hit_rate) %>% round(digits = 0), 
      "% hit rate (inhibition >= 75%)"
    )
  ) +
  scale_color_manual(
    name = "", 
    values = c("#258A25", "#E85D75")
  ) +
  scale_shape_manual(
    values = c("Active" = 17, "Inactive" = 16),
    name = ""
  ) +
  geom_segment(
    aes(x = 6.5, y = -Inf, xend = 6.5, yend = 0.6),
    color = "#3277BC",
    linetype = "dashed",
    size = 0.65
  ) +
  geom_segment(
    aes(x = 6.5, y = 0.6, xend = Inf, yend = 0.6),
    color = "#3277BC",
    linetype = "dashed",
    size = 0.65
  ) +
  guides(
    color = guide_legend(title = "Experimental\nresult"),
    shape = guide_legend(title = "Experimental\nresult")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size),
    axis.title = element_text(size = axis_title_size),
    plot.title = element_text(hjust = 0.5, size = plot_title_size),
    legend.position = c(0.9, 0.88),
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    legend.box.background = element_rect(color = "gray", linewidth = 0.5),
    plot.margin = margin(0.6, 0, 0.6, 0, "cm")
  )
pl_pwkrr_new_cmp


df_pred_active %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct >= 6.5,
    .data$pwkrr_two_stage_uncertainty <= 0.6
  ) %>% 
  dplyr::count(.data$exp_result)

df_pred_active__new_cmp %>% 
  dplyr::filter(
    .data$pwkrr_two_stage_pred_pAct >= 6.5,
    .data$pwkrr_two_stage_uncertainty <= 0.6
  ) %>% 
  dplyr::count(.data$exp_result)


pl_var_kinase <- 
  df_pred_active %>%
  dplyr::filter(
    nearest_train_tc == 1
  ) %>% 
  dplyr::group_by(
    .data$hgnc_symbol
  ) %>% 
  dplyr::mutate(
    median_var = median(.data$pwkrr_two_stage_uncertainty)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    hgnc_symbol = forcats::fct_reorder(.data$hgnc_symbol, desc(median_var))
  ) %>% 
  ggplot(
    aes(x = hgnc_symbol, y = pwkrr_two_stage_uncertainty)
  ) +
  geom_violin(
    color = "#C35F2C",
    fill = "#F3CCBD",
    alpha = 0.5,
    trim = FALSE
  ) +
  stat_summary(
    fun = "median", 
    geom = "point", 
    shape = 18, 
    size = 2, 
    color = "#C35F2C") +
  labs(
    y = "Model uncertainty estimate",
    x = "Kinase"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = axis_txt_size-3.5),
    axis.text.y = element_text(size = axis_txt_size),
    axis.title = element_text(size = axis_title_size)
  )
pl_var_kinase


# Per kinase hit rate ----------------------------
kinases <- df_pred_active$hgnc_symbol %>% unique()

df_hit_rate_at_75_per_kin <- 
  purrr::map_dfr(
    kinases,
    function(kin){
      thr <- 75
      
      df_tmp <-
        df_pred_active %>%
        dplyr::filter(
          .data$hgnc_symbol == kin
        ) 
      
      df_count <- 
        df_tmp %>% 
        dplyr::count(true_class = .data$exp_result)
      
      if (!("active" %in% df_count$true_class)){
        df_count <- 
          df_count %>% 
          dplyr::bind_rows(
            tibble::tribble(
              ~true_class, ~n, 
              "active", 0
            )
          )
      }
      
      if (!("inactive" %in% df_count$true_class)){
        df_count <- 
          df_count %>% 
          dplyr::bind_rows(
            tibble::tribble(
              ~true_class, ~n, 
              "inactive", 0
            )
          )
      }
      
      hit_rate_perc <- 
        (
          ((df_count %>% 
              dplyr::filter(.data$true_class == "active") %>%
              dplyr::pull(.data$n)) / nrow(df_tmp) )* 100
        ) %>% 
        round(digits = 2) 
      
      median_var <- df_tmp$pwkrr_two_stage_uncertainty %>% median()
      
      n <- df_tmp %>% nrow()
      
      df_res <- 
        tibble::tribble(
          ~kin, ~inh_thr, ~hit_rate, ~median_var, ~N,
          kin, thr, hit_rate_perc, median_var, n
        )
      
      return(df_res)
    }
  )

df_metrics <- 
  df_hit_rate_at_75_per_kin %>%
  dplyr::summarize(
    Spearman = cor(median_var, hit_rate, method = "spearman") %>% round(digits = 2),
    Pearson = cor(median_var, hit_rate, method = "pearson") %>% round(digits = 2)
  ) 
df_metrics

pl_per_kin_hit_rate <- 
  df_hit_rate_at_75_per_kin %>% 
  ggplot2::ggplot(
    aes(x = median_var, y = hit_rate)
  ) +
  geom_point(
    color = "#337CA7",
    aes(size = N)
  ) +
  geom_text(
    aes(label = kin), 
    color = "#105176", 
    size = 3,
    nudge_y = 5
  ) +
  geom_text(
    data = df_metrics,
    aes(x = 0.603, y = 108, label = paste("Pearson = ", Pearson), hjust = 0),
    size = 3.5
  ) + 
  geom_text(
    data = df_metrics,
    aes(x = 0.603, y = 103, label = paste("Spearman = ", Spearman), hjust = 0),
    size = 3.5
  ) + 
  geom_smooth(
    method = "lm",
    alpha = 0.1,
    color = "gray",
    size = 0.6
  ) + 
  scale_y_continuous(
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    x = "Median model uncertainty estimate",
    y = "Hit rate, %"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = axis_txt_size-3.5),
    axis.title = element_text(size = axis_title_size),
    legend.position = c(0.08, 0.2),
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    legend.box.background = element_rect(color = "gray", linewidth = 0.5)
  )
pl_per_kin_hit_rate


# Final figure
pl_hr + 
pl_npv + 
pl_pwkrr_all +
pl_pwkrr_new_cmp + 
pl_per_kin_hit_rate +
(
  pl_tc + 
    pl_var_kinase +
    patchwork::plot_layout(ncol = 1) 
) +
  patchwork::plot_layout(ncol = 2)

ggsave(
  "fig4.pdf",
  height = 15,
  width = 11
)
