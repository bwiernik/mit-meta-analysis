library(tidyverse)
library(lme4)
library(metafor)
library(readxl)
library(here)
library(broom.mixed)
library(patchwork)
set.seed(20160904)

n_eff <- function(n1, n2) {
  4 * (n1 * n2) / (n1 + n2)
}

sd_pooled <- function(n1, sd1, n2, sd2) {
  sqrt(
    ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  )
}

gppc <- function(n1, pre1, post1, sd1, n2, pre2, post2, sd2) {
  df <- n1 + n2 - 2
  cp <- sqrt(2 / df) * exp(lgamma(df / 2) - lgamma((df - 1) / 2))
  g <- cp * ((post1 - pre1) - (post2 - pre2)) / sd_pooled(n1, sd1, n2, sd2)
  return(g)
}

v_gppc <- function(n1, pre1, post1, sd1, n2, pre2, post2, sd2, rho = .90) {
  delta <- gppc(n1, pre1, post1, sd1, n2, pre2, post2, sd2)
  df <- n1 + n2 - 2
  cp <- sqrt(2 / df) * exp(lgamma(df / 2) - lgamma((df - 1) / 2))
  n_eff_inv <- 2 * (1 - rho) * ((n1 + n2) / (n1 * n2))
  vg <- cp^2 *
    n_eff_inv *
    ((n1 + n2 - 2) / (n1 + n2 - 4)) *
    (1 + delta^2 / n_eff_inv) -
    delta^2
  return(vg)
}

m_diff <- function(n1, pre1, post1, sd1) {
  df <- 2 * n1 - 2
  cp <- sqrt(2 / df) * exp(lgamma(df / 2) - lgamma((df - 1) / 2))
  g <- cp * (post1 - pre1) / sd1
  return(g)
}

v_m_diff <- function(n1, pre1, post1, sd1, rho = .90) {
  delta <- m_diff(n1, pre1, post1, sd1)
  df <- 2 * n1 - 2
  cp <- sqrt(2 / df) * exp(lgamma(df / 2) - lgamma((df - 1) / 2))
  vg <- cp^2 *
    2 * (1 - rho) / n1 +
    delta^2 / (2 * n1 - 2)
  return(vg)
}

ma_dat_all <-
  read_excel(here("data", "mit-meta-analsysis-data.xlsx"),
             "STUDIESALL"
  ) %>%
  filter(!is.na(`1st auth`)) %>%
  mutate(.id = paste(`1st auth`, `Year`, sep = "_"),
         .before = `Study type`)
ma_dat_ipd <-
  read_excel(here("data", "mit-meta-analsysis-data.xlsx"),
             "STUDIESIPD",
             na = c("", "NaN")
  ) %>%
  filter(!is.na(`1st auth`)) %>%
  mutate(.id = paste(`1st auth`, `Year`, sep = "_"),
         .before = `1st auth`)
categ_scheme <-
  read_excel(here("data", "mit-meta-analsysis-data.xlsx"),
             "categScheme"
  )

categ_scheme <-
  categ_scheme %>%
  mutate(across(where(is.character), tolower))

ma_dat_all <-
  ma_dat_all %>%
  select(`.id`, `Study type`,
         IPD,
         contains("MPO"),
         contains("N|"),
         contains("IdentityI"),
         `MIT variant`:`Orig MIT`,
         `Test battery`:last_col()) %>%
  mutate(across(where(is.character), tolower)) %>%
  rename(mpo_mean_groupMIT = `MPO|µ|groupMIT`,
         mpo_sd_groupMIT = `MPO|SD|groupMIT`,
         mpo_mean_groupC = `MPO|µ|groupC`,
         mpo_sd_groupC = `MPO|SD|groupC`
         )

ma_dat_ipd <- ma_dat_ipd %>%
  select(`.id`,
         Dodgy:`#`,
         Age,
         MPO,
         `MIT variant`,
         `Orig MIT`,
         `Test battery`:last_col()) %>%
  mutate(across(where(is.character), tolower),
         .person = `#`
         ) %>%
  left_join(categ_scheme,
            by = c(`Test battery` = "TEST BATTERY",
                   `Task / subtest` = "TASK/SCALE/SUBTEST")) %>%
  mutate(ability = factor(ABILITY, levels = c("repetition",
                                              "articulatory agility",
                                              "auditory comprehension",
                                              "cognitive-executive skills",
                                              "everyday communication",
                                              "grammatical form",
                                              "naming",
                                              "overall language performance",
                                              "phrase length",
                                              "speech-motor planning",
                                              "spontaneous speech",
                                              "syllable production",
                                              "verbal expression",
                                              "word production",
                                              "written comprehension",
                                              "writing"
                                              )),
         broader = factor(`BROADER CATEGORY`, levels = c("language expression",
                                                         "aphasia severity",
                                                         "communication",
                                                         "domain-general function",
                                                         "language comprehension",
                                                         "speech-motor planning"))
  ) %>%
  mutate(
    broader = forcats::fct_relabel(broader, stringr::str_to_title),
    broader = forcats::fct_recode(
      broader,
      "Non-Communicative<br>Language Expression" = "Language Expression",
      "Domain-General<br>Functioning" = "Domain-General Function",
      "Language<br>Comprehension" = "Language Comprehension",
      "Speech-Motor<br>Planning" = "Speech-Motor Planning"
    )
  ) %>%
  mutate(Note = if_else(is.na(Note) & broader == "aphasia severity",
                        "validated test",
                        Note),
         unvalidated = recode(Note,
                              "validated test" = "validated",
                              "unvalidated test" = "unvalidated")
  ) %>%
  mutate(unvalidated = factor(if_else(
    stringr::str_detect(`Test battery`, "trained"),
    paste0(unvalidated, `Test battery`),
    unvalidated),
    levels = c("validated", "unvalidated__untrained", "unvalidated__trained")
  )) %>%
  mutate(modified_mit = factor(`Orig MIT`,
                               levels = c(1, 0),
                               labels = c("Original MIT", "Modified MIT")
  ))

ma_dat_group <-
  filter(ma_dat_all,
         ! .id %in% ma_dat_ipd$.id,
         ! .id %in% c("lim_2013")
  ) %>%
  mutate(n_eff = n_eff(`N|groupMIT`, `N|groupC`),
         gppc = gppc(`N|groupMIT`,
                     `t1|raw|µ|groupMIT`,
                     `t2|raw|µ|groupMIT`,
                     `t1|raw|SD|groupMIT`,
                     `N|groupC`,
                     `t1|raw|µ|groupC`,
                     `t2|raw|µ|groupC`,
                     `t1|raw|SD|groupC`),
         v_gppc = v_gppc(`N|groupMIT`,
                         `t1|raw|µ|groupMIT`,
                         `t2|raw|µ|groupMIT`,
                         `t1|raw|SD|groupMIT`,
                         `N|groupC`,
                         `t1|raw|µ|groupC`,
                         `t2|raw|µ|groupC`,
                         `t1|raw|SD|groupC`,
                         rho = .90),
         m_diff_c = m_diff(`N|groupC`,
                           `t1|raw|µ|groupC`,
                           `t2|raw|µ|groupC`,
                           `t1|raw|SD|groupC`),
         v_m_diff_c = v_m_diff(`N|groupC`,
                               `t1|raw|µ|groupC`,
                               `t2|raw|µ|groupC`,
                               `t1|raw|SD|groupC`,
                               rho = .90)
         ) %>%
  left_join(categ_scheme,
            by = c(`Test battery` = "TEST BATTERY",
                   `Task / subtest` = "TASK/SCALE/SUBTEST")) %>%
  mutate(ability = factor(ABILITY, levels = c("repetition", "auditory comprehension", "everyday communication", "naming")),
         broader = factor(`BROADER CATEGORY`, levels = c("language expression", "communication", "language comprehension"))
  ) %>%
  mutate(
    broader = forcats::fct_relabel(broader, stringr::str_to_title),
    broader = forcats::fct_recode(
      broader,
      "Non-Communicative<br>Language Expression" = "Language Expression",
      "Language<br>Comprehension" = "Language Comprehension"
    )
  ) %>%
  mutate(Note = if_else(is.na(Note) & broader == "aphasia severity",
                        "validated test",
                        Note),
         unvalidated = recode(Note,
                              "validated test" = "validated",
                              "unvalidated test" = "unvalidated")
  ) %>%
  mutate(unvalidated = factor(if_else(
    stringr::str_detect(`Test battery`, "trained"),
    paste0(unvalidated, `Test battery`),
    unvalidated),
    levels = c("validated", "unvalidated__untrained", "unvalidated__trained")
  )) %>%
  mutate(mpo_mean = (mpo_mean_groupMIT + mpo_mean_groupC) / 2,
         mpo_diff = mpo_mean_groupMIT - mpo_mean_groupC)

pomp_z_mod_all <- ma_dat_ipd %>%
  select(.id, .person, t1_POMP) %>%
  lmer(t1_POMP ~ 1 + (1| .id / .person), data = .)
pomp_z_mn_all <- tidy(pomp_z_mod_all)$estimate[1]
pomp_z_sd_all <- tidy(pomp_z_mod_all)$estimate[2]

pomp_z_mod_noSamplemax <- ma_dat_ipd %>%
  filter(is.na(POMP_sample_max)) %>%
  select(.id, .person, t1_POMP) %>%
  lmer(t1_POMP ~ 1 + (1| .id / .person), data = .)
pomp_z_mn_noSamplemax <- tidy(pomp_z_mod_noSamplemax)$estimate[1]
pomp_z_sd_noSamplemax <- tidy(pomp_z_mod_noSamplemax)$estimate[2]

ma_dat_ipd <-
  ma_dat_ipd %>%
  mutate(t1_z = (t1_POMP - pomp_z_mn_all) / pomp_z_sd_all,
         t2_z = (t2_POMP - pomp_z_mn_all) / pomp_z_sd_all,
         t1_z_noSamplemax = (t1_POMP - pomp_z_mn_noSamplemax) / pomp_z_sd_noSamplemax,
         t2_z_noSamplemax = (t2_POMP - pomp_z_mn_noSamplemax) / pomp_z_sd_noSamplemax,
         ) %>%
  mutate(diff_z = if_else(
    .id == "naeser_1985" & ability == "auditory comprehension",
    t2_raw - t1_raw,
    NA_real_
  )) %>%
  mutate(diff_z = ifelse(
    !is.na(diff_z),
    diff_z,
    t2_z - t1_z
  )) %>%
  mutate(diff_z_noSamplemax = ifelse(
    !is.na(diff_z),
    diff_z,
    t2_z_noSamplemax - t1_z_noSamplemax
  )) %>%
  mutate(diff_z = ifelse(
    !is.na(diff_z),
    diff_z,
    diff_POMP / pomp_z_sd_all
  )) %>%
  mutate(diff_z_noSamplemax = ifelse(
    !is.na(diff_z_noSamplemax),
    diff_z_noSamplemax,
    diff_POMP / pomp_z_sd_noSamplemax
  ))

dat_lim_2013 <-
  filter(ma_dat_all,
         .id == "lim_2013")

# RCT meta-analyses
ma_rct_overall <- rma.mv(
  gppc ~ 1, v_gppc,
  random = ~ 1 | .id,
  data = ma_dat_group,
  test = "t"
)
ma_rct_overall_res <-
  tidy(ma_rct_overall, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_rct_overall)[[1]][2,1],
      statistic = ma_rct_overall$QE,
      df.residual = df.residual(ma_rct_overall),
      p.value = ma_rct_overall$QEp,
      conf.low = confint(ma_rct_overall)[[1]][2,2],
      conf.high = confint(ma_rct_overall)[[1]][2,3]
    )
  )

ma_rct_broader <- rma.mv(
  gppc ~ 0 + broader * unvalidated, v_gppc,
  random = ~ broader | .id,
  data = ma_dat_group,
  struct = "CS",
  test = "t"
)
ma_rct_broader_res <-
  tidy(ma_rct_broader, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_rct_broader)[[1]][[1]][2,1],
      statistic = ma_rct_broader$QE,
      df.residual = df.residual(ma_rct_broader),
      p.value = ma_rct_broader$QEp,
      conf.low = confint(ma_rct_broader)[[1]][[1]][2,2],
      conf.high = confint(ma_rct_broader)[[1]][[1]][2,3]
    ),
    tibble(
      term = "rho",
      type = "ran_pars",
      estimate = confint(ma_rct_broader)[[2]][[1]][1,1],
      conf.low = confint(ma_rct_broader)[[2]][[1]][1,2],
      conf.high = confint(ma_rct_broader)[[2]][[1]][1,3]
    )
  )

ma_rct_overall_mpo <- rma.mv(
  gppc ~ mpo_mean + mpo_diff, v_gppc,
  random = ~ broader | .id,
  data = ma_dat_group,
  struct = "CS",
  test = "t"
)
ma_rct_overall_mpo_res <-
  tidy(ma_rct_overall_mpo, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_rct_overall_mpo)[[1]][[1]][2,1],
      statistic = ma_rct_overall_mpo$QE,
      df.residual = df.residual(ma_rct_overall_mpo),
      p.value = ma_rct_overall_mpo$QEp,
      conf.low = confint(ma_rct_overall_mpo)[[1]][[1]][2,2],
      conf.high = confint(ma_rct_overall_mpo)[[1]][[1]][2,3]
    ),
    tibble(
      term = "rho",
      type = "ran_pars",
      estimate = confint(ma_rct_overall_mpo)[[2]][[1]][1,1],
      conf.low = confint(ma_rct_overall_mpo)[[2]][[1]][1,2],
      conf.high = confint(ma_rct_overall_mpo)[[2]][[1]][1,3]
    )
  )

ma_rct_broader_mpo <- rma.mv(
  gppc ~ 0 + mpo_mean_groupMIT + mpo_mean_groupC + broader * unvalidated, v_gppc,
  random = ~ broader | .id,
  data = ma_dat_group,
  struct = "CS",
  test = "t"
)
ma_rct_broader_mpo_res <-
  tidy(ma_rct_broader_mpo, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_rct_broader_mpo)[[1]][[1]][2,1],
      statistic = ma_rct_broader_mpo$QE,
      df.residual = df.residual(ma_rct_broader_mpo),
      p.value = ma_rct_broader_mpo$QEp,
      conf.low = confint(ma_rct_broader_mpo)[[1]][[1]][2,2],
      conf.high = confint(ma_rct_broader_mpo)[[1]][[1]][2,3]
    ),
    tibble(
      term = "rho",
      type = "ran_pars",
      estimate = confint(ma_rct_broader_mpo)[[2]][[1]][1,1],
      conf.low = confint(ma_rct_broader_mpo)[[2]][[1]][1,2],
      conf.high = confint(ma_rct_broader_mpo)[[2]][[1]][1,3]
    )
  )

## For RCT, only repetition/language expression had unvalidated measures

# IPD meta-analyses
ma_ipd_overall <- lmer(
  formula = diff_z ~ 1 + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_overall_res <- tidy(ma_ipd_overall, conf.int = TRUE, conf.method = "profile")

ma_ipd_overall_noSamplemax <- lmer(
  formula = diff_z_noSamplemax ~ 1 + (1 | .id / .person),
  data = filter(ma_dat_ipd, is.na(POMP_sample_max))
)
ma_ipd_overall_noSamplemax_res <- tidy(ma_ipd_overall_noSamplemax, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader <- lmer(
  formula = diff_z ~ 0 + broader + unvalidated + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_broader_res <- tidy(ma_ipd_broader, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader_mpo <- lmer(
  formula = diff_z ~ 0 + broader + unvalidated + MPO + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_broader_mpo_res <- tidy(ma_ipd_broader_mpo, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader_modMIT <- lmer(
  formula = diff_z ~ 0 + broader + unvalidated + modified_mit + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_broader_modMIT_res <- tidy(ma_ipd_broader_modMIT, conf.int = TRUE, conf.method = "profile")

rct_plot <- ma_dat_group %>%
  ggplot() +
  aes(x = gppc, shape = unvalidated, color = .id) +
  ggdist::stat_dist_slabinterval(
    data = cbind(
      broader = factor(levels(ma_dat_group$broader), levels = levels(ma_dat_group$broader)),
      as.data.frame(predict(ma_rct_broader,
                            newmods = matrix(c(1, 0, 0, 0, 0,
                                               0, 1, 0, 0, 0,
                                               0, 0, 1, 0, 0),
                                             byrow = TRUE, ncol = 5)))
    ),
    aes(y = fct_rev(broader),
        dist = distributional::dist_student_t(
          df.residual(ma_rct_broader), pred, se
        ),
        slab_alpha = after_stat(1 - abs(1 - 2 * cdf))
    ),
    fill = see::material_colors("indigo"),
    fill_type = "gradient",
    inherit.aes = FALSE, height = .35
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(y = stage(start = fct_rev(broader), after_stat = y - .25)),
    groupOnX = FALSE,
    cex = 3, priority = "ascending",
    na.rm = TRUE, size = 2.5
  ) +
  see::theme_modern(plot.title.space = 5) +
  guides(shape ="none", color = "none") +
  see::scale_color_material_d() +
  scale_y_discrete(name = NULL) +
  ggdist::scale_slab_alpha_continuous(range = c(.15, .9), guide = "none") +
  labs(x = NULL,
       title = "Randomised Control Trials: *g<sub>ppc</sub>* by Domain"
  ) +
  theme(plot.title = ggtext::element_markdown(face = "bold"),
        axis.text.y = ggtext::element_markdown(),
        legend.position = "bottom",
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line(size = .25)
        ) +
  coord_cartesian(xlim = c(-2, 5))

ipd_plot <- ma_dat_ipd %>%
  group_by(.id, broader, unvalidated) %>%
  summarize(gppc = mean(diff_z, na.rm = TRUE), .groups = "drop") %>%
  ggplot() +
  aes(y = fct_rev(broader), x = gppc, shape = unvalidated, color = .id) +
  ggdist::stat_dist_slabinterval(
    data = tidy(ma_ipd_broader) %>%
      filter(str_detect(term, "broader")) %>%
      mutate(broader = as_factor(str_replace(term, "broader", ""))),
    aes(y = fct_rev(broader),
        dist = distributional::dist_student_t(
          df.residual(ma_ipd_broader), estimate, std.error
        ),
        slab_alpha = after_stat(1 - abs(1 - 2 * cdf))
    ),
    fill = colorspace::darken(see::material_colors("green")),
    fill_type = "gradient",
    inherit.aes = FALSE, height = .35
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(y = stage(start = fct_rev(broader), after_stat = y - .25)),
    groupOnX = FALSE,
    cex = 1.5, priority = "ascending",
    na.rm = TRUE, size = 2.5
  ) +
  see::theme_modern(plot.title.space = 5) +
  scale_shape_discrete(labels = c("Validated", "Unvalidated (untrained)", "Unvalidated (trained)")) +
  guides(shape = guide_legend(NULL), color = "none") +
  see::scale_color_material_d() +
  scale_y_discrete(name = NULL) +
  labs(x = NULL,
       title = "Case Reports: *g<sub>pp</sub>* by Domain") +
  theme(
    plot.title = ggtext::element_markdown(face = "bold"),
    axis.text.y = ggtext::element_markdown(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(),
    panel.grid.minor.x = element_line(size = .25)
  ) +
  coord_cartesian(xlim = c(-2, 5))

plot_combo <- rct_plot / ipd_plot + plot_layout(heights = c(1, 2))

write_csv(ma_dat_group, file = "data_rct.csv")
write_csv(ma_dat_ipd, file = "data_ipd.csv")
ggsave("plots.pdf", plot_combo, device = cairo_pdf, height = 10, width = 7)
ggsave("plots.svg", plot_combo, device = svg, height = 10, width = 7)
ggsave("plots.png", plot_combo, device = png, type="cairo", dpi = 300, height = 10, width = 7)

ma_control_group_change <- rma.mv(
  m_diff_c ~ 0 + broader * unvalidated, v_m_diff_c,
  random = ~ broader | .id,
  data = ma_dat_group,
  struct = "CS",
  test = "t"
)
ma_control_group_change_res <-
  tidy(ma_control_group_change, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_control_group_change)[[1]][[1]][2,1],
      statistic = ma_control_group_change$QE,
      df.residual = df.residual(ma_control_group_change),
      p.value = ma_control_group_change$QEp,
      conf.low = confint(ma_control_group_change)[[1]][[1]][2,2],
      conf.high = confint(ma_control_group_change)[[1]][[1]][2,3]
    ),
    tibble(
      term = "rho",
      type = "ran_pars",
      estimate = confint(ma_control_group_change)[[2]][[1]][1,1],
      conf.low = confint(ma_control_group_change)[[2]][[1]][1,2],
      conf.high = confint(ma_control_group_change)[[2]][[1]][1,3]
    )
  )

ma_t1_broader_mpo <- lmer(
  formula = t1_z ~ 0 + broader + unvalidated + MPO + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_t1_broader_mpo_res <- tidy(ma_t1_broader_mpo, conf.int = TRUE, conf.method = "profile")

