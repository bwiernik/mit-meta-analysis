library(tidyverse)
library(lme4)
library(metafor)
library(readxl)
library(here)
library(broom.mixed)

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
         contains("N ("),
         contains("Identity ("),
         `MIT variant`:`Orig MIT`,
         `TEST BATTERY OR QUESTIONNAIRE`:time_lapse_T1_T3) %>%
  mutate(across(where(is.character), tolower))

ma_dat_ipd <- ma_dat_ipd %>%
  select(`.id`,
         Dodgy:`#`,
         Age,
         MPO,
         Treatment,
         `MIT variant`,
         Category:`MIT:T3 (follow-up) (raw)`,
         POMP_sample_max) %>%
  mutate(across(where(is.character), tolower),
         .person = `#`
         ) %>%
  filter(! `Task / subtest` %in% c("naming (responsive ~)", "naming (confrontation ~)")) %>%
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
                                              "written comprehension"
                                              )),
         broader = factor(`BROADER CATEGORY`, levels = c("language expression",
                                                         "aphasia severity",
                                                         "communication",
                                                         "domain-general function",
                                                         "language comprehension",
                                                         "speech-motor planning")),
         unvalidated = factor(Note, levels = c("validated test", "unvalidated test"))
  )

ma_dat_ipd_summ <-
  filter(ma_dat_all,
         .id %in% ma_dat_ipd$.id)

ma_dat_ipd <-
  left_join(
    ma_dat_ipd,
    select(ma_dat_ipd_summ, .id, `Orig MIT`)
  )

ma_dat_group <-
  filter(ma_dat_all,
         ! .id %in% ma_dat_ipd$.id,
         ! .id %in% c("bonakdarpour_2003",
                      "stahl_2013",
                      "tabei_2016",
                      "lim_2013")
  ) %>%
  mutate(n_eff = n_eff(`N (group MIT)`, `N (group ctrl)`),
         gppc = gppc(`N (group MIT)`,
                     `Mn: group MIT: T1 (raw)`,
                     `Mn: group MIT: T2 (raw)`,
                     `SD: group MIT: T1 (raw)`,
                     `N (group ctrl)`,
                     `Mn:group Ctrl: T1 (raw)`,
                     `Mn:group Ctrl: T2 (raw)`,
                     `SD:group Ctrl: T1 (raw)`),
         v_gppc = v_gppc(`N (group MIT)`,
                         `Mn: group MIT: T1 (raw)`,
                         `Mn: group MIT: T2 (raw)`,
                         `SD: group MIT: T1 (raw)`,
                         `N (group ctrl)`,
                         `Mn:group Ctrl: T1 (raw)`,
                         `Mn:group Ctrl: T2 (raw)`,
                         `SD:group Ctrl: T1 (raw)`,
                         rho = .90)
         ) %>%
  left_join(categ_scheme,
            by = c(`TEST BATTERY OR QUESTIONNAIRE` = "TEST BATTERY",
                   `TASK/SUBTEST` = "TASK/SCALE/SUBTEST")) %>%
  mutate(ability = factor(ABILITY, levels = c("repetition", "auditory comprehension", "everyday communication", "naming")),
         broader = factor(`BROADER CATEGORY`, levels = c("language expression", "communication", "language comprehension")),
         unvalidated = factor(Note, levels = c("validated test", "unvalidated test"))
  )

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

dat_tabei_2016 <-
  filter(ma_dat_all,
         .id == "tabei_2016")

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

ma_rct_ability <- rma.mv(
  gppc ~ ability * unvalidated, v_gppc,
  random = ~ ability | .id,
  data = ma_dat_group,
  struct = "CS",
  test = "t"
)
ma_rct_ability_res <-
  tidy(ma_rct_ability, conf.int = TRUE) %>%
  bind_rows(
    tibble(
      term = "tau",
      type = "ran_pars",
      estimate = confint(ma_rct_ability)[[1]][[1]][2,1],
      statistic = ma_rct_ability$QE,
      df.residual = df.residual(ma_rct_ability),
      p.value = ma_rct_ability$QEp,
      conf.low = confint(ma_rct_ability)[[1]][[1]][2,2],
      conf.high = confint(ma_rct_ability)[[1]][[1]][2,3]
    ),
    tibble(
      term = "rho",
      type = "ran_pars",
      estimate = confint(ma_rct_ability)[[2]][[1]][1,1],
      conf.low = confint(ma_rct_ability)[[2]][[1]][1,2],
      conf.high = confint(ma_rct_ability)[[2]][[1]][1,3]
    )
  )

ma_rct_broader <- rma.mv(
  gppc ~ broader * unvalidated, v_gppc,
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

## For RCT, only repetition/language expression had unvalidated measures

# IPD meta-analyses
ma_ipd_overall <- lmer(
  formula = diff_z ~ 1 + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_overall_res <- tidy(ma_ipd_overall, conf.int = TRUE, conf.method = "profile")

ma_ipd_ability <- lmer(
  formula = diff_z ~ 1 + ability + unvalidated + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_ability_res <- tidy(ma_ipd_ability, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader <- lmer(
  formula = diff_z ~ 1 + broader + unvalidated + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_broader_res <- tidy(ma_ipd_broader, conf.int = TRUE, conf.method = "profile")

ma_ipd_overall_noSamplemax <- lmer(
  formula = diff_z_noSamplemax ~ 1 + (1 | .id / .person),
  data = filter(ma_dat_ipd, is.na(POMP_sample_max))
)
ma_ipd_overall_noSamplemaxl_res <- tidy(ma_ipd_overall_noSamplemax, conf.int = TRUE, conf.method = "profile")

ma_ipd_ability_noSamplemax <- lmer(
  formula = diff_z_noSamplemax ~ 1 + ability + unvalidated + (1 | .id / .person),
  data = filter(ma_dat_ipd, is.na(POMP_sample_max))
)
ma_ipd_ability_noSamplemax_res <- tidy(ma_ipd_ability_noSamplemax, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader_noSamplemax <- lmer(
  formula = diff_z_noSamplemax ~ 1 + broader + unvalidated + (1 | .id / .person),
  data = filter(ma_dat_ipd, is.na(POMP_sample_max))
)
ma_ipd_broader_noSamplemax_res <- tidy(ma_ipd_broader_noSamplemax, conf.int = TRUE, conf.method = "profile")

ma_ipd_broader_otherMods <- lmer(
  formula = diff_z ~ 1 + broader + unvalidated + MPO + Age + `Orig MIT` + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_broader_otherMods_res <- tidy(ma_ipd_broader_otherMods, conf.int = TRUE, conf.method = "profile")

ma_ipd_trained <- lmer(
  formula = diff_z ~ 1 + trained + (1 | .id / .person),
  data = filter(ma_dat_ipd, stringr::str_detect(`Test battery`, "trained")) %>%
    mutate(trained = factor(`Test battery`, levels = c("__untrained", "__trained")))
)
ma_ipd_trained_res <- tidy(ma_ipd_trained, conf.int = TRUE, conf.method = "profile")

ma_ipd_target <- lmer(
  formula = diff_z ~ 1 + `TARGET SYNDROME` + (1 | .id / .person),
  data = ma_dat_ipd
)
ma_ipd_target_res <- tidy(ma_ipd_target, conf.int = TRUE, conf.method = "profile")


write_csv(ma_dat_group, file = "data_rct.csv")
write_csv(ma_dat_ipd, file = "data_ipd.csv")

# Add MPO as covariate (months post-onset), maybe age?

# $`Haro-Martínez_2017`
# -- maybe percentiles?

# $Belin_1996
# -- drop case 2

# $`Al-Janabi_2014`
# -- two sets of stimuli for each person

# $Akanuma_2016


# Analyses:
#
# - Overall
# - Only validated tests
# -- each BroaderCategory
# - Only non-validated tests
# - Non-standard MIT protocols
# - Study type
# -- group
# -- ipd
# - Control
# -- controlled
# -- uncontrolled
