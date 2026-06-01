#!/usr/bin/env Rscript

# Replicate the returns-to-education post-estimation pieces of
# Heckman, Humphries, and Veramendi (2018) in R.
#
# The original replication estimates the dynamic factor model in C++ and then
# post-processes simulated potential outcomes in Stata. This script ports the
# returns-to-education post-estimation logic for wage and PV wage to R.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(haven)
  library(tidyr)
})

repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
  starts <- unique(c(getwd(), dirname(normalizePath(script_path, mustWork = FALSE))))

  for (start in starts) {
    here <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(here, "reference", "2015166data", "replication"))) {
        return(here)
      }
      parent <- dirname(here)
      if (identical(parent, here)) break
      here <- parent
    }
  }
  stop("Could not find repository root containing reference/2015166data/replication.")
}

root <- repo_root()
rep_dir <- file.path(root, "reference", "2015166data", "replication")
out_dir <- file.path(root, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

na_if_code <- function(x, code = -9999) {
  if (is.numeric(x)) {
    x[x == code] <- NA_real_
  }
  x
}

mean_if_enough <- function(x, min_n = 1L) {
  n <- sum(!is.na(x))
  if (n < min_n) return(NA_real_)
  mean(x, na.rm = TRUE)
}

nonmissing_ntile <- function(x, n) {
  out <- rep(NA_integer_, length(x))
  ok <- !is.na(x)
  out[ok] <- dplyr::ntile(x[ok], n)
  out
}

read_observed_raw <- function() {
  var_path <- file.path(rep_dir, "ana", "hhvneo_varlist3.txt")
  raw_path <- file.path(rep_dir, "ana", "hhvneo_v3.raw")
  vars <- readLines(var_path, warn = FALSE)
  vars <- vars[nzchar(vars)]

  dat <- read.table(
    raw_path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    na.strings = "-9999",
    col.names = vars,
    check.names = FALSE
  )
  as_tibble(dat)
}

load_simulation_data <- function() {
  candidates <- c(
    file.path(rep_dir, "ana", "simulation_data.dta"),
    file.path(rep_dir, "postestimation_code", "simulation_data.dta"),
    file.path(root, "data", "simulation_data.dta")
  )
  sim_path <- candidates[file.exists(candidates)][1]
  if (is.na(sim_path)) return(NULL)

  message("Reading simulated potential outcomes from: ", sim_path)
  read_dta(sim_path) |>
    zap_labels() |>
    mutate(across(everything(), na_if_code))
}

merge_observed_covariates <- function(sim) {
  observed_path <- file.path(rep_dir, "data", "hhv-neo_v3.dta")
  needed <- c("id", "tuition4y17", "tuition4y22")
  if (!file.exists(observed_path) || all(needed[-1] %in% names(sim))) return(sim)

  observed <- read_dta(observed_path) |>
    zap_labels() |>
    mutate(across(everything(), na_if_code)) |>
    select(any_of(needed))

  merged <- left_join(sim, observed, by = "id", suffix = c("", "_observed"))
  for (v in needed[-1]) {
    observed_v <- paste0(v, "_observed")
    if (observed_v %in% names(merged)) {
      merged[[v]] <- coalesce(merged[[v]], merged[[observed_v]])
    }
  }
  select(merged, -any_of(paste0(needed[-1], "_observed")))
}

add_dynamic_schooling <- function(dat) {
  if (!"sim_fac_c" %in% names(dat) && "sim_fac0" %in% names(dat)) {
    dat$sim_fac_c <- dat$sim_fac0
  }
  if (!"sim_fac_nc" %in% names(dat) && "sim_fac1" %in% names(dat)) {
    dat$sim_fac_nc <- dat$sim_fac1
  }
  for (m in 1:4) {
    vfac <- paste0("sim_educD", m, "_Vfac")
    vfac0 <- paste0("sim_educD", m, "_Vfac0")
    vfac1 <- paste0("sim_educD", m, "_Vfac1")
    if (!vfac %in% names(dat) && all(c(vfac0, vfac1) %in% names(dat))) {
      dat[[vfac]] <- dat[[vfac0]] + dat[[vfac1]]
    }
  }

  dat |>
    mutate(
      quant_c = if_else(sim_flag == 1, nonmissing_ntile(sim_fac_c, 2), NA_integer_),
      quant_nc = if_else(sim_flag == 1, nonmissing_ntile(sim_fac_nc, 2), NA_integer_),
      dec_c = if_else(sim_flag == 1, nonmissing_ntile(sim_fac_c, 10), NA_integer_),
      dec_nc = if_else(sim_flag == 1, nonmissing_ntile(sim_fac_nc, 10), NA_integer_),
      sim_sch = case_when(
        sim_educD2 == 0 & sim_educD1 == 0 & sim_flag == 1 ~ 1,
        sim_educD2 == 1 & sim_educD1 == 0 & sim_flag == 1 ~ 2,
        sim_educD3 == 0 & sim_educD1 == 1 & sim_flag == 1 ~ 3,
        sim_educD4 == 0 & sim_educD3 == 1 & sim_educD1 == 1 ~ 4,
        sim_educD4 == 1 & sim_educD3 == 1 & sim_educD1 == 1 ~ 5,
        TRUE ~ NA_real_
      ),
      sim_D1 = sim_educD1,
      sim_D2 = if_else(sim_educD1 == 1, NA_real_, sim_educD2),
      sim_D3 = if_else(sim_educD1 == 0, NA_real_, sim_educD3),
      sim_D4 = if_else(sim_educD1 == 0 | sim_educD3 == 0, NA_real_, sim_educD4)
    ) |>
    mutate(
      sim_educD1_Vtot = if_else(sim_flag == 1, sim_educD1_Vobs + sim_educD1_Vfac, NA_real_),
      sim_educD2_Vtot = if_else(sim_flag == 1, sim_educD2_Vobs + sim_educD2_Vfac, NA_real_),
      sim_educD3_Vtot = if_else(sim_flag == 1, sim_educD3_Vobs + sim_educD3_Vfac, NA_real_),
      sim_educD4_Vtot = if_else(sim_flag == 1, sim_educD4_Vobs + sim_educD4_Vfac, NA_real_),
      sim_educD1_absnormVtot = abs(sim_educD1_Vtot) / sd(sim_educD1_Vtot, na.rm = TRUE),
      sim_educD2_absnormVtot = abs(sim_educD2_Vtot) / sd(sim_educD2_Vtot, na.rm = TRUE),
      sim_educD3_absnormVtot = abs(sim_educD3_Vtot) / sd(sim_educD3_Vtot, na.rm = TRUE),
      sim_educD4_absnormVtot = abs(sim_educD4_Vtot) / sd(sim_educD4_Vtot, na.rm = TRUE)
    )
}

add_treatment_effects <- function(dat, outcomes = c("wage", "PVwage")) {
  for (x in outcomes) {
    required <- c(
      paste0("sim_", x, "_sch", 1:5),
      "sim_educD2_shprob",
      "sim_educD3_shprob",
      "sim_educD4_shprob"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0) {
      stop("Missing variables for ", x, ": ", paste(missing, collapse = ", "))
    }

    ev4 <- paste0(x, "_EV4")
    dat[[ev4]] <- dat$sim_educD4_shprob * dat[[paste0("sim_", x, "_sch5")]] +
      (1 - dat$sim_educD4_shprob) * dat[[paste0("sim_", x, "_sch4")]]

    dat[[paste0(x, "_TE1")]] <- dat$sim_educD3_shprob * dat[[ev4]] +
      (1 - dat$sim_educD3_shprob) * dat[[paste0("sim_", x, "_sch3")]] -
      dat$sim_educD2_shprob * dat[[paste0("sim_", x, "_sch2")]] -
      (1 - dat$sim_educD2_shprob) * dat[[paste0("sim_", x, "_sch1")]]
    dat[[paste0(x, "_TE2")]] <- dat[[paste0("sim_", x, "_sch2")]] - dat[[paste0("sim_", x, "_sch1")]]
    dat[[paste0(x, "_TE3")]] <- dat[[ev4]] - dat[[paste0("sim_", x, "_sch3")]]
    dat[[paste0(x, "_TE4")]] <- dat[[paste0("sim_", x, "_sch5")]] - dat[[paste0("sim_", x, "_sch4")]]

    dat[[paste0(x, "_OPT1")]] <- -dat$sim_educD3_shprob * dat[[paste0("sim_", x, "_sch3")]] +
      dat$sim_educD3_shprob * dat[[ev4]]
    dat[[paste0(x, "_DIR1")]] <- dat[[paste0("sim_", x, "_sch3")]] -
      dat$sim_educD2_shprob * dat[[paste0("sim_", x, "_sch2")]] -
      (1 - dat$sim_educD2_shprob) * dat[[paste0("sim_", x, "_sch1")]]
    dat[[paste0(x, "_OPT3")]] <- -dat$sim_educD4_shprob * dat[[paste0("sim_", x, "_sch4")]] +
      dat$sim_educD4_shprob * dat[[paste0("sim_", x, "_sch5")]]
    dat[[paste0(x, "_DIR3")]] <- dat[[paste0("sim_", x, "_sch4")]] - dat[[paste0("sim_", x, "_sch3")]]

    dat[[paste0(x, "_OR1")]] <- dat[[paste0(x, "_OPT1")]] / dat[[paste0(x, "_TE1")]]
    dat[[paste0(x, "_OR3")]] <- dat[[paste0(x, "_OPT3")]] / dat[[paste0(x, "_TE3")]]

    dat[[paste0(x, "_TTE1")]] <- dat[[paste0("sim_", x, "_sch2")]] - dat[[paste0("sim_", x, "_sch1")]]
    dat[[paste0(x, "_TTE2")]] <- dat[[paste0("sim_", x, "_sch3")]] - dat[[paste0("sim_", x, "_sch1")]]
    dat[[paste0(x, "_TTE3")]] <- dat[[paste0("sim_", x, "_sch4")]] - dat[[paste0("sim_", x, "_sch3")]]
    dat[[paste0(x, "_TTE4")]] <- dat[[paste0("sim_", x, "_sch5")]] - dat[[paste0("sim_", x, "_sch4")]]
  }
  dat
}

summarize_treatment_effects <- function(dat, outcomes = c("wage", "PVwage"), min_n = 200L) {
  margin_labels <- c(
    "1" = "HS graduate vs dropout",
    "2" = "GED vs dropout",
    "3" = "Attend college vs HS graduate",
    "4" = "College graduate vs some college"
  )

  bind_rows(lapply(outcomes, function(outcome) {
    bind_rows(lapply(1:4, function(m) {
      te <- dat[[paste0(outcome, "_TE", m)]]
      d <- dat[[paste0("sim_D", m)]]
      amte_var <- dat[[paste0("sim_educD", m, "_absnormVtot")]]
      tibble(
        outcome = outcome,
        margin = m,
        margin_label = unname(margin_labels[as.character(m)]),
        estimand = c("ATE_star", "ATE", "TT", "TUT", "AMTE"),
        n = c(
          sum(dat$sim_flag == 1 & !is.na(te), na.rm = TRUE),
          sum(dat$sim_flag == 1 & !is.na(d) & !is.na(te), na.rm = TRUE),
          sum(dat$sim_flag == 1 & d == 1 & !is.na(te), na.rm = TRUE),
          sum(dat$sim_flag == 1 & d == 0 & !is.na(te), na.rm = TRUE),
          sum(dat$sim_flag == 1 & !is.na(d) & amte_var < 0.02 & !is.na(te), na.rm = TRUE)
        ),
        estimate = c(
          mean_if_enough(te[dat$sim_flag == 1], min_n),
          mean_if_enough(te[dat$sim_flag == 1 & !is.na(d)], min_n),
          mean_if_enough(te[dat$sim_flag == 1 & d == 1], min_n),
          mean_if_enough(te[dat$sim_flag == 1 & d == 0], min_n),
          mean_if_enough(te[dat$sim_flag == 1 & !is.na(d) & amte_var < 0.02], min_n)
        )
      )
    }))
  }))
}

summarize_option_values <- function(dat, outcomes = c("wage", "PVwage"), min_n = 200L) {
  bind_rows(lapply(outcomes, function(outcome) {
    bind_rows(lapply(c(1, 3), function(m) {
      d <- dat[[paste0("sim_D", m)]]
      amte_var <- dat[[paste0("sim_educD", m, "_absnormVtot")]]
      bind_rows(lapply(c("DIR", "OPT"), function(part) {
        x <- dat[[paste0(outcome, "_", part, m)]]
        tibble(
          outcome = outcome,
          margin = m,
          component = recode(part, DIR = "direct", OPT = "continuation"),
          estimand = c("ATE_star", "ATE", "TT", "TUT", "AMTE"),
          n = c(
            sum(dat$sim_flag == 1 & !is.na(x), na.rm = TRUE),
            sum(dat$sim_flag == 1 & !is.na(d) & !is.na(x), na.rm = TRUE),
            sum(dat$sim_flag == 1 & d == 1 & !is.na(x), na.rm = TRUE),
            sum(dat$sim_flag == 1 & d == 0 & !is.na(x), na.rm = TRUE),
            sum(dat$sim_flag == 1 & !is.na(d) & amte_var < 0.02 & !is.na(x), na.rm = TRUE)
          ),
          estimate = c(
            mean_if_enough(x[dat$sim_flag == 1], min_n),
            mean_if_enough(x[dat$sim_flag == 1 & !is.na(d)], min_n),
            mean_if_enough(x[dat$sim_flag == 1 & d == 1], min_n),
            mean_if_enough(x[dat$sim_flag == 1 & d == 0], min_n),
            mean_if_enough(x[dat$sim_flag == 1 & !is.na(d) & amte_var < 0.02], min_n)
          )
        )
      }))
    }))
  }))
}

add_prte_policy <- function(dat, seed = 2018L) {
  set.seed(seed)
  alpha3 <- -0.256
  alpha4 <- -0.016

  dat <- dat |>
    mutate(
      eps3 = if_else(sim_flag == 1, rnorm(n()), NA_real_),
      eps4 = if_else(sim_flag == 1, rnorm(n()), NA_real_),
      nu3 = if_else(sim_flag == 1, sim_educD3_Vobs + sim_educD3_Vfac0 + sim_educD3_Vfac1 + eps3, NA_real_),
      nu4 = if_else(sim_flag == 1, sim_educD4_Vobs + sim_educD4_Vfac0 + sim_educD4_Vfac1 + eps4, NA_real_),
      alt_educD3 = if_else(sim_flag == 1, as.numeric(nu3 > 0), NA_real_),
      alt_educD4 = if_else(sim_flag == 1, as.numeric(nu4 > 0), NA_real_)
    )

  sd_tuition4y17 <- sd(dat$tuition4y17[dat$sim_flag == 1], na.rm = TRUE)
  sd_tuition4y22 <- sd(dat$tuition4y22[dat$sim_flag == 1], na.rm = TRUE)

  dat |>
    mutate(
      alt_count1_educD3 = if_else(sim_flag == 1, as.numeric((nu3 - alpha3 * sd_tuition4y17) > 0), NA_real_),
      alt_count1_educD4 = if_else(sim_flag == 1, as.numeric((nu4 - alpha4 * sd_tuition4y22) > 0), NA_real_),
      prte_weight_educ3 = if_else(
        sim_flag == 1 & sim_educD1 == 1 & !is.na(alt_educD3) &
          (alt_count1_educD3 - alt_educD3) != 0,
        alt_count1_educD3 - alt_educD3,
        NA_real_
      )
    )
}

summarize_prte <- function(dat, outcomes = c("wage", "PVwage")) {
  bind_rows(lapply(outcomes, function(outcome) {
    te3 <- dat[[paste0(outcome, "_TE3")]]
    weight <- dat$prte_weight_educ3
    switchers <- dat$sim_flag == 1 & dat$sim_educD1 == 1 & !is.na(weight) & weight != 0

    tibble(
      outcome = outcome,
      policy = "1 SD tuition reduction at college-entry margin",
      n_switchers = sum(switchers, na.rm = TRUE),
      prte = mean(te3[switchers] * weight[switchers], na.rm = TRUE),
      prte_graduate = mean((te3 * weight)[switchers & dat$alt_educD4 == 1], na.rm = TRUE),
      prte_non_graduate = mean((te3 * weight)[switchers & dat$alt_educD4 == 0], na.rm = TRUE)
    )
  }))
}

plot_te_by_margin <- function(te_summary) {
  te_summary |>
    filter(estimand %in% c("ATE", "TT", "TUT", "AMTE")) |>
    ggplot(aes(x = factor(margin), y = estimate, fill = estimand)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    facet_wrap(~ outcome, scales = "free_y") +
    scale_x_discrete(labels = c("D1", "D2", "D3", "D4")) +
    labs(x = "Education decision margin", y = "Average treatment effect", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

observed_returns_fallback <- function() {
  message("simulation_data.dta was not found. Creating observed-data returns summary only.")

  dat <- read_dta(file.path(rep_dir, "data", "hhv-neo_v3.dta")) |>
    zap_labels() |>
    mutate(across(everything(), na_if_code))

  if (!"newwage" %in% names(dat) && "wage" %in% names(dat)) {
    dat <- mutate(dat, newwage = wage)
  }

  dat <- dat |>
    mutate(
      schooling = case_when(
        educ5 == 1 ~ "HS dropout",
        educ5 == 2 ~ "GED",
        educ5 == 3 ~ "HS graduate",
        educ5 == 4 ~ "Some college",
        educ5 == 5 ~ "College graduate",
        collegedeg == 1 ~ "College graduate",
        somecollege == 1 | anycollege == 1 ~ "Some college",
        ged == 1 | du_ged == 1 ~ "GED",
        hsc == 1 | du_hs == 1 ~ "HS graduate",
        TRUE ~ "HS dropout"
      ),
      schooling = factor(
        schooling,
        levels = c("HS dropout", "GED", "HS graduate", "Some college", "College graduate")
      ),
      log_wage = if_else(newwage > 0, log(newwage), NA_real_)
    )

  observed <- dat |>
    filter(!is.na(schooling)) |>
    group_by(schooling) |>
    summarise(
      n = n(),
      mean_log_wage = mean(log_wage, na.rm = TRUE),
      mean_wage = mean(newwage, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      log_wage_return_vs_hs_dropout = mean_log_wage - mean_log_wage[schooling == "HS dropout"],
      wage_return_vs_hs_dropout = mean_wage - mean_wage[schooling == "HS dropout"]
    )

  write.csv(observed, file.path(out_dir, "observed_returns_to_education_fallback.csv"), row.names = FALSE)
  message("Wrote observed-data fallback table to output/observed_returns_to_education_fallback.csv")
  invisible(observed)
}

main <- function() {
  sim <- load_simulation_data()

  if (is.null(sim)) {
    observed_returns_fallback()
    message(
      "\nFor the dynamic model tables, first generate ana/simulation_data.dta with the original ",
      "simulation step, then rerun this script."
    )
    return(invisible(NULL))
  }

  outcomes <- c("wage", "PVwage")
  sim2 <- sim |>
    merge_observed_covariates() |>
    add_dynamic_schooling() |>
    add_treatment_effects(outcomes = outcomes) |>
    add_prte_policy()

  te_summary <- summarize_treatment_effects(sim2, outcomes = outcomes)
  option_summary <- summarize_option_values(sim2, outcomes = outcomes)
  prte_summary <- summarize_prte(sim2, outcomes = outcomes)

  write.csv(te_summary, file.path(out_dir, "dynamic_returns_te_summary.csv"), row.names = FALSE)
  write.csv(option_summary, file.path(out_dir, "dynamic_returns_option_value_summary.csv"), row.names = FALSE)
  write.csv(prte_summary, file.path(out_dir, "dynamic_returns_prte_summary.csv"), row.names = FALSE)
  saveRDS(sim2, file.path(out_dir, "dynamic_returns_postestimation_data.rds"))

  ggsave(
    file.path(out_dir, "dynamic_returns_te_by_margin.png"),
    plot_te_by_margin(te_summary),
    width = 9,
    height = 5,
    dpi = 300
  )

  message("Wrote dynamic returns outputs to: ", out_dir)
  invisible(list(te = te_summary, option_values = option_summary, prte = prte_summary))
}

main()
