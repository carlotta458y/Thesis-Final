# =======================================================================
# Assuming that there is an initital allocation of patients between the surgery types in a decision scenario
#
# =======================================================================




setwd("")

source("Accuracy_Model.R")

# Get value of moving whole cohort to see if the AI-based emdical device should be applied above or below threshold
get_value_move_all <- function (NMB_wl, NMB_nl, wouldleaker_share){
  wouldleaker_share * NMB_wl + (1- wouldleaker_share) * NMB_nl
}

# Recomputing total value with initial allocation
get_total_NMB_confirming <- function (se, sp, NMB_wl, NMB_nl, wouldleaker_share, upgrade, compliance, price, q){
  if(upgrade){
    total_NMB <- (1 -q) * compliance * (wouldleaker_share * se * NMB_wl + (1- wouldleaker_share) * (1-sp) * NMB_nl) -
      q * ( compliance * (wouldleaker_share * (1 -se) * NMB_wl + (1- wouldleaker_share) * sp * NMB_nl))
  }else{
    total_NMB <- (1-q) * ( compliance * (wouldleaker_share * (1 -se) * NMB_wl + (1- wouldleaker_share) * sp * NMB_nl)) -
      q * compliance * (wouldleaker_share * se * NMB_wl + (1- wouldleaker_share) * (1-sp) * NMB_nl)
  }
  if(price != 0){
    total_NMB - price
  }else{
    total_NMB
  }

}

get_confirming_threshold <- function (se_base, sp_base, NMB_wl, NMB_nl,wouldleaker_share, upgrade, compliance){
  if(upgrade){
    q <- (wouldleaker_share * se_base * NMB_wl + (1 - wouldleaker_share)* (1-sp_base) * NMB_nl)/ (wouldleaker_share * NMB_wl + (1- wouldleaker_share) * NMB_nl)

  }else{
    q <- (wouldleaker_share * (1 -se_base) * NMB_wl + (1 - wouldleaker_share)* sp_base * NMB_nl)/ (wouldleaker_share * NMB_wl + (1- wouldleaker_share) * NMB_nl)

  }
  q
}


confirming_cell <- function (q_threshold, value_tool, value_move_all){
  use_below <- value_move_all > 0

  ifelse(is.na(q_threshold) | q_threshold < 0 | q_threshold
    > 1,
         # the value never crosses zero within the cohort, so it keeps the sign it starts with
         ifelse(value_tool > 0, "any q", "never"),
         ifelse(use_below,
                sprintf("q < %.2f", q_threshold),
                sprintf("q > %.2f", q_threshold)))
}

build_confirming_table <- function (results, p){
  lam_sorted <- sort(unique(results$lambda))
  ref        <- results[results$lambda == lam_sorted[1], ]

  tbl <- data.frame(scenario = ref$scenario,
                    type     = ref$type,
                    decision = ref$decision)
  raw <- tbl

  for(lam in lam_sorted){
    res <- results[results$lambda == lam, ]
    res <- res[match(ref$scenario, res$scenario), ]

    q_star <- mapply(get_confirming_threshold,
                     NMB_wl            = res$NMB_wl,
                     NMB_nl            = res$NMB_nl,
                     wouldleaker_share = res$pi,
                     upgrade           = res$type == "Upgrading",
                     MoreArgs = list(se_base    = p$sensitivity_base,
                                     sp_base    = p$specificity_base,
                                     compliance = p$compliance))

  # what the tool is worth when nobody is on the destination surgery yet
  value_tool <- mapply(get_total_NMB_confirming,
  NMB_wl            = res$NMB_wl,
  NMB_nl            = res$NMB_nl,
  wouldleaker_share = res$pi,
  upgrade           = res$type == "Upgrading",
  MoreArgs = list(se = p$sensitivity_base, sp = p$specificity_base,
  compliance = p$compliance, price = 0, q = 0))

  # what moving the whole cohort would be worth, as the benchmark without the tool
  value_move_all <- get_value_move_all(res$NMB_wl,
  res$NMB_nl, res$pi)

  tbl[[sprintf("WTP_%dk", lam / 1000)]] <- confirming_cell(q_star, value_tool, value_move_all)

  raw[[sprintf("q_%dk",         lam / 1000)]] <- unname(q_star)
  raw[[sprintf("tool_%dk",      lam / 1000)]] <- unname(value_tool)
  raw[[sprintf("move_all_%dk",  lam / 1000)]] <- unname(value_move_all)
  }

  list(table = tbl, raw = raw)
  }



  if (sys.nframe() == 0) {
  raw_data <- load_raw()
  p        <- build_params(raw_data)
  results  <- run_scenario_accuracy(p)

  confirming <- build_confirming_table(results, p)
  print(confirming$table)
    stopifnot(all(abs(confirming$raw$nmb_q0_20k -
                        results$total_NMB_base[results$lambda
                                                 == 20000]) < 1e-9))

    write_xlsx(confirming, "../Model/results_confirming.xlsx")

  }