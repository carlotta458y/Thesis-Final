# =======================================================================
# Accuracy and headroom analysis
#
# This file adds no health-economic modelling. Every cost, utility,
# transition probability and cohort trace it uses is produced by the engine
# in Markov_Model.R, from the same input parameter sheet. What is new
# here is arithmetic on the accuracy and compliance rate.
#
# =======================================================================


setwd("")
source("Markov_Model.R")


# Checks the NMB sign for would-leakers and non-leakers. Based on this comparison, makes the verdict
# whether using the tool and therefore classifying patients is useful -> cost-effective
# Returns  list with one element, `verdict`, one of:
#            "accuracy matters"  - the two patient types disagree, so telling
#                                  them apart creates value
#            "move all" / "move none"
#                                - the NMBs share a sign, the same decision is
#                                  right for the whole cohort, the tool adds
#                                  nothing
#            "inverted"          - the types disagree in the opposite
#                                  direction
is_accuracy_valuable <- function (NMB_wl, NMB_nl, upgrade){
  if (upgrade) {
    se_useful <- NMB_wl > 0
    sp_useful <- NMB_nl < 0
  } else {
    se_useful <- NMB_wl < 0
    sp_useful <- NMB_nl > 0
  }
  list(
    verdict = if(se_useful && sp_useful) {
      "accuracy matters"} else if(!se_useful && !sp_useful){
      "inverted"
    } else if (NMB_wl > 0){
      "whole group should receive destination surgery-> accuracy doesn't matter"
    }else{
      "whole goup should stay in source surgery -> accuracy doesn't matter"
    }
  )

}

# returns the average value (net benefit) the AI-based medical device creates per patient from the patient cohort, considering compliance
# Returns  single number - total NMB
get_total_NMB <- function (se, sp, NMB_wl, NMB_nl, wouldleaker_share, upgrade, compliance, price){
  if(upgrade){
    total_NMB <- compliance * (wouldleaker_share * se * NMB_wl + (1- wouldleaker_share) * (1-sp) * NMB_nl) # TP and FP (TP is the amount of would-leakers that get moved and FP are the amount of non leakers that get moved) -> the other stay at source so no impact
  }else{
    total_NMB <- compliance * (wouldleaker_share * (1 -se) * NMB_wl + (1- wouldleaker_share) * sp * NMB_nl) # FN and TN (FN are the amount of would-leakers that get moved and TN are the amount of non leakers that get moved)
  }
  if(price != 0){
    total_NMB - price
  }else{
    total_NMB
  }

}


# Sets the total NMB to zero and solves for either sensitivity holding specificity at its base value
# or the other way around.
# Returns  list of four:
#   se_threshold  sensitivity at which the decision breaks even, sp held at base
#   sp_threshold  specificity at which it breaks even, se held at base
#   value_se       = compliance * pi * NMB_wl
#   value_sp       = compliance * (1-pi) * (-NMB_nl)
getaccuracy_threshold <- function (se_base, sp_base, NMB_wl, NMB_nl,wouldleaker_share, upgrade, verdict, compliance){

  if(upgrade){
    se_threshold <- -((1- wouldleaker_share) * (1-sp_base) * NMB_nl)/ (wouldleaker_share  * NMB_wl)
    sp_threshold <-  1 + ((wouldleaker_share * se_base * NMB_wl) / ((1- wouldleaker_share)  * NMB_nl))
    value_se <- compliance * (wouldleaker_share  * NMB_wl) # this is the slope of the se in the total NMB formula
    value_sp <- -compliance* ((1- wouldleaker_share)  * NMB_nl)
  }else{
    se_threshold <- ( (1- wouldleaker_share) * sp_base * NMB_nl) / (wouldleaker_share  * NMB_wl) +1
    sp_threshold <- -(wouldleaker_share * (1 -se_base) * NMB_wl)/ ((1- wouldleaker_share)  * NMB_nl)
    value_se <- -compliance *(wouldleaker_share  * NMB_wl)
    value_sp <- compliance *((1- wouldleaker_share)  * NMB_nl)
  }
  if(verdict != "accuracy matters"){
      se_threshold <- NA_real_
      sp_threshold <- NA_real_
  }
  list(se_threshold = se_threshold, sp_threshold = sp_threshold, value_se = value_se, value_sp = value_sp)
}

# Returns  list of two:
#   value_compliance      net benefit at base accuracy, before compliance
#                         scaling and before price
#   compliance_threshold  compliance at which that benefit covers the price
get_compliance_value <- function (se_base, sp_base, NMB_wl, NMB_nl, wouldleaker_share, upgrade, price){
  if(upgrade){
    total_NMB <- wouldleaker_share * se_base * NMB_wl + (1 - wouldleaker_share) * (1 - sp_base) * NMB_nl
  }else{
    total_NMB <- wouldleaker_share * (1 - se_base) * NMB_wl + (1 - wouldleaker_share) * sp_base * NMB_nl

  }
  compliance_threshold <- if(total_NMB > 0) price / total_NMB else NA_real_
  list(value_compliance = total_NMB, compliance_threshold = compliance_threshold)

}


# --- Stage 4: scenarios -----------------------------------------------------

run_scenario_accuracy <- function(p) {
  rows <- list()

  for (sc in seq_len(nrow(p$scenarios))) {

    source <- p$scenarios$source[sc]
    destination <- p$scenarios$destination[sc]


    # The would-leaker division is defined by the riskier of the two arms
    scenario_pair   <- unique(c(source, destination))
    riskier_surgery <- scenario_pair[which.max(p$F5[scenario_pair])]
    wouldleaker_share <- if (p$scenarios$type[sc] == "Upgrading") {
      unname(p$F5[riskier_surgery])
    } else {
      p$scenarios$pi_wouldleak[sc]
    }
    stopifnot(wouldleaker_share > 0, wouldleaker_share < 1)


    ### Source surgery of the scenario pair -> doing the analysis ------------------------------------
    pe_source_wl <- endoleak_prob(p, source, riskier_surgery, "leaker")
    pe_source_nl <- endoleak_prob(p, source, riskier_surgery, "never")


    source_trace_wl <- run_EVAR(p, source, pe_source_wl, can_shrink = shrink_allowed(p, "leaker"))
    source_trace_nl <- run_EVAR(p, source, pe_source_nl, can_shrink = shrink_allowed(p, "never"))

    trans_cost_source_wl <- get_transition_cost(p, source, pe_source_wl, can_shrink = shrink_allowed(p, "leaker"))
     trans_cost_source_nl <- get_transition_cost(p, source, pe_source_nl, can_shrink = shrink_allowed(p, "never"))

    endoleaks_source <- get_endoleak_count(p, source, pe_source_wl, can_shrink = shrink_allowed(p, "leaker")) + get_endoleak_count(p, source, pe_source_nl, can_shrink = shrink_allowed(p, "never"))

    total_cost_source_wl <- p$cost_surgery[source] + trans_cost_source_wl + get_state_cost(p, source_trace_wl)
    total_cost_source_nl <- p$cost_surgery[source] + trans_cost_source_nl + get_state_cost(p, source_trace_nl)


    ### Destination surgery of the scenario pair-----------------------------------------

    pe_destination_wl <- endoleak_prob(p, destination, riskier_surgery, "leaker")
    pe_destination_nl <- endoleak_prob(p, destination, riskier_surgery, "never")

    destination_trace_wl <- run_EVAR(p, destination, pe_destination_wl, can_shrink = shrink_allowed(p, "leaker"))
    destination_trace_nl <- run_EVAR(p, destination, pe_destination_nl, can_shrink = shrink_allowed(p, "never"))

    trans_cost_destination_wl <- get_transition_cost(p, destination, pe_destination_wl, can_shrink = shrink_allowed(p, "leaker"))
    trans_cost_destination_nl <- get_transition_cost(p, destination, pe_destination_nl, can_shrink = shrink_allowed(p, "never"))

    endoleaks_destination <- get_endoleak_count(p, destination, pe_destination_wl, can_shrink = shrink_allowed(p, "leaker")) + get_endoleak_count(p, destination, pe_destination_nl, can_shrink = shrink_allowed(p, "never"))

    total_cost_destination_wl <- p$cost_surgery[destination] + trans_cost_destination_wl + get_state_cost(p, destination_trace_wl)
    total_cost_destination_nl <- p$cost_surgery[destination] + trans_cost_destination_nl + get_state_cost(p, destination_trace_nl)


    ### Detas ---------------------------------------------------------
    delta_cost_wl    <- total_cost_destination_wl - total_cost_source_wl
    delta_cost_nl <- total_cost_destination_nl - total_cost_source_nl
   delta_utility_wl <- get_QALY(p, destination_trace_wl) - get_QALY(p, source_trace_wl)
    delta_utility_nl <- get_QALY(p, destination_trace_nl) - get_QALY(p, source_trace_nl)
    delta_leaks   <- endoleaks_destination - endoleaks_source

    delta_cost_source <- total_cost_source_wl - total_cost_source_nl
    delta_cost_destination <- total_cost_destination_wl - total_cost_destination_nl


  isupgrade <- p$scenarios$type[sc] =="Upgrading"
  for(ln in names(p$lambdas)){

    lam <- p$lambdas[[ln]]
    NMB_wl <- lam*delta_utility_wl - delta_cost_wl # cost effectiveness of moving a would-leaker
    NMB_nl <- lam* delta_utility_nl - delta_cost_nl #cost effectiveness of moving a non leaker
    total_NMB <- get_total_NMB(p$sensitivity_base, p$specificity_base, NMB_wl , NMB_nl, wouldleaker_share, isupgrade, p$compliance, p$cost_spro )
    acc <- is_accuracy_valuable(NMB_wl, NMB_nl, isupgrade)
    acc_threshold <- getaccuracy_threshold(p$sensitivity_base, p$specificity_base, NMB_wl, NMB_nl, wouldleaker_share, isupgrade, acc$verdict, p$compliance)
    comp <- get_compliance_value(p$sensitivity_base, p$specificity_base, NMB_wl, NMB_nl, wouldleaker_share, isupgrade, p$cost_spro)
    headroom_floor <- get_total_NMB(p$sensitivity_floor, p$specificity_floor, NMB_wl, NMB_nl, wouldleaker_share, isupgrade, p$compliance_floor, p$cost_spro )
    headroom_base <- total_NMB
    headroom_max <- get_total_NMB(1,1,NMB_wl, NMB_nl, wouldleaker_share, isupgrade, 1, p$cost_spro)

    rows[[length(rows) + 1]] <- data.frame(
      lambda = lam,
      scenario = sc,
      type     = p$scenarios$type[sc],
      decision =  paste(names(p$surgery_code)[source], "->",
                        names(p$surgery_code)[destination]),
      pi       = unname(wouldleaker_share),
      dcost_wl = unname(delta_cost_wl),
      dcost_nl = unname(delta_cost_nl),
      dqaly_wl = unname(delta_utility_wl),
      dqaly_nl = unname(delta_utility_nl),
      NMB_wl = unname(NMB_wl),
      NMB_nl = unname(NMB_nl),
      total_NMB_base = unname(total_NMB),
      verdict =  acc$verdict,
      se_threshold = acc_threshold$se_threshold,
      sp_threshold = acc_threshold$sp_threshold,
      se_value = acc_threshold$value_se,
      sp_value = acc_threshold$value_sp,
    compliance_value = comp$value_compliance,
      compliance_threshold = comp$compliance_threshold,
    headroom_floor = headroom_floor,
      headroom_base = headroom_base,
    headroom_max = headroom_max,
      delta_cost_source = delta_cost_source,
    delta_cost_destination = delta_cost_destination)
  }




  }

  do.call(rbind, rows)
}


run_model_accuracy <- function(raw, overrides = list(), curve_mult = list(), time_mult = list()) {
  run_scenario_accuracy(build_params(raw, overrides, curve_mult, time_mult))
}




if (sys.nframe() == 0) {
  raw <- load_raw()
  results <- run_model_accuracy(raw)


  # --- Validation checks ---------
  checks <- list()
  add_check <- function(name, pass, detail = "") {
    checks[[length(checks) + 1]] <<- data.frame(check = name,
                                                pass = pass, detail = detail)
    cat(sprintf("[%s] %s%s\n",
                if (pass) "PASS" else "FAIL",
                name,
                if (nzchar(detail)) paste0(" - ", detail) else
                  ""))

    if (!pass) stop("VALIDATION FAILED: ", name, " - ", detail)
  }

# 1. Zero compliance => the tool changes nothing => total NMB is exactly 0
nmb_zero_compl <- mapply(
    function(nwl, nnl, pi, type)
      get_total_NMB(se = 0.85, sp = 0.85,
                    NMB_wl = nwl, NMB_nl = nnl,
                    wouldleaker_share = pi,
                    upgrade = (type == "Upgrading"),
                    compliance = 0, price = 0),
    results$NMB_wl, results$NMB_nl, results$pi, results$type)

  add_check("zero_compliance_zero_delta",
            all(abs(nmb_zero_compl) < 1e-12),
            paste("max |NMB| =",
                  signif(max(abs(nmb_zero_compl)), 3)))


  # 2. Identity decision (source == destination) => all deltas exactly 0
  raw_id <- raw
  raw_id$scenarios <- data.frame(scenario = 0, type =
  "Upgrading",
  source = "Standard",
  destination = "Standard",
  pi_wouldleak = 0.4)
  p_id       <- build_params(raw_id)
  res_id_acc <- run_scenario_accuracy(p_id)
  res_id_sim <- run_scenarios(p_id)

  id_acc <- c(res_id_acc$dcost_wl, res_id_acc$dcost_nl,
  res_id_acc$dqaly_wl, res_id_acc$dqaly_nl,
  res_id_acc$delta_leaks)
  id_sim <- c(res_id_sim$delta_cost, res_id_sim$delta_qaly,
  res_id_sim$delta_leaks_per_100)

  add_check("identity_decision_zero_delta_acc", all(abs(id_acc)
  < 1e-10),  paste("max |delta| =", signif(max(abs(id_acc)),
  3)))
  add_check("identity_decision_zero_delta_sim", all(abs(id_sim)  < 1e-10),
  paste("max |delta| =", signif(max(abs(id_sim)), 3)))






write_xlsx(results, "../Model/results_temp.xlsx")
}
