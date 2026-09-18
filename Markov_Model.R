# ---------------------------------------------------------------------------
# Markov model - Cost-effectiveness analysis of the AI-based medical device
#
# Model code for: Economic evaluation of AI-based medical devices , MSc Management of Technology,
# TU Delft, 2026. Author: Carlotta Lichtenauer.
#
# The model compares two surgery types within a decision scenario which are considered for
# treating AAA patients. The outcome of the AI-based medical device identifies patients which would
# get an endoleak or not. Based on this patients either get moved within the decision scenario
# or not. The model estimates whether using the tool and movinf flagged patients is cost-effective.
# The model runs over a life time horizon of 25 years and is conducted from a hospital payer perspective.

# INPUT: excel sheet with cost, time, follow-up and utility data
# OUTPUT: excel sheet which evaluates each ecsion scneario per WTP threshold

# The code is organised in four stages so that the whole analysis can be
# re-run with perturbed inputs for the deterministic (DSA_Markov.R) and
# probabilistic (PSA_MArkov.R) sensitivity analyses:
# Inputs      load_raw()      read the workbook once
# Parameters  build_params()  raw sheets + overrides -> everything the engine needs
# Engine      pure, take `p` as first argument
# Scenarios   run_scenarios() the scenario loop
#             run_model()     build_params + run_scenarios
#
# ---------------------------------------------------------------------------

library(readxl)
library(writexl)
setwd("")


# --- Stage 1: workbook ------------------------------------------------------

load_raw <- function(path = "../Model/parameters.xlsx") {
  list(
    fup        = read_excel(path, sheet = "follow_up"),
    general    = read_excel(path, sheet = "general"),
    life_table = read_excel(path, sheet = "life_table"),
    scenarios  = read_excel(path, sheet = "scenarios"),
    events     = read_excel(path, sheet = "event_probs"),
    costs      = read_excel(path, sheet = "costs"),
    utilities  = read_excel(path, sheet = "utilities"),
    time       = read_excel(path, sheet = "time")
  )
}


# --- Stage 2: parameters ----------------------------------------------------

# overrides  : list with all border values of the sensitivity analysis -> costs, utilities, event probs
# curve_mult : in sensitivity analysis takes the mutliplier for all follow up curves -> scales all parameters per surgery per follow up moment the same
# time_mult  : in sensitivity analysis takes the multiplier for all time values -> scales all parameters per surgery type the same
build_params <- function(raw, overrides = list(), curve_mult = list(), time_mult = list()) {


  lookup <- function(df, name, what) {
    if (!is.null(overrides[[name]])) return(overrides[[name]])
    v <- df$value[match(name, df$parameter)]
    if (is.na(v)) stop(what, " parameter not found: ", name)
    v
  }
  general_of <- function(n) lookup(raw$general,   n, "general")
  cost_of    <- function(n) lookup(raw$costs,     n, "cost")
  event_of   <- function(n) lookup(raw$events,    n, "event")
  utility_of <- function(n) lookup(raw$utilities, n, "utility")

  time_tbl <- raw$time[!is.na(raw$time$parameter), ]

  # Doesn't run when time_mult is null -> so doesn't execute in the base case
  for (pname in names(time_mult)) {
    rows <- time_tbl$parameter == pname
    if (!any(rows)) stop("time_mult parameter not in time
  sheet: ", pname)
    time_tbl$value[rows] <- time_tbl$value[rows] * time_mult[[pname]]
  }

  # gets the time value per surgery type
  time_of <- function(param, surgery = NA) {
    d <- time_tbl[time_tbl$parameter == param, ]
    if (!is.na(surgery)) d <- d[d$surgery == surgery, ]
    v <- d$value[!is.na(d$value)]
    if (length(v) != 1)
      stop("time parameter not found: ", param,
           if (!is.na(surgery)) paste0(" / ", surgery) else
             "")
    v
  }

  #Initializing the time grid of the model
  horizon <- general_of("horizon_years_scenario")
  time_data <- general_of("ref_time_wouldleak_years")
  visit_days       <- c(30, 365 * seq_len(horizon))
  k_ref <- which(visit_days == time_data * 365)
  n_cycles         <- length(visit_days)
  visit_time_years <- visit_days / 365
  cycle_len_years  <- diff(c(0, visit_time_years))
  t_mid            <- visit_time_years - cycle_len_years / 2

  # Surgery type -> integer code, used for indexing throughout
  surgery_code <- c(Standard = 1, FEVAR = 2, Endoanchor = 3, OSR = 4)
  n_surgery    <- length(surgery_code)



  ## Follow-up curves
  # Input data are CUMULATIVE probabilities of the event per visit.
  fup <- raw$fup[!is.na(raw$fup$Parameter), ]

  # Function for reading an extrapolating none existing values
  # Turns cumulative survival data into hazard function
  read_cumulative_prob_curve <- function(parameter, surgery) {
    d <- fup[fup$Parameter == parameter & fup$surgery == surgery &
               !is.na(fup$cum_value), ]
    if (nrow(d) == 0) return(rep(0, n_cycles))
    data_days <- d$days
    values    <- d$cum_value
    stopifnot(all(diff(values) >= 0), all(values >= 0), all(values < 1))
    # add a day-0 anchor so the 30-day point can be interpolated
    if (!any(data_days == 0)) {
      data_days <- c(0, data_days)
      values    <- c(0, values)
    }

    hazard_data <- -log(1 - values) # cumulative incidence F -> cumulative hazard H

    # Scale whole follow-up curves by a multiplier, applied to all surgery types
    # at once so the contrast between surgeries is preserved and only the level of the hazard moves
    mult <- curve_mult[[parameter]]
    if(!is.null(mult)) hazard_data <- hazard_data* mult
    #this only works within the time frame the moment we are outside the data time frame this is flat lining
    hazard_curve <- (approx(data_days, hazard_data, xout = visit_days, rule = 2)$y)

    n <- length(data_days)
    past_data <- visit_days > data_days[n]
    if (any(past_data)) {
      # get the slope at the last values we have of the in time data -> this is essentially the hazard ratio
      constant_hazard_ratio <- (hazard_data[n] - hazard_data[n - 1]) / (data_days[n] - data_days[n - 1])
      hazard_curve[past_data] <- hazard_data[n] + constant_hazard_ratio * (visit_days[past_data] - data_days[n])
    }
    #translate back to cumulative probability (inverse of survival)
    1 - exp(-hazard_curve)
  }

  # Build the cumulative matrix with 4 surgery type rows and columns per follow up visist
  build_curve_matrix <- function(parameter) {
    m <- matrix(0, n_surgery, n_cycles)
    for (s in names(surgery_code)) {
      m[surgery_code[[s]], ] <- read_cumulative_prob_curve(parameter, s)
    }
    m
  }

  # The cumulative data matrices
  endoleak_matrix <- build_curve_matrix("p_endoleak_cum")
  reint_cum       <- build_curve_matrix("p_reint_annual")
  death_endo_cum  <- build_curve_matrix("p_death_endo_annual")
  shrink_cum      <- build_curve_matrix("p_sac_shrinking")
  death_osr_cum <- build_curve_matrix("p_death_initial_OSR")
  rupture_cum <- build_curve_matrix("p_rupture_annual")


  # Turn the cumulative data matrices into per cycle probability matrices
  reint_matrix              <- matrix(0, n_surgery, n_cycles)
  death_endo                <- matrix(0, n_surgery, n_cycles)
  shrink_sac                <- matrix(0, n_surgery, n_cycles)
  death_osr                <- matrix(0, n_surgery, n_cycles)
  rupture_matrix          <- matrix(0, n_surgery, n_cycles)
  for (k in seq_len(n_surgery)) {
    reint_matrix[k, ]              <- cumulative_to_cycle_prob(reint_cum[k, ])
    death_endo[k, ]                <- cumulative_to_cycle_prob(death_endo_cum[k, ])
    shrink_sac[k, ]                <- cumulative_to_cycle_prob(shrink_cum[k, ])
    death_osr[k,]                 <- cumulative_to_cycle_prob(death_osr_cum[k,])
    rupture_matrix[k,]            <- cumulative_to_cycle_prob(rupture_cum[k,])

  }

  # Peri-operative mortality: single value per surgery, at 30 days
  # Read straight off cum_value -> if DSA applies then mutliplication has been applied
  p_perio <- numeric(n_surgery)
  mult_perio <- curve_mult[["p_perio"]]
  for (s in names(surgery_code)) {

    v <- fup$cum_value[fup$Parameter == "p_perio" & fup$surgery == s]
    v <- v[!is.na(v)]
    stopifnot(length(v) == 1)
    # during sensitivity analysis the perio mortality is mutliplied
    # however not the cumulative value but the hazard
    if(!is.null(mult_perio)){
      hazard <- -log( 1 - v) * mult_perio
      v <- 1 - exp(-hazard)
    }
    p_perio[surgery_code[[s]]] <- v
  }

  ## Background mortality ---------------------------------------------------
  # Visit moments are never longer than one year, so a single life-table lookup per cycle start is enough
  start_age          <- general_of("start_age")
  age_at_cycle_start <- start_age + c(0, visit_time_years[-n_cycles])
  #qx is a conditional probability
  qx <- raw$life_table$qx_annual[match(floor(age_at_cycle_start), raw$life_table$age)]
  stopifnot(!any(is.na(qx)))
  stopifnot(floor(age_at_cycle_start) ==
              floor(start_age + visit_time_years - 1e-9))

  death_base_per_cycle <- 1 - (1 - qx)^cycle_len_years
  death_base_con <- matrix(rep(death_base_per_cycle, each = n_surgery),
                           n_surgery, n_cycles)

  ## Scenarios --------------------------------------------------------------
  scenarios <- raw$scenarios
  scenarios$source      <- surgery_code[scenarios$source]
  scenarios$destination <- surgery_code[scenarios$destination]

  # The elicit pi_wouldleak values needs also be scaled during the senitivity analysis the same amount as the scale on the endoleak probability
  # Assumption that pi_wouldleak is the cumulative
  # pi_wouldleak gets its won multiplier since it is an elicit value and therefore has more uncertainty
  mult_pi <- curve_mult[["pi_wouldleak"]]
  if (!is.null(mult_pi)) {
    H <- -log(1 - scenarios$pi_wouldleak)
    H <- H * mult_pi #Again multiply the hazard not the cumulative
    scenarios$pi_wouldleak <- 1 - exp(-H)
  }

  # Decision rule of high and low ERI of the AI-based medical device
  is_upgrade <- scenarios$type == "Upgrading"
  scenarios$flagged_s    <- ifelse(is_upgrade, scenarios$destination, scenarios$source) #High ERI cases will recive the destination surgery in upgrade scenarios
  scenarios$notflagged_s <- ifelse(is_upgrade, scenarios$source,      scenarios$destination) #Low ERI cases will receive the source surgery in upgrade scenarios


  ## Costs ------------------------------------------------------------------

  # costs are based on time in oeprating room and time by the clinican (opportunity cost)
  cost_min <- cost_of("c_operating_room") + cost_of("c_clinican") / 60

  procedure_cost <- function(surgery, device) {
    time_of("time of procedure", surgery) * cost_min +
      time_of("ICU stay",      surgery) * cost_of("c_ICU") +
      time_of("Hospital stay", surgery) *
        cost_of("c_hospital") +
      device
  }

  cost_surgery <- c(
    Standard   = procedure_cost("Standard", cost_of("c_endurant")),
    FEVAR      = procedure_cost("FEVAR", cost_of("c_FEVAR")),
    Endoanchor = procedure_cost("Endoanchor", cost_of("c_endurant") + cost_of("c_helifx")),
    OSR        = procedure_cost("OSR",        0))
  stopifnot(identical(names(cost_surgery),
                      names(surgery_code)))

  # Reintervention cost
  cost_emergency_OSR <- procedure_cost("OSR emergency", 0)
  cost_emergency_reintervention <- procedure_cost("EVAR emergency", cost_of("c_endurant"))

  cost_echo    <- cost_of("c_echo")
  cost_CT      <- cost_of("c_CT")
  cost_follow_up <- cost_of("c_follow_up_appointment")


  # Follow-up intensity per sac state
  cost_fu_expanding <- time_of("f_ct_expansion")   * cost_CT + cost_follow_up
  cost_fu_stable    <- time_of("f_echo_stable")   * cost_echo + cost_follow_up
  cost_fu_shrinking <- time_of("f_echo_shrinking") * cost_echo + cost_follow_up

  states <- c("sac_shrinking", "sac_stable", "sac_expansion",
              "Rupture", "PostLeak", "OSR", "Dead")

  # State-occupancy cost vector
  state_costs <- c(
    sac_shrinking = cost_fu_shrinking,
    sac_stable    = cost_fu_stable,
    sac_expansion = cost_fu_expanding,
    Rupture = 0,  # rutpure state istself doesn't have a cost just moving from the rupture state has the reintervention cost
    PostLeak      = cost_fu_shrinking, #same follow up costs for shrinking -> every 5 years
    OSR           = cost_fu_shrinking, # OSR has less follow up moments
    Dead = 0)

  transition_cost_matrix <- matrix(0, 7, 7, dimnames = list(states, states))
  transition_cost_matrix["Rupture", "OSR"]            <- cost_emergency_OSR
  transition_cost_matrix["Rupture", "PostLeak"]       <- cost_emergency_reintervention
  #transition_cost_matrix[, "Rupture"]                 <- cost_rupture -> gettting a rupture doesn't have a cost in hospital perspective
  transition_cost_matrix["Dead", "Dead"]              <- 0

  ## Utilities
  state_utils <- c(
    sac_shrinking = utility_of("u_shrinking"),
    sac_stable    = utility_of("u_stable"),
    sac_expansion = utility_of("u_expansion"),
    Rupture       = utility_of("u_rupture"),
    PostLeak      = utility_of("u_postleak"),
    OSR           = utility_of("u_osr"),
    Dead          = utility_of("u_dead"))

  list(
    states           = states,
    n_cycles         = n_cycles,
    visit_time_years = visit_time_years,
    cycle_len_years  = cycle_len_years,
    t_mid            = t_mid,
    surgery_code     = surgery_code,
    n_surgery        = n_surgery,

    endoleak_matrix = endoleak_matrix,
    k_ref = k_ref,
    reint_matrix     = reint_matrix,
    death_endo       = death_endo,
    shrink_sac       = shrink_sac,
    death_base_con   = death_base_con,
    rupture_matrix  = rupture_matrix,
    death_osr    = death_osr,
    p_perio          = p_perio,
    r_OSR            = event_of("r_osr_conversion"),

    cost_surgery           = cost_surgery,
    cost_CT               = cost_CT,
    cost_spro              = general_of("price_spro"),
    compliance              = general_of("compliance"),
    state_costs            = state_costs,
    transition_cost_matrix = transition_cost_matrix,
    state_utils            = state_utils,

    disc_cost   = general_of("disc_cost"),
    disc_qaly   = general_of("disc_qaly"),
    sensitivity_base = general_of("sens_base"),
    specificity_base = general_of("spec_base"),
    sensitivity_floor = general_of("sens_floor"),
    specificity_floor = general_of("spec_floor"),
    compliance_floor = general_of("compliance_floor"),
    lambdas     = c(low  = general_of("lambda_low"),
                    mid  = general_of("lambda_mid"),
                    high = general_of("lambda_high")),

    scenarios = scenarios,
    F5        = endoleak_matrix[, k_ref]
  )
}


# --- Stage 3: engine --------------------------------------------------------


# Cumulative incidence -> per-cycle probability, conditional on still being at
# risk in that cycle.
# Returns  vector of the same length; entry k is the probability of the event
# in cycle k given it has not happened before k. Returns 0 where
# fewer than 1e-12 of the cohort remains at risk.
cumulative_to_cycle_prob <- function(cumulative_prob) {
  previous <- c(0, cumulative_prob[-length(cumulative_prob)])
  current  <- 1 - previous
  ifelse(current <= 1e-12, 0, (cumulative_prob - previous) / current)
}


# Transition probabilities out of every state for surgery `s` at one visit.
#   pe - the endoleak probability for THIS visit only
# Returns  7 x 7 matrix, rows = from-state, cols = to-state, dimnames p$states.
# Each row sums to 1. Cycle 1 is special: peri-operative mortality replaces background and endoleak-related death.
make_transition_matrix <- function(p, s, visit, pe, can_shrink ) {
  osr <- p$surgery_code[["OSR"]]
  d_base <- p$death_base_con[s, visit]
  # In cycle 1 peri-operative mortality replaces the background/endoleak death
  if (visit == 1) d_base <- p$p_perio[s]
  d_endo    <- p$death_endo[s, visit]
  d_osr     <-  1 - (1 - d_base) * (1 - p$death_osr[osr, visit])
  p_reint   <- p$reint_matrix[s, visit]
  p_shrink  <- if (can_shrink[visit]) p$shrink_sac[s, visit] else 0
  p_rupture <- p$rupture_matrix[s,visit]


  p_rupture_osr <- p$rupture_matrix[osr, visit]
  p_reint_osr <- p$reint_matrix[osr, visit]


  m <- matrix(0, 7, 7, dimnames = list(p$states, p$states))

  m["sac_shrinking", "sac_shrinking"] <- 1 - d_base
  m["sac_shrinking", "Dead"]          <- d_base

  m["sac_stable", "sac_shrinking"] <- (1 - d_base) * (1 - pe) * p_shrink
  m["sac_stable", "sac_expansion"] <- (1 - d_base) * pe
  m["sac_stable", "Dead"]          <- d_base
  m["sac_stable", "sac_stable"]    <- (1 - d_base) * (1 - pe) * (1 - p_shrink)

  p_die_exp <- 1 - (1 - d_base) * (1 - d_endo)
  m["sac_expansion", "PostLeak"]      <- (1 - p_die_exp) * p_reint
  m["sac_expansion", "Rupture"]       <- (1 - p_die_exp) * p_rupture
  m["sac_expansion", "Dead"]          <- p_die_exp
  m["sac_expansion", "sac_expansion"] <- (1 - p_die_exp) * (1 - p_reint - p_rupture)


  # Rupture is just a routing state -> which moves patients to emergency states -> however if a patient has died from rupture that has happened already from expaning to death
  # no one can die in this state => death because of aneurysm related prob is already applied before
  m["Rupture", "OSR"]      <- p$r_OSR
  m["Rupture", "PostLeak"] <- 1 - p$r_OSR

  m["PostLeak", "PostLeak"] <- 1 - d_base
  m["PostLeak", "Dead"]     <- d_base
  # TODO: possible re-leak after reintervention (p_releak_annual, default 0)

  #TODO: OSR can be from the get go or after emegerency -> different mortality rate
  m["OSR", "Dead"] <- d_osr
  m["OSR", "Rupture"] <- (1 - d_osr) * p_rupture_osr
  m["OSR", "PostLeak"] <- (1 - d_osr) * p_reint_osr
  m["OSR", "OSR"]  <- (1 - d_osr) * (1- p_rupture_osr - p_reint_osr)

  m["Dead", "Dead"] <- 1
  m

}

# Starting cohort distribution at t = 0.
# Returns  named vector of length 7 summing to 1. Scenarios beginning with
# open repair (s = 4) start in "OSR"; all others start in "sac_stable"
start_state <- function(p, s) {
  # Starting cohort distribution over the states. If a scenario starts in OSR this is made explicit
  st <- setNames(numeric(7), p$states)
  st[if (s == 4) "OSR" else "sac_stable"] <- 1
  st
}

# Trace matrix
# trace has n_cycles + 1 rows: one per cycle boundary, not one per cycle
# Row 1 is the t=0 starting distribution; row k+1 is the state after visit k
# The half-cycle correction averages each cycle's start and end, so it needs
# both boundaries -> trace[k] and trace[k+1] for cycle k
# Returns  trace matrix, (p$n_cycles + 1) x 7, columns named by p$states.
# Row 1 is t = 0; row k+1 is the distribution after cycle k.
run_EVAR <- function(p, s, pe, can_shrink = TRUE) {
  state <- start_state(p, s)
  trace <- matrix(0, p$n_cycles + 1, 7, dimnames = list(NULL, p$states))
  trace[1, ] <- state
  for (visit in seq_len(p$n_cycles)) {
    m <- make_transition_matrix(p, s, visit, pe[visit], can_shrink)
    state <- as.vector(state %*% m)
    trace[visit + 1, ] <- state
  }
  trace
}

# Discounted total of a state-weighted quantity over the whole trace, with
# half-cycle correction. Shared by QALYs and state costs, which differ only in
# the weights and the discount rate.
# Returns  single number, discounted to t = 0.
accumulate_state <- function(p, trace, weights, disc_rate) {
  total <- 0
  for (k in seq_len(p$n_cycles)) {
    # Half-cycle correction: people leave good states gradually during a cycle,
    # so the average of start and end describes time spent better
    average_occupancy <- (trace[k, ] + trace[k + 1, ]) / 2
    this_cycle <- sum(average_occupancy * weights)
    # scale by cycle length, then discount from the cycle midpoint
    total <- total + this_cycle * p$cycle_len_years[k] / (1 + disc_rate)^p$t_mid[k]
  }
  total
}

get_QALY       <- function(p, trace) accumulate_state(p, trace, p$state_utils, p$disc_qaly)
get_state_cost <- function(p, trace) accumulate_state(p, trace, p$state_costs, p$disc_cost)

# Costs attached to transitions (reintervention, conversion, rupture)
# Returns  single number
get_transition_cost <- function(p, s, pe, can_shrink = TRUE) {
  total <- 0
  state <- start_state(p, s)
  dead_i <- which(p$states == "Dead")
  tcm <- p$transition_cost_matrix
  tcm["sac_expansion", "PostLeak"] <- unname(p$cost_surgery[s]) # elective reintervention is always the same costs as initial surgery type
  for (visit in seq_len(p$n_cycles)) {
    m <- make_transition_matrix(p, s, visit, pe[visit], can_shrink)
    # flow[i, j] = who sits in i x their chance of moving to j
    flow  <- outer(as.vector(state), rep(1, 7)) * m
    total <- total + sum(tcm * flow) / (1 + p$disc_cost)^p$t_mid[visit]
    state <- as.vector(state %*% m)
    # Adding the 30 day CT scan cost to the patients that haven't died at that moment
    if(visit == 1 && s!= 4){
      total<- total + (1-state[dead_i])* p$cost_CT / (1 + p$disc_cost)^p$t_mid[visit]

    }

  }
  total
}

# Expected number of new endoleaks per patient: everyone crossing
# sac_stable -> sac_expansion, summed over cycles.
# Returns  single number
get_endoleak_count <- function(p, s, pe, can_shrink = TRUE) {
  state <- start_state(p, s)
  total <- 0
  for (visit in seq_len(p$n_cycles)) {
    m <- make_transition_matrix(p, s, visit, pe[visit], can_shrink)
    total <- total + state[2] * m[2, 3]
    state <- as.vector(state %*% m)
  }
  unname(total)
}

# Which cycles fall inside the observed follow-up window, i.e. up to the
# reference time the would-leaker split is defined on.
# Returns  logical vector of length p$n_cycles.
in_data_window <- function (p){
  seq_len(p$n_cycles) <= p$k_ref
}


# Cycle probability of developing an endoleak depending on whether a patient was defined as a "leaker" or "never" leaker
# WIthin the the time window -> only "leakers" can get an endoleak
# Outside the time window -> no more distinction between "leakers" and "never" leakers -> same porbbaility of group that hasn't leaked yet to get an endoleak
# Note: "never" leaker obviously higher since many "leakers" already have gotten an endoleak during the time window
# Returns  per-cycle endoleak probability, length p$n_cycles, clamped to [0, 1].
endoleak_prob <- function (p,s, risk_surgery, type = c("leaker", "never")){
  type <- match.arg(type)
  inw <- in_data_window(p)
  pe <- numeric((p$n_cycles))
  F_s <- p$endoleak_matrix[s,]

  if(type == "leaker"){
    F_ref <- p$endoleak_matrix[risk_surgery, p$k_ref] # note: for upgrade F_ref is the same as the wouldleaker_share but not the same for downgrade scenario
    stopifnot(F_ref > 0)
    pe[inw] <- cumulative_to_cycle_prob(F_s[inw] / F_ref) #Conditional probability of getting an endoleak based on the would leaker definition of the riskier ssurgery
  }

  pe[!inw] <- cumulative_to_cycle_prob(F_s)[!inw] # not a conditional anymore -> based on the whole cohort
  stopifnot(all(pe >= -1e-12), all(pe <= 1 + 1e-12))
  pmin(pmax(pe, 0), 1)

}

# once we are out of the data window also earlier defined would leakers can enter the sac shrinking state
# Returns  logical vector, length p$n_cycles - whether a patient of this type
# may enter sac_shrinking in each cycle. Would-leakers cannot shrink
# inside the data window, by definition of the split.
shrink_allowed <- function(p, type = c("leaker", "never")) {
  type <- match.arg(type)
  no_shrinking <- (type == "leaker") & in_data_window(p)
  !no_shrinking
}




# --- Stage 4: Run the model per decision scenario -----------------------------------------------------

run_scenarios <- function(p) {
  rows <- list()

  for (sc in seq_len(nrow(p$scenarios))) {

    s <- p$scenarios$source[sc]

    # Layer 1 of the model logic
    # The cohort is split into would-leakers and non-leakers before the tool is apploed
    # The would-leaker division is defined by the riskier of the two arms
    scenario_pair   <- unique(c(p$scenarios$flagged_s[sc], p$scenarios$notflagged_s[sc]))
    riskier_surgery <- scenario_pair[which.max(p$F5[scenario_pair])]
    wouldleaker_share <- if (p$scenarios$type[sc] == "Upgrading") {
      unname(p$F5[riskier_surgery])
    } else {
      p$scenarios$pi_wouldleak[sc]
    }
    stopifnot(wouldleaker_share > 0, wouldleaker_share < 1)



  # Layer 2 of the model logic
  # Computes what changes for one patient type when it is moved between
  # the two surgeries of the decision pair. That needs four cells:
  # would-leaker and non-leaker, each under the source and under the
  # destination surgery. All four are run_EVAR calls, but they are spread
  # over the two arms below


    ### Comparator: static EVAR, no tool -> patient receive source surgery ------------------------------------
    pe_static_wl <- endoleak_prob(p, s, riskier_surgery, "leaker")
    pe_static_nl <- endoleak_prob(p, s, riskier_surgery  , "never")


    static_trace <- wouldleaker_share * run_EVAR(p, s, pe_static_wl, can_shrink = shrink_allowed(p, "leaker")) + (1 - wouldleaker_share) * run_EVAR(p, s, pe_static_nl, can_shrink = shrink_allowed(p, "never"))

    trans_cost_static <- wouldleaker_share  * get_transition_cost(p, s, pe_static_wl, can_shrink = shrink_allowed(p, "leaker")) + (1 - wouldleaker_share) * get_transition_cost(p, s, pe_static_nl, can_shrink = shrink_allowed(p, "never"))

    endoleaks_old <- wouldleaker_share * get_endoleak_count(p, s, pe_static_wl, can_shrink = shrink_allowed(p, "leaker")) + (1 - wouldleaker_share) * get_endoleak_count(p, s, pe_static_nl,  can_shrink = shrink_allowed(p, "never"))

    total_cost_static <- p$cost_surgery[s] + trans_cost_static + get_state_cost(p, static_trace)

    ### Intervention: AI-based medical device -----------------------------------------

    # Lookup vectors mapping the current surgery to what each arm receives
    transition_flagged    <- c(1, 2, 3, 4)
    transition_notFlagged <- c(1, 2, 3, 4)
    transition_flagged[s]    <- p$scenarios$flagged_s[sc]
    transition_notFlagged[s] <- p$scenarios$notflagged_s[sc]

    new_s_flagged    <- transition_flagged[s]
    new_s_notflagged <- transition_notFlagged[s]

    # TP and FN are known would-leakers, so their leak probability is
    # conditional on the surgery they end up receiving.
    pe_TP <- endoleak_prob(p, new_s_flagged, riskier_surgery,    "leaker")
    pe_FN <- endoleak_prob(p, new_s_notflagged, riskier_surgery,  "leaker")
    pe_TN <- endoleak_prob(p, new_s_notflagged, riskier_surgery,  "never")
    pe_FP <- endoleak_prob(p, new_s_flagged, riskier_surgery, "never")

    # Layer 3 of the model logic
    # Sensitivity and specificity are applied to classify patients
    # Results in four classifications groups
    TP <- wouldleaker_share       * p$sensitivity_base
    FN <- wouldleaker_share       * (1 - p$sensitivity_base)
    TN <- (1 - wouldleaker_share) * p$specificity_base
    FP <- (1 - wouldleaker_share) * (1 - p$specificity_base)

    trace_TP <- run_EVAR(p, new_s_flagged,    pe_TP, can_shrink = shrink_allowed(p, "leaker"))
    trace_FN <- run_EVAR(p, new_s_notflagged, pe_FN, can_shrink = shrink_allowed(p, "leaker"))
    trace_TN <- run_EVAR(p, new_s_notflagged, pe_TN, can_shrink = shrink_allowed(p, "never"))
    trace_FP <- run_EVAR(p, new_s_flagged,   pe_FP, can_shrink = shrink_allowed(p, "never"))

    trace_SP <- TP * trace_TP + FN * trace_FN + TN * trace_TN + FP * trace_FP

    # Reintervention transition cost
    # Also for TN and FP since they can get an endoleak after the data window is over
    trans_cost_SP <-
      TP * get_transition_cost(p, new_s_flagged,    pe_TP, can_shrink = shrink_allowed(p, "leaker"))+
        FN * get_transition_cost(p, new_s_notflagged, pe_FN, can_shrink = shrink_allowed(p, "leaker"))+
    TN * get_transition_cost(p, new_s_notflagged,    pe_TN, can_shrink = shrink_allowed(p, "never"))+
     FP * get_transition_cost(p, new_s_flagged, pe_FP, can_shrink = shrink_allowed(p, "never"))

    # Surgery cost, weighted by classification group size
    cost_surgery_SP <-
      TP * p$cost_surgery[new_s_flagged]    +
      FP * p$cost_surgery[new_s_flagged]    +
      FN * p$cost_surgery[new_s_notflagged] +
      TN * p$cost_surgery[new_s_notflagged]

    endoleaks_new <-
      TP * get_endoleak_count(p, new_s_flagged,    pe_TP, can_shrink = shrink_allowed(p, "leaker")) +
        FN * get_endoleak_count(p, new_s_notflagged, pe_FN, can_shrink = shrink_allowed(p, "leaker"))+
    TN * get_endoleak_count(p, new_s_notflagged,    pe_TN, can_shrink = shrink_allowed(p, "never")) +
      FP * get_endoleak_count(p, new_s_flagged, pe_FP, can_shrink = shrink_allowed(p, "never"))

    total_cost_SP <- cost_surgery_SP + trans_cost_SP + get_state_cost(p, trace_SP)

    ### Incrementals ---------------------------------------------------------
    delta_cost    <- p$compliance* (total_cost_SP - total_cost_static)
    delta_utility <- p$compliance* (get_QALY(p, trace_SP) - get_QALY(p, static_trace))
    delta_leaks   <- p$compliance* (endoleaks_old - endoleaks_new)

    for (ln in names(p$lambdas)) {
      lam      <- p$lambdas[[ln]]
      headroom <- lam * delta_utility - delta_cost   # max justifiable price
      rows[[length(rows) + 1]] <- data.frame(
        scenario   = sc,
        type       = p$scenarios$type[sc],
        lambda     = lam,
        delta_cost = unname(delta_cost),
        delta_qaly = delta_utility,
        delta_leaks_per_100 = delta_leaks * 100,
        NMB        = unname(headroom - p$cost_spro),
        headroom   = unname(headroom))
    }
  }

  do.call(rbind, rows)
}


run_model <- function(raw, overrides = list(), curve_mult = list(), time_mult = list()) {
  run_scenarios(build_params(raw, overrides, curve_mult, time_mult))
}


# --- Base case --------------------------------------------------------------
# Guarded so that sourcing this file from DSA_Markov.R does not overwrite
# results.xlsx or run the base case twice.

if (sys.nframe() == 0) {
  raw <- load_raw()
  results <- run_model(raw)
  write_xlsx(results, "../Model/results.xlsx")
}