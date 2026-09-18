setwd("")
source("Markov_Model.R")

raw  <- load_raw()
base <- run_model(raw)

#since the models is probabilisitc need to ensure results don't change with every run
set.seed(42)


sens_analysis <- read_excel("../Model/parameters_v3.xlsx", sheet = "sensitivity_analysis")
names(sens_analysis) <- trimws(names(sens_analysis))

sens_analysis <- sens_analysis[!is.na(sens_analysis$parameter) & !is.na(sens_analysis$standard_error) &
                                 !is.na(sens_analysis$base_value) & !is.na(sens_analysis$include)   & sens_analysis$include == 1, ]


# ---------------------------------------------------------------------------
# Convert the sheet's mean + standard error into the
# parameters R's distribution functions expect
# ---------------------------------------------------------------------------


# draw with: rbeta(n, shape1, shape2)
# returns value between 0 and 1
# centerd around the mean
fit_beta <- function (mean, s_e){
  if (mean <= 0 || mean >= 1) stop("beta: mean must be strictly between 0 and 1")
  if (s_e <= 0) stop("beta: s_e must be positive")
  # the beta only exists if the variance stays below mean * (1 - mean)
  if (s_e^2 >= mean * (1 - mean))
    stop("beta: variance must be less than mean * (1 - mean)")
  # shape1 + shape2, the effective sample size the mean and s_e imply
  total_shape <- mean * (1 - mean) / s_e^2 - 1
  list(shape1 = mean * total_shape,
       shape2 = (1 - mean) * total_shape)
}



# draw with: rgamma(n, shape, scale)
# often right skewed -> long tail to the right
# most patients have the similar kind of costs or time values -> costs and time are symmetric > most people get similar kind
# just a few outliers -> long tail on the right
# not bound between 0 and 1
fit_gamma <- function(mean, s_e){
  if (mean <= 0) stop("gamma: mean must be positive")
  if (s_e  <= 0) stop("gamma: s_e must be positive")
  scale <- s_e^2 / mean
  list(shape = mean / scale,
       scale = scale)
}


# applied on mutliplier
# draw with: rlnorm(n, meanlog, sdlog)
fit_lnorm <- function(mean, s_e){
  if (mean <= 0) stop("lnorm: mean must be positive")
  if (s_e  <= 0) stop("lnorm: s_e must be positive")
  # variance on the log scale; meanlog is shifted by half of it so that the
  # mean on the natural scale comes out equal to `mean`
  var_log <- log((s_e^2 + mean^2) / mean^2)
  list(meanlog = log(mean) - 1/2 * var_log,
       sdlog   = sqrt(var_log))
}



# Which run_model() argument does a drawn value belong to?
# follow_up rows are multipliers on a whole curve, time rows are multpliers for all surgery types
# everything else replaces one number
destination_of <- function(sheet) {
  switch (sheet, follow_up= , scenarios = "curve_mult", time = "time_mult","overrides")
}


# Where the workbook keeps the official spelling of a parameter name, so the
# sensitivity sheet can be checked against it
known_names <- function(sheet, raw) {
  switch(sheet,
         costs       = raw$costs$parameter,
         utilities   = raw$utilities$parameter,
         event_probs = raw$events$parameter,
         general     = raw$general$parameter,
         follow_up   = raw$fup$Parameter,
         time        = raw$time$parameter,
         scenarios = names(raw$scenarios),
         stop("unknown sheet '", sheet, "'"))
}


# gets the right distribution type
fit_one <- function(distribution, mean, s_e) {
  switch(distribution,
         beta  = fit_beta (mean, s_e),
         gamma = fit_gamma(mean, s_e),
         lnorm = fit_lnorm(mean, s_e),
         stop("unknown distribution '", distribution, "'"))
}


# one random value from an already-fitted parameter
draw_random <- function(distribution, p) {
  switch(distribution,
         beta  = rbeta (1, shape1  = p$shape1,  shape2 = p$shape2),
         gamma = rgamma(1, shape   = p$shape,   scale  = p$scale),
         lnorm = rlnorm(1, meanlog = p$meanlog, sdlog  = p$sdlog),
         stop("unknown distribution '", distribution, "'"))
}



# Get the distirbution shape of every parameter
fit_all <- function(sens, raw) {

  if (anyDuplicated(sens$parameter))
    stop("duplicate parameter names in the sensitivity sheet: ",
         paste(unique(sens$parameter[duplicated(sens$parameter)]), collapse = ", "))

  fitted <- vector("list", nrow(sens))

  for (i in seq_len(nrow(sens))) {
    name         <- sens$parameter[i]
    sheet        <- sens$sheet[i]
    distribution <- sens$distribution[i]

    # build_params() silently ignores an override whose name it never looks up,
    # so a typo here would make the parameter simply not vary. Catch it now.
    if (!name %in% known_names(sheet, raw))
      stop("parameter '", name, "' is not in the '", sheet, "' sheet of the workbook")

    fitted[[i]] <- list(
      parameter    = name,
      distribution = distribution,
      arg          = destination_of(sheet),
      params       = fit_one(distribution, sens$base_value[i], sens$standard_error[i])
    )
  }

  names(fitted) <- sens$parameter
  fitted
}


# One PSA iteration's worth of inputs: a single draw for EVERY parameter at once
draw_params <- function(fitted) {
  drawn <- list(overrides = list(), curve_mult = list(), time_mult = list())
  for (f in fitted) {
    drawn[[f$arg]][[f$parameter]] <- draw_random(f$distribution, f$params)
  }
  drawn
}

#------- RUN PSA -------------

fitted <- fit_all(sens_analysis, raw)
n_sim <- 1000
results <- vector("list", n_sim)

for(i in seq_len(n_sim)){
  drawn <- draw_params(fitted)
  results[[i]]<- run_model(raw, overrides = drawn$overrides, curve_mult = drawn$curve_mult, time_mult = drawn$time_mult)
}

psa <- do.call(rbind, results)

#------ PLOT CE plane ------------

plot_ce_plane <- function(psa, scenario, lambdas = sort(unique(psa$lambda))) {
  lam_ref <- psa$lambda[psa$scenario == scenario][1]
  d <- psa[psa$scenario == scenario & psa$lambda == lam_ref, ]
  stopifnot(nrow(d) > 0)



  op <- par(mar = c(4.5, 4.5, 3, 2))
  on.exit(par(op))

  # symmetric-ish limits with breathing room
  xr <- range(c(d$delta_qaly, 0)); xr <- xr + c(-1, 1) * diff(xr) * 0.08
  yr <- range(c(d$delta_cost,   0)); yr <- yr + c(-1, 1) * diff(yr) * 0.08

  plot(d$delta_qaly, d$delta_cost, xlim = xr, ylim = yr,
       pch = 16, col = rgb(0, 0, 0, 0.35), cex = 0.6,
       xlab = "Incremental QALYs",
       ylab = "Incremental cost (EUR)",
       main = sprintf("Scenario %d", scenario))

  abline(h = 0, v = 0, col = "grey50")                 # quadrant axes
  wtp_col <- c("#2a78d6", "#eb6834", "#1baf7a")[seq_along(lambdas)]
  for (i in seq_along(lambdas))
    abline(a = 0, b = lambdas[i], col = wtp_col[i], lty = 2, lwd = 2)


  # mean point
  points(mean(d$delta_qaly), mean(d$delta_cost),
         pch = 18, col = "blue", cex = 1.8)

  # proportion cost-effective at lambda, shown in corner
  p_ce <- sapply(lambdas, function(l) mean(d$delta_cost <= l *
    d$delta_qaly))
  legend("topleft", bty = "n", cex = 0.8, lty = 2, lwd = 2, col =
    wtp_col,
         legend = sprintf("lambda = %s  ->  P(CE) = %.0f%%",
                          format(lambdas, big.mark = ","), 100 * p_ce))

}

pdf("../Model/PSA_ce_plane.pdf", width = 7, height = 6)
for (sc in sort(unique(psa$scenario))) plot_ce_plane(psa, sc)
dev.off()
