setwd("")

source("Accuracy_Model.R")
source("Confirming_Model.R")
raw     <- load_raw()
p       <- build_params(raw)
results <- run_scenario_accuracy(p)

lambdas <- p$lambdas



# ------- Cost-effective region in ROC space -------


# returning the accuracy surface with the accruacy threshold line
nmb_surface <- function(row, grid, price = 0) {
  outer(grid, grid, function(se, sp) {
    if (row$type == "Upgrading") {
      # compute the total value for all combinations of se and sp
      row$se_value * se - row$sp_value * (1 - sp) - price
    } else {
      -row$se_value * (1 - se) + row$sp_value * sp- price
    }
  })
}

axis_grid  <- seq(0, 1, length.out = 101)
col_ce     <- rgb(0.22, 0.55, 0.34, 0.55)
col_not_ce <- rgb(0.75, 0.22, 0.17, 0.40)
col_no_acc_ce <-  rgb(0.42, 0.60, 0.48, 0.22)
col_no_acc_not_ce <- rgb(0.72, 0.48, 0.46, 0.22)

plot_ROC <- function(row, p, price = 0) {

  surface <- nmb_surface(row, axis_grid, price)
  acc_matters <- !grepl("doesn't matter", row$verdict, fixed = TRUE)
  colours <- if(acc_matters) c(col_not_ce, col_ce) else c(col_no_acc_not_ce, col_no_acc_ce  )
  # paint the square by the sign of the total NMB
  image(axis_grid, axis_grid, ifelse(surface > 0,
                                     1, 0),
        col = colours, zlim = c(0,
                                              1), useRaster = TRUE,
        xlab = "sensitivity", ylab =
          "specificity",
        main = sprintf("S%d: %s", row$scenario,
                       row$decision))
  # the break-even line
  contour(axis_grid, axis_grid, surface, levels =
    0,
          add = TRUE, lwd = 2, drawlabels =
            FALSE)

  # SmatPlan-Pro base accuracy
  if(acc_matters){
  points(p$sensitivity_base, p$specificity_base,
         pch = 19, cex = 1.2)}
  box()
}
dir.create("../Model/plots", showWarnings = FALSE)

for(lam in lambdas){
  res    <- results[results$lambda == lam, ]
  png(file.path("../Model/plots",
                sprintf("roc_frontiers_%dk.png", lam/ 1000)),
      width = 2000, height = 1500, res = 150)

  par(mfrow = c(3, 4), mar = c(3.5, 3.5, 2.5, 1),
      mgp = c(2.2, 0.7, 0), oma = c(0, 0, 2.5, 0))

  for (i in seq_len(nrow(res)))
    plot_ROC(res[i, ], p, p$cost_spro)

  mtext(sprintf("Cost-effective region in ROC space | lambda = %dk EUR/QALY | price = %d | dot =  %d%%/%d%%",
                lam / 1000, p$cost_spro,
                100 * p$sensitivity_base, 100 *
                  p$specificity_base),
        side = 3, line = 0.8, outer = TRUE, cex =
          0.9)

  dev.off()
}

# ------- Value of accuracy (sensitivity & specificity) -------


col_se <- "#2E6E8E"
col_sp <- "#8E5E2E"

plot_accuracy_value <- function(res, lam) {
  vals <- t(as.matrix(res[, c("se_value", "sp_value")]))
  colnames(vals) <- sprintf("S%d: %s", res$scenario, res$decision)

  barplot(vals, beside = TRUE, las = 2, cex.names = 0.75, cex.axis =
    0.85,
          col = c(col_se, col_sp),
          ylab = "EUR of NMB per unit of accuracy, per patient",
          main = sprintf("Value of accuracy | lambda = %dk EUR/QALY",
                         lam / 1000))
  legend("topleft",
         c("value of sensitivity (dNMB/dse)", "value of specificity
  (dNMB/dsp)"),
         fill = c(col_se, col_sp), bty = "n")
  abline(h = 0)
}

for (lam in lambdas) {
  res <- results[results$lambda == lam, ]
  png(file.path("../Model/plots", sprintf("value_accuracy_%dk.png", lam / 1000)),
      width = 1600, height = 950, res = 150)
  par(mar = c(9, 6.5, 3, 1), mgp = c(4.5, 0.7, 0))
  plot_accuracy_value(res, lam)
  dev.off()
}


# ------- Value of compliance -------

col_lam <- c("#BFD7EA", "#5B8FB9", "#1F4E79")

plot_compliance_value <- function(results, lambdas) {
  lam_sorted <- sort(unique(lambdas))
  ref  <- results[results$lambda == lam_sorted[1], ]

  # rows = lambda, cols = scenario
  vals <- t(sapply(lam_sorted, function(lam) {
    res <- results[results$lambda == lam, ]
    res$compliance_value[match(ref$scenario, res$scenario)]
  }))
  colnames(vals) <- sprintf("S%d: %s", ref$scenario, ref$decision)

  barplot(vals, beside = TRUE, las = 2, cex.names = 0.75, cex.axis =
    0.85,
          col = col_lam[seq_along(lam_sorted)],
          ylab = "EUR of NMB per unit of compliance, per patient",
          main = "Value of compliance across WTP thresholds")
  legend("topleft", sprintf("lambda = %dk EUR/QALY", lam_sorted / 1000),
         fill = col_lam[seq_along(lam_sorted)], bty = "n")
  abline(h = 0)
}

png(file.path("../Model/plots", "value_compliance_all_lambdas.png"),
    width = 1600, height = 950, res = 150)
par(mar = c(9, 6.5, 3, 1), mgp = c(4.5, 0.7, 0))
plot_compliance_value(results, lambdas)
dev.off()


# ------- Headroom price -------

col_floor   <- "#BFD7EA"
col_base    <- "#5B8FB9"
col_max <- "#1F4E79"

plot_headroom <- function(res, p, lam) {
  hm <- t(as.matrix(res[, c("headroom_floor", "headroom_base",
                            "headroom_max")]))
  colnames(hm) <- sprintf("S%d: %s", res$scenario, res$decision)

  barplot(hm, beside = TRUE, las = 2, cex.names = 0.75, cex.axis = 0.85,
          col = c(col_floor, col_base, col_max),
          ylab = "max cost-effectiveness price per patient (EUR)",
          main = sprintf("Headroom price | WTP = %dk EUR/QALY", lam /
            1000))
  abline(h = 0)
  # abline(h = p$cost_spro, lty = 2, col = "red")   # current asking price
  legend("topleft", c("floor (min accuracy & adoption)", "base case",
                      "max (perfect & full adoption)"),
         fill = c(col_floor, col_base, col_max), bty = "n")
  # legend("topright", sprintf("asking price = %d", p$cost_spro),
  #        lty = 2, col = "red", bty = "n")
}

for (lam in lambdas) {
  res <- results[results$lambda == lam, ]
  png(file.path("../Model/plots", sprintf("headroom_%dk.png", lam / 1000)),
      width = 1600, height = 950, res = 150)
  par(mar = c(9, 6.5, 3, 1), mgp = c(4.5, 0.7, 0))
  plot_headroom(res, p, lam)
  dev.off()
}


# ------- Cost-effectiveness plane -------

# cohort-level incrementals (run_scenarios lives in Markov_Model.R)
cohort <- run_scenarios(p)
cohort$decision <- paste(names(p$surgery_code)[p$scenarios$source], "->",

names(p$surgery_code)[p$scenarios$destination])[cohort$scenario]
cohort$delta_cost_incl_price <- cohort$delta_cost + p$cost_spro

col_up   <- "#2E6E8E"
col_down <- "#8E2E4E"

  pad_range <- function(x) {          # always keep  origin in

    r <- range(c(0, x))
    r + c(-1, 1) * 0.15 * diff(r)
  }

  plot_ce_plane <- function(res, p, lam) {
    is_up <- res$type == "Upgrading"

    plot(res$delta_qaly, res$delta_cost_incl_price,
    xlim = pad_range(res$delta_qaly),
    ylim = pad_range(res$delta_cost_incl_price),
    pch = ifelse(is_up, 19, 17),
    col = ifelse(is_up, col_up, col_down), cex = 1.4,
    xlab = "incremental QALYs per patient",
    ylab = "incremental cost per patient",
    main = sprintf("Cost-effectiveness plane | lambda = %dk EUR/QALY", lam / 1000))

    mtext(sprintf("%d%% sensitivity %d%% specificity | compliance %.2f | price =  %d",
    100 * p$sensitivity_base, 100 *
    p$specificity_base,
    p$compliance, p$cost_spro), side = 3, line =
    0.3, cex = 0.85)
    abline(h = 0, v = 0, col = "grey50")
    abline(a = 0, b = lam, lty = 2, col = "grey40")   # NMB = 0,  price already in y
    text(res$delta_qaly, res$delta_cost_incl_price,
    labels = sprintf("S%d", res$scenario), pos = 3, cex =
    0.8)
    legend("topleft",
    c("Upgrading", "Downgrading",
    sprintf("break-even (NMB = 0) at %dk/QALY", lam /
    1000)),
    pch = c(19, 17, NA), lty = c(NA, NA, 2),
    col = c(col_up, col_down, "grey40"), bty = "n")
  }

  for (lam in lambdas) {
    res <- cohort[cohort$lambda == lam, ]
    png(file.path("../Model/plots", sprintf("ce_plane_%dk.png", lam /
    1000)),
    width = 1200, height = 1000, res = 150)
    par(mar = c(4.5, 4.5, 4, 1))
    plot_ce_plane(res, p, lam)
    dev.off()
  }










