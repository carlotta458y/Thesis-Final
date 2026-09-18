setwd("")
source("Markov_Model.R")

raw  <- load_raw()
base <- run_model(raw)

LAMBDA_MID<- build_params(raw)$lambdas[["mid"]]

sens_analysis <- read_excel("../Model/parameters_v3.xlsx", sheet = "sensitivity_analysis")
names(sens_analysis) <- trimws(names(sens_analysis))

sens_analysis <- sens_analysis[!is.na(sens_analysis$parameter) & !is.na(sens_analysis$include)   & sens_analysis$include == 1, ]

# runs the model with one prameter changed
# if follow up then the whole curve it mulitplied
# if time then all the time parameter for all surgeries is multiplied
run_deterministic <- function(raw, parameter, sheet, value) {
  arg <- setNames(list(value), parameter)
  if (sheet %in% c("follow_up", "scenarios")) run_model(raw, curve_mult = arg)
  else if (sheet == "time") run_model(raw, time_mult = arg)
  else                      run_model(raw, overrides  = arg)
}

out <- list()
for(i in seq_len(nrow(sens_analysis))){
  for(border in c("low", "high")){
    value <- sens_analysis[[border]][i]
    r <- run_deterministic(raw, sens_analysis$parameter[i], sens_analysis$sheet[i], value)
    r$parameter <- sens_analysis$parameter[i]
    r$border <- border
    r$value <- value

    out[[length(out) +1 ]] <- r

  }
}

dsa_long <- do.call(rbind, out)

d <- dsa_long[dsa_long$lambda == LAMBDA_MID, ]

lo <- d[d$border == "low",  c("scenario", "parameter", "NMB")]
hi <- d[d$border == "high", c("scenario", "parameter", "NMB")]
names(lo)[3] <- "NMB_at_low"
names(hi)[3] <- "NMB_at_high"

tornado <- merge(lo, hi, by = c("scenario", "parameter"))

b <- base[base$lambda == LAMBDA_MID, c("scenario", "NMB")]
names(b)[2] <- "base"

tornado <- merge(tornado, b, by = "scenario")

tornado$span <- abs(tornado$NMB_at_high - tornado$NMB_at_low)
tornado <- tornado[order(tornado$scenario, -tornado$span), ]

tornado <- tornado[tornado$span > 1, ]

plot_tornado <- function(tornado, scenario, top_n = 10) {
  t1 <- tornado[tornado$scenario == scenario, ]
  t1 <- head(t1, top_n)
  t1 <- t1[order(t1$span), ]

  base_val <- t1$base[1]
  xr <- range(c(t1$NMB_at_low, t1$NMB_at_high, base_val))
  xr <- xr + c(-1, 1) * diff(xr) * 0.08     # 8% breathing room

  op <- par(mar = c(4.5, 11, 3, 2))
  # wide left margin for labels
  on.exit(par(op))

  plot(NA, xlim = xr, ylim = c(0.4, nrow(t1) + 0.6), yaxt = "n",
       xlab = "NMB (EUR) at lambda = 50,000", ylab = "",
       main = sprintf("Scenario %d", scenario))

  for (k in seq_len(nrow(t1))) {
    rect(min(t1$NMB_at_low[k],  base_val), k - 0.35,
         max(t1$NMB_at_low[k],  base_val), k + 0.35, col = "#9ecae1", border = NA)
    rect(min(t1$NMB_at_high[k], base_val), k - 0.35,
         max(t1$NMB_at_high[k], base_val), k + 0.35, col = "#fc9272", border = NA)
  }
  abline(v = base_val, lwd = 2)
  axis(2, at = seq_len(nrow(t1)), labels = t1$parameter, las = 1, cex.axis = 0.8)
  legend("topright", fill = c("#9ecae1", "#fc9272"), bty = "n", cex = 0.8,
         legend = c("low", "high"))
}

pdf(sprintf("DSA_tornado_%dk.pdf", LAMBDA_MID /
  1000), width = 9, height = 6, pointsize = 20)
for (sc in sort(unique(tornado$scenario))) plot_tornado(tornado, sc)
dev.off()

