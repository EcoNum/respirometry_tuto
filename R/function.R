# Functions respirometry
# Copyright (c) 2018, Ph. Grosjean (Philippe.Grosjean@umons.ac.be)

correct_monitoring <- function(dates, values, calib.dates, calib.values, extrapolate = FALSE) {
  if (isTRUE(extrapolate)) rule <- 2 else rule <- 1
  # Approximate the values in the series to the manual dates
  series.values <- approx(dates, values, xout = calib.dates, rule = rule)$y
  # Calculate deltas between series measurements and manual values
  deltas <- calib.values - series.values
  corr <- data.frame(dates = calib.dates, values = calib.values,
                     measures = series.values, deltas = deltas)
  # Interpolate linearly these deltas
  all_deltas <- approx(calib.dates, deltas, xout = dates, rule = rule)$y
  # Apply the correction
  res <- values + all_deltas
  structure(res, correction = corr, dates = dates, deltas = all_deltas, class = c("corrected", "numeric"))
}

plot.corrected <- function(x, y, ...) {
  if (missing(y)) y <- x - attr(x, "deltas")
  dates <- attr(x, "dates")
  range <- range(x, y, na.rm = TRUE)
  plot(dates, y, type = "l", col = "gray", ylim = range)
  lines(dates, x, type = "l", col = "red")
  points(attr(x, "correction")$dates, attr(x, "correction")$values, col = "red")
}


# Oxygen balance on a chosen time interval from respirometer data
respirometry <- function(data, series, pos, n = 1, mass = 1, vol.respi = 1.3,
                         ref.respi = 0, main = "Variation d'oxygène en respiromètre", ...) {
  if (!is.integer(pos))
    stop("'pos' must be a vector of integers")
  if ((length(pos) %% 2) != 0)
    stop("'pos' must be a vector of even length (pairs of starts and stops)")
  if (length(pos) < n * 2)
    stop("Cannot take 'n' = ", n, " item when only ", length(pos), " elements in 'pos'")
  posn <- pos[c(n * 2 - 1, n * 2)]
  measure <- data[posn[1]:posn[2], c("Time", series)]
  names(measure) <- c("Time", "O2")

  plot(measure$Time, measure$O2, type = "l", xlab = "Temps", ylab = "[O2] (mg/L)",
       main = main, ...)

  reg <- lm(O2 ~ Time, data = measure)
  print(summary(reg))
  res <- coef(reg)["Time"] * 3600 * vol.respi

  abline(coef = coef(reg), col = "red", lwd = 2)

  res_corr <- res - ref.respi

  res_corr_mass <- res_corr / mass
  res <- data.frame(
    time0 = as.POSIXct(measure$Time[1]),
    time1 = as.POSIXct(measure$Time[nrow(measure)]),
    respi = res_corr_mass,
    ref.respi = ref.respi,
    ref = main)
  attr(res, "metadata") <- list(n = n, pos = posn, mass = mass,
                                vol.respi = vol.respi, reg = reg)
  res
}
