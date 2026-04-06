library(tidyterra)
library(minpack.lm) 

fit_annual_curve <- function(growing_seasons, ndvi_obs, boot_int = 50, ci = .95){
  
  elmore_curve <- function(m, t) {
    m[1] + (m[2] - m[7]*t) * ((1/(1+exp((m[3]-t)/m[4]))) - (1/(1+exp((m[5]-t)/m[6]))))
  }
  
  safe_get <- function(x, l, i) {
    if (length(x) < i || length(x) > l) {
      return(NA)  # too few elements OR too many
    } else if (is.na(x[i])) {
      return(NA)
    } else {
      return(x[i])
    }
  }
  
  len_zero <- function(x){
    if ( length(x) == 0){
      return(c(NA))
      } else {
        return(x)
        }
    }

  get_transitions_doys <- function(bootstrapped_curves, doy, CI = ci, year){
  
   results <- vector("list", ncol(bootstrapped_curves))
  
  for(i in seq(ncol(bootstrapped_curves))){
    
    d1 <- diff(bootstrapped_curves[, i])              
    d2 <- diff(d1)                         
    d1_aligned <- d1[-length(d1)]          
    K      <- d2 / (1 + d1_aligned^2)^(3/2)   
    Kprime <- diff(K)                          
    Kpp    <- diff(Kprime)                     
    doy_Kp <- doy[1:length(Kprime)]
    
    # Local extrema of K' are where K'' changes sign
    sc         <- which(diff(sign(Kpp)) != 0)
    directions <- diff(sign(Kpp))[sc]
    
    # Local maxima of K' (K'' goes + to -) → Greenup, Maturity 
    maxima_doys <- doy_Kp[sc[directions < 0]]
    # Local minima of K' (K'' goes - to +) → Senescence, Dormancy
    minima_doys <- doy_Kp[sc[directions > 0]]
      
    doy_d2 <- doy[1:length(d2)]
    concavity_doys <- doy_d2[which(diff(sign(d2)) != 0)]
    
    concavity_doys <- len_zero(concavity_doys)
    maxima_doys <- len_zero(maxima_doys)
    minima_doys <- len_zero(minima_doys)
    
    peak_doy_boot <- doy[which.max(bootstrapped_curves[, i])]
    
    concavity_spring <- len_zero(concavity_doys[concavity_doys < peak_doy])
    concavity_fall   <- len_zero(concavity_doys[concavity_doys > peak_doy])

    maxima_spring <- len_zero(maxima_doys[maxima_doys < peak_doy])
    minima_fall   <- len_zero(minima_doys[minima_doys > peak_doy])
    
    data.frame(
      Greenup_doy = safe_get(maxima_spring, 2, 1),
      SOS_doy = safe_get(concavity_spring, 1, 1),
      Maturity_doy = safe_get(maxima_spring, 2, 2),
      Senescence_doy = safe_get(minima_fall, 2, 1),
      EOS_doy = safe_get(concavity_fall, 1, 1),
      Dormancy_doy = safe_get(minima_fall, 2, 2),
      Peak_NDVI_doy = safe_get(peak_doy_boot, 1, 1)) -> results[[i]]
  }
   
   do.call(rbind, results) %>%
  summarise(
    across(
      c(Greenup_doy, SOS_doy, Maturity_doy, Senescence_doy,
        EOS_doy, Dormancy_doy, Peak_NDVI_doy),
      list(
        mean  = ~mean(.x, na.rm = TRUE),
        sd    = ~sd(.x,   na.rm = TRUE),
        ci_lo = ~quantile(.x, probs = (1 - CI) / 2,       na.rm = TRUE),
        ci_hi = ~quantile(.x, probs = 1 - (1 - CI) / 2,  na.rm = TRUE)
      ),
      .names = "{.fn}_boot_{.col}"
    )
  ) %>%
  mutate(season = year)

  }
  
  growing_seasons %>% 
    pull({{ ndvi_obs }}) -> y_vals

  smooth.spline(growing_seasons$im_doy, y_vals, df = 7) -> annual_spline

  seq(min(growing_seasons$im_doy), max(growing_seasons$im_doy), by = 1) -> x_daily

  x_daily[x_daily <= 180] -> x_daily_spring
  x_daily[x_daily > 180] -> x_daily_fall

  diff(predict(annual_spline, x_daily_spring)$y) -> dy_spring
  diff(predict(annual_spline, x_daily_fall)$y) -> dy_fall

  x_daily_spring[which.max(dy_spring)] -> sos_doy_est
  x_daily_fall[which.min(dy_fall)] -> eos_doy_est

  x0 <- growing_seasons$im_doy
  y0 <- y_vals

  m1_idx <- which(x0 < 90 | x0 > 330)
  m1 <- quantile(y0[m1_idx], 0.50)
  m2 <- quantile(y0, 0.99) - m1
  m3 <- sos_doy_est
  m4 <- 10
  m5 <- eos_doy_est
  m6 <- 10
  m7 <- 4
  
  fit1 <- nlsLM(y0 ~ a + (b - g*x0) * ((1/(1+exp((c-x0)/d))) - (1/(1+exp((e-x0)/f)))),
                start = list(a=m1, b=m2, c=m3, d=m4, e=m5, f=m6, g=m7))

  p1 <- fit1$m$getPars()
  w1 <- max(abs(residuals(fit1)))^2 - abs(residuals(fit1))^2

  fit2 <- nlsLM(y0 ~ a + (b - g*x0) * ((1/(1+exp((c-x0)/d))) - (1/(1+exp((e-x0)/f)))),
                start = list(a=p1[1], b=p1[2], c=p1[3], d=p1[4], e=p1[5], f=p1[6], g=p1[7]),
                weights = w1)

  p2 <- fit2$m$getPars()
  w2 <- max(abs(residuals(fit2)))^2 - abs(residuals(fit2))^2

  fit3 <- nlsLM(y0 ~ a + (b - g*x0) * ((1/(1+exp((c-x0)/d))) - (1/(1+exp((e-x0)/f)))),
                start = list(a=p2[1], b=p2[2], c=p2[3], d=p2[4], e=p2[5], f=p2[6], g=p2[7]),
                weights = w2)

  p3 <- fit3$m$getPars()
  w3 <- max(abs(residuals(fit3)))^2 - abs(residuals(fit3))^2

  fitted_vals <- elmore_curve(p3, x0)
  resids <- y0 - fitted_vals

  pred_mat <- matrix(NA, nrow = max(x0) - min(x0), ncol = boot_int)

  doy <- seq(min(x0), max(x0))

  pred_mat <- matrix(NA, nrow = length(doy), ncol = boot_int)

  for (j in seq(boot_int)) {
    yb <- fitted_vals + sample(resids, replace = TRUE)
    try({
      fit_b <- nlsLM(yb ~ a + (b - g*x0) * ((1/(1+exp((c-x0)/d))) - (1/(1+exp((e-x0)/f)))),
                     start = list(a=p3[1], b=p3[2], c=p3[3], d=p3[4], e=p3[5], f=p3[6], g=p3[7]),
                     weights = w3)
      pred_mat[, j] <- elmore_curve(fit_b$m$getPars(), doy)
    }, silent = TRUE)
  }

  fitted_curve <- elmore_curve(p3, doy)

    d1 <- diff(fitted_curve)
    d2 <- diff(d1)
    d1_aligned <- d1[-length(d1)]
    K  <- d2 / (1 + d1_aligned^2)^(3/2)
    Kprime <- diff(K)
    Kpp    <- diff(Kprime)
    doy_Kp <- doy[1:length(Kprime)]

    # Local extrema of K' are where K'' changes sign
    sc  <- which(diff(sign(Kpp)) != 0)
    directions <- diff(sign(Kpp))[sc]

    # Local maxima of K' (K'' goes + to -) → Greenup, Maturity
    maxima_doys <- doy_Kp[sc[directions < 0]]
    # Local minima of K' (K'' goes - to +) → Senescence, Dormancy
    minima_doys <- doy_Kp[sc[directions > 0]]

    doy_d2           <- doy[1:length(d2)]
    concavity_doys   <- doy_d2[which(diff(sign(d2)) != 0)]
    
    concavity_doys <- len_zero(concavity_doys)
    maxima_doys <- len_zero(maxima_doys)
    minima_doys <- len_zero(minima_doys)
    
    peak_doy <- doy[which.max(fitted_curve)]

    concavity_spring <- len_zero(concavity_doys[concavity_doys < peak_doy])
    concavity_fall   <- len_zero(concavity_doys[concavity_doys > peak_doy])

    maxima_spring <- len_zero(maxima_doys[maxima_doys < peak_doy])
    minima_fall   <- len_zero(minima_doys[minima_doys > peak_doy])
    
   tryCatch(get_transitions_doys(pred_mat, doy, year = unique(growing_seasons$season)),
    error = function(e) {
      warning("Bootstrap transition extraction failed for UID ", 
              unique(growing_seasons$UID), ": ", e$message)
      NULL}) -> boot_transitions 
    
      data.frame(
        UID = unique(growing_seasons$UID),
        flag = unique(growing_seasons$flag),
        season = unique(growing_seasons$season),
        meth = "Elmore",
        R2 = 1 - (sum((y0 - elmore_curve(p3, x0))^2 / sum((y0 - mean(y0))^2))),
        R = sqrt(1 - (sum((y0 - elmore_curve(p3, x0))^2 / sum((y0 - mean(y0))^2)))),
        RMSE = sqrt(mean((y0 - elmore_curve(p3, x0))^2)),
        Greenup_doy = safe_get(maxima_spring, 2, 1),
        SOS_doy = safe_get(concavity_spring, 1, 1),
        Maturity_doy = safe_get(maxima_spring, 2, 2),
        Senescence_doy = safe_get(minima_fall, 2, 1),
        EOS_doy = safe_get(concavity_fall, 1, 1),
        Dormancy_doy = safe_get(minima_fall, 2, 2),
        Peak_NDVI_doy = peak_doy) -> output_df
      
        if (!is.null(boot_transitions)) {
          output_df %>% 
            left_join(boot_transitions, by = "season")
          } else {
            output_df 
            }
}

# fit_annual_curve(seasons[[3]], ndvi_obs = ndvi_normalized) 

