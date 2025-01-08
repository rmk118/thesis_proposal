

mytheme <- theme_classic() + #define custom theme for ggplots
  theme(
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    text = element_text(size = 13))

mytheme_light <- theme_light() + #define custom theme for ggplots
  theme(
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    text = element_text(size = 13))

# Fake_crustaceans --------------------------------------------------------

fake_crustaceans <- function(L50 = 100, # length at 50% maturity on ref var scale
                             slope = 5, # slope parameter for logistic maturity
                             n = 1000, # number of crustaceans sampled
                             # mean of reference variable, e.g., carapace width in mm
                             x_mean = 105,
                             # standard deviation of reference variable
                             x_sd = 20,
                             allo_params = c(1.2, # immature slope parameter
                                             0.1, # immature intercept parameter
                                             1.2, # mature slope parameter
                                             0.1),# mature intercept parameter
                             error_scale = 20) # SD of errors
{
  
  
  # Create normal distribution of carapace widths for a given n, mean, and SD
  fake_crustaceans <- data.frame(x = stats::rnorm(n = n, mean = x_mean, sd = x_sd))
  
  # Add probability of maturity for each individual crab
  # based on a logistic distribution with given location (L50) and
  # shape (slope of the logistic curve) parameters
  fake_crustaceans$prob_mat <- stats::plogis(fake_crustaceans$x, L50, slope)
  
  # Based on the probabilities of maturity,
  # use a binomial distribution to assign each crab a maturity status
  # (0 = immature, 1 = mature)
  mature_vec <- stats::rbinom(n, 1, fake_crustaceans$prob_mat)
  
  # Add vector of maturities to data frame of x-vars and maturity probabilities
  fake_crustaceans$mature <- as.factor(mature_vec)
  
  err_sd <- fake_crustaceans %>%
    dplyr::summarise(
      range = max(.data$x, na.rm = TRUE) - min(.data$x, na.rm = TRUE)
    ) %>%
    dplyr::mutate(err_sd = .data$range * 0.01 / error_scale) %>%
    dplyr::pull(err_sd)
  
  err <- stats::rnorm(n = n, sd = err_sd)
  fake_crustaceans$errs <- exp(err)
  
  a0 <- allo_params[1] # Immature slope parameter
  b0 <- allo_params[2] # Immature intercept parameter
  a1 <- allo_params[3] # Mature slope parameter
  b1 <- allo_params[4] # Immature intercept parameter
  
  fake_crustaceans <- fake_crustaceans %>%
    #if crab is immature, use immature parameters
    dplyr::mutate(y = dplyr::case_when(
      .data$mature == 0 ~ b0 * (.data$x ^ (a0)) * .data$errs,
      #if crab is mature, use mature parameters
      .data$mature == 1 ~ b1 * (.data$x ^ (a1)) * .data$errs),
      log_x = log(.data$x),  #find log of x
      log_y = log(.data$y)  #find log of x
    )
  
  fake_crustaceans <- fake_crustaceans %>% dplyr::select(-"errs")
  
  return(fake_crustaceans)
  
}



# Broken-stick Stevens ----------------------------------------------------

broken_stick_stevens <- function(dat,
                                 xvar,
                                 yvar,
                                 lower = NULL,
                                 upper = NULL,
                                 verbose = FALSE) {
  
  stevens <- dat %>% dplyr::arrange(.data[[xvar]])
  
  xraw <- stevens[[xvar]]
  yraw <- stevens[[yvar]]
  
  if (is.null(lower)) {
    lower <- stats::quantile(xraw, 0.2)
  }
  
  if (is.null(upper)) {
    upper <- stats::quantile(xraw, 0.8)
  }
  
  left_x <- (xraw <= lower) # T/F vector
  low_ndx <- sum(left_x) # largest group 1 point
  right_x <- (xraw >= upper) # T/F vector
  high_ndx <- (length(xraw) - sum(right_x)) + 1 # smallest group 2 point
  min_x <- xraw[low_ndx] # lowest T value
  min_y <- yraw[low_ndx] # lowest T value
  
  stevens$xvar <- xraw
  stevens$yvar <- yraw
  
  lm0 <- stats::lm(yvar ~ xvar, data = stevens)
  rss0 <- stats::anova(lm0)[[2, 2]] # residual sum of squares
  ms0 <- stats::anova(lm0)[[3]] # mean squared error
  F0 <- ms0[1] / ms0[2] # F value
  n0 <- dim(stevens)[1]
  rss_min <- rss0
  mse0 <- mean(lm0$residuals ^ 2)
  
  # assign group membership
  # 1 = left line, 2= right line
  memb <- rep(1, nrow(stevens))
  memb_low <- (xraw <= min_x) # T/F list if less than low range
  memb_high <- (yraw > min_y) # T/F list if GT than high range
  memb[memb_low] <- 1 # assign 1 to those < low
  memb[memb_high] <- 2 # assign 2 to those > high
  memb_sum1 <- summary(as.factor(memb))
  stevens$prior <- memb
  stevens$group <- memb
  
  run <- 0
  
  while (min_x < upper) {
    run <- run + 1
    
    # Left regression
    lm1 <- stats::lm(
      I(yvar[memb == 1] - min_y) ~ 0 + I(xvar[memb == 1] - min_x),
      data = stevens
    )
    b1 <- stats::coef(lm1)[[1]]
    a1 <- min_y - (b1 * min_x)
    df1 <- stats::anova(lm1)[[1]]
    rss1 <- stats::anova(lm1)[[2, 2]]
    ms1 <- stats::anova(lm1)[[3]]
    
    # Right regression
    lm2 <- stats::lm(
      I(yvar[memb == 2] - min_y) ~ 0 + I(xvar[memb == 2] - min_x),
      data = stevens
    )
    b2 <- stats::coef(lm2)[[1]]
    a2 <- min_y - (b2 * min_x)
    df2 <- stats::anova(lm2)[[1]]
    rss2 <- stats::anova(lm2)[[2, 2]]
    ms2 <- stats::anova(lm2)[[3]]
    
    # calculate combined RSS and F
    rss_pool <- rss1 + rss2 # add residual sum of squares
    ms_diff <- (rss0 - rss_pool) / 2
    ms_pool <- rss_pool / (n0 - 4)
    F2 <- ms_diff / ms_pool
    F2_p <- 1 - stats::pf(F2,
                          df1 = 2,
                          df = n0 - 4,
                          lower.tail = F)
    
    if (run == 1 |
        (rss_pool < rss_min)) {
      # Run 1 OR pooled RSS
      rss_min <- rss_pool
      joint_x <- min_x
      joint_y <- min_y
      a1_1 <- a1 # reset old values
      a2_1 <- a2
      b1_1 <- b1
      b2_1 <- b2
    }
    
    # next point
    low_ndx <- low_ndx + 1
    min_x <- stevens$xvar[low_ndx]
    min_y <- stevens$yvar[low_ndx]
    memb_low <- stevens$xvar <= min_x # T/F list if less than low range
    memb_high <- stevens$xvar > min_x # T/F list if GT than high range
    memb[memb_low] <- 1 # assign 1 to those < low
    memb[memb_high] <- 2 # assign 2 to those > high
  } # end loop
  
  SM50 <- joint_x
  
  memb_low <- stevens$xvar <= joint_x # T/F list if less than low range
  memb_high <- stevens$xvar > joint_x # T/F list if GT than high range
  memb[memb_low] <- 1 # assign 1 to those < low
  memb[memb_high] <- 2 # assign 2 to those > high
  stevens$group <- memb
  memb_sum2 <- summary(as.factor(stevens$group))
  n_tot <- sum(memb_sum2)
  
  output <- list(
    data = stevens %>% dplyr::select(-c("xvar", "yvar", "prior")),
    SM50 = SM50,
    imm_slope = b1_1,
    imm_int = a1_1,
    mat_slope = b2_1,
    mat_int = a2_1,
    F_val = F2,
    p_val = F2_p
  )
  if (verbose == TRUE) {
    return(output)
  }
  else
    return(SM50)
  
  
}


# Two-line Stevens --------------------------------------------------------

two_line_stevens <- function(dat,
                             xvar,
                             yvar,
                             lower = NULL,
                             upper = NULL,
                             verbose = FALSE,
                             bps = "even",
                             num_bps = 100) {
  stevens <- dat %>% dplyr::arrange(.data[[xvar]])
  
  xraw <- stevens[[xvar]]
  yraw <- stevens[[yvar]]
  
  if (is.null(lower)) {
    lower <- stats::quantile(xraw, 0.2)
  }
  
  if (is.null(upper)) {
    upper <- stats::quantile(xraw, 0.8)
  }
  
  left_x <- (xraw <= lower) # T/F vector
  low_ndx <- sum(left_x) # largest group 1 point
  right_x <- (xraw >= upper) # T/F vector
  high_ndx <- (length(xraw) - sum(right_x)) + 1 # smallest group 2 point
  min_x <- xraw[low_ndx] # lowest T value
  min_y <- yraw[low_ndx] # lowest T value
  
  stevens$xvar <- xraw
  stevens$yvar <- yraw
  
  lm0 <- stats::lm(yvar ~ xvar, data = stevens)
  rss0 <- stats::anova(lm0)[[2, 2]] # residual sum of squares
  ms0 <- stats::anova(lm0)[[3]] # mean squared error
  F0 <- ms0[1] / ms0[2] # F value
  n0 <- dim(stevens)[1]
  rss_min <- rss0
  mse0 <- mean(lm0$residuals ^ 2)
  
  # assign group membership
  # 1 = left line, 2= right line
  memb <- rep(1, nrow(stevens))
  memb_low <- (xraw <= min_x) # T/F list if less than low range
  memb_high <- (yraw > min_y) # T/F list if GT than high range
  memb[memb_low] <- 1 # assign 1 to those < low
  memb[memb_high] <- 2 # assign 2 to those > high
  memb_sum1 <- summary(as.factor(memb))
  stevens$group <- memb
  
  #### Loop
  
  if (bps == "obs") {
    mse <- rep(0, n0)
    
    for (i in 1:n0) {
      piecewise1 <- stats::lm(
        yvar ~ xvar * (xvar < xvar[i]) + xvar * (xvar >= xvar[i]),
        data = stevens)
      mse[i] <- mean(piecewise1$residuals ^ 2)
    }
    
    ### find breakpoint (bp) that gives lowest MSE
    bp_ind <- which(mse == min(mse))
    bp <- stevens$xvar[bp_ind] # this is not necessarily where the lines cross
    
  }
  
  if (bps == "even") {
    ## determine increment for loop
    steps <- seq(lower, upper, l = num_bps)
    
    #### Loop
    mse <- rep(0, num_bps)
    for (i in 1:num_bps) {
      piecewise1 <- stats::lm(yvar ~ xvar * (xvar < steps[i]) +
                                xvar * (xvar >= steps[i]), data = stevens)
      mse[i] <- mean(piecewise1$residuals ^ 2)
    }
    
    ### find breakpoint (bp) that gives lowest MSE
    bp_ind <- which(mse == min(mse))
    bp <- steps[bp_ind] # this is not necessarily where the lines cross
  }
  
  if (length(bp) > 1) {
    bp <- stats::median(bp)
  }
  
  ## rerun piecewise regression at best bp
  piecewise2 <- stats::lm(yvar ~ xvar * (xvar < bp) + xvar * (xvar > bp),
                          data = stevens)
  
  pw_vals <- stats::coef(piecewise2)
  pw_vals[which(is.na(pw_vals))] <- 0
  a_lo <- pw_vals[1] + pw_vals[3]
  b_lo <- pw_vals[2] + pw_vals[5]
  a_hi <- pw_vals[1] + pw_vals[4]
  b_hi <- pw_vals[2]
  
  jx <- as.numeric((a_lo - a_hi) / (b_hi - b_lo)) #the point where 2 lines meet
  
  ####  Reassign group membership
  memb_pw <- rep(1, n0)
  memb_pw[stevens$xvar >= bp] <- 2
  stevens$group <- memb_pw
  
  output <- list(
    data = stevens,
    breakpoint = bp,
    intersection = jx,
    imm_slope = b_lo,
    imm_int = a_lo,
    mat_slope = b_hi,
    mat_int = a_hi
  )
  
  if (verbose == TRUE) {
    return(output)
  }
  else
    return(c(breakpoint = bp, intersection = jx))
  
}



# REGRANS -----------------------------------------------------------------

regrans_fun <- function(dat,
                        xvar,
                        yvar,
                        lower = NULL,
                        upper = NULL,
                        verbose = FALSE,
                        n_tries = 100) {
  
  x <- dat[[xvar]]
  y <- dat[[yvar]]
  
  if (is.null(lower)) {
    lower <- stats::quantile(x, 0.2)
  }
  
  if (is.null(upper)) {
    upper <- stats::quantile(x, 0.8)
  }
  
  changept_choices <- seq(lower, upper, length.out = n_tries)
  
  help_fun <- function(i)
  {
    x2star <- (x - i) * as.numeric(x > i)
    fit <- stats::lm(y ~ x + x2star)
    sum_sq <- stats::anova(fit)["Residuals", "Sum Sq"]
    return(sum_sq)
  }
  
  breakpt <- sapply(changept_choices, help_fun)
  
  breakpt <- data.frame(changept = changept_choices, sum_sq = breakpt)
  
  if (verbose == TRUE) {
    return(breakpt)
  }
  else {
    out <- breakpt[which.min(breakpt$sum_sq), "changept"]
    return(out)
  }
  
}


# Somerton method ---------------------------------------------------------

somerton_fun <- function(
    dat, # data.frame with columns corresponding to xvar, yvar
    xvar, # X variable
    yvar, # Y variable
    trans = "none", # transformation to apply
    lower = NULL, # lower bound of unknown range
    upper = NULL, # upper bound of unknown range
    max_iter = 50 # maximum number of iterations
) {
  
  if (is.null(lower)) {
    lower <- stats::quantile(dat[[xvar]], 0.2)
  }
  
  if (is.null(upper)) {
    upper <- stats::quantile(dat[[xvar]], 0.8)
  }
  
  if (trans == "log") {
    dat$xvar <- log(dat[[xvar]])
    dat$yvar <- log(dat[[yvar]])
  }
  else if (trans == "std") {
    dat$xvar <- scale(dat[[xvar]])
    dat$xvar <- scale(dat[[yvar]])
  }
  else {
    dat$xvar <- dat[[xvar]]
    dat$yvar <- dat[[yvar]]
  }
  
  
  df <- dat %>%
    mutate(group = case_when(
      xvar < lower ~ "juv",
      xvar > upper ~ "adult",
      .default = NA))
  
  df$temp_group <- df$group
  
  
  rsq_vec <- rep(NA, max_iter)
  RSS_vec <- rep(NA, max_iter)
  
  for (i in 1:max_iter) {
    
    # fit known juveniles
    juv_df <- df[df$temp_group=="juv",]
    juv_fit <- lm(yvar ~ xvar, juv_df) # fit linear model
    rss_juv <- anova(juv_fit)[[2]][2] # juvenile model residual sum of squares
    
    # fit known adults
    ad_df <- df[df$temp_group=="adult",]
    ad_fit <- lm(yvar ~ xvar, ad_df) # fit linear model
    rss_ad <- anova(ad_fit)[[2]][2] # adult model residual sum of squares
    
    RSS <- rss_juv + rss_ad # add residual sum of squares
    TSS <- sum((df$yvar - mean(df$yvar))^2) # total sum of squares
    rsq <- 1-(RSS/TSS)
    
    df <- df %>%
      mutate(
        pred_juv = predict(juv_fit, newdata = data.frame(xvar)), # predict yvar based on juvenile model
        pred_ad = predict(ad_fit, newdata = data.frame(xvar)), # predict yvar based on adult model
        resid_juv = abs(yvar - pred_juv), # juvenile residuals
        resid_ad = abs(yvar - pred_ad), # adult residuals
        # temp_group = if_else(resid_juv < resid_ad, "juv", "adult"))
        # Option 1: all points can be reclassified to either maturity stage
        temp_group = case_when(
          is.na(group) & resid_juv < resid_ad ~ "juv",
          is.na(group) & resid_juv >= resid_ad ~ "adult",
          .default = group
        ))
    
    rsq_vec[i] = rsq
    RSS_vec[i] = RSS
  }
  
  df <- df %>% rename(init_group = group, pred_mat =  temp_group) %>%
    mutate(pred_mat_num = if_else(pred_mat == "adult", 1, 0))
  
  output <- list(data = df, rsq = rsq_vec, RSS = RSS_vec, juv_mod = juv_fit, adult_mod = ad_fit)
  return(output)
}


# infl_pt -----------------------------------------------------------------

infl_pt <- function(dat, x, y, log = FALSE, plot = FALSE) {
  # find the ratio between the two morphometric variables
  if (isTRUE(log)) {
    ratio <- log(dat[[y]])/log(dat[[x]])
  }
  else {
    ratio <- dat[[y]]/dat[[x]]
  }
  
  
  # compute a kernel density estimate (essentially a smoothed histogram)
  density_test <- stats::density(ratio)
  
  # convert into a data frame
  density_test <- data.frame(x = density_test$x, y = density_test$y)
  
  # find the local minimum between the two peaks
  density_test$is_min <- splus2R::peaks(
    x = -density_test$y, span = 3, strict = FALSE)
  
  min <- density_test %>%
    dplyr::filter(.data$is_min == TRUE) %>%
    dplyr::pull(x)
  
  min <- stats::median(min)
  
  # optionally visualize the density plot with minimum
  if (plot == TRUE) {
    print(ggplot2::ggplot() +
            ggplot2::geom_line(aes(x = density_test$x, y = density_test$y)) +
            ggplot2::geom_vline(xintercept = min, lty = "dashed") +
            labs(x = "Ratio", y = NULL) +
            ggplot2::theme_light())
  }
  
  # return the minimum ratio, equivalent to the slope of a line
  # separating the two clouds of points
  if (is.na(min)) {
    warning("No local minimum detected.")
  }
  
  return(min)
}
