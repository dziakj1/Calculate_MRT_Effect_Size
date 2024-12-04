# Adapted by John Dziak from Tianchen Qian

calculate_mrt_effect_size <- function(
    data, # long format;
    id,
    outcome,
    treatment,
    time,
    rand_prob = NULL,
    availability = NULL,
    covariates = NULL,
    loess_degree = 1,
    loess_span = 1/4,
    standardize = TRUE,
    do_bootstrap = TRUE,
    code_for_no_prompt = 0,
    code_for_yes_prompt = 1,
    boot_replications = 1000,
    confidence_alpha=.05) {
  
  all_subject_ids <- sort(unique(data[,id]))
  n_participants <- length(all_subject_ids)
  
  all_decision_points <- sort(unique(data[,time]))
  n_decision_points <- length(all_decision_points)
  
  estimate_mrt_effect_size_inner <- function(all_ids,
                                             selected_subjects) {
    
    bt <- numeric(n_decision_points)
    n0 <- numeric(n_decision_points)
    n1 <- numeric(n_decision_points)
    s0 <- numeric(n_decision_points)
    s1 <- numeric(n_decision_points)
    spoolt <- numeric(n_decision_points)
    
    t_axis <- 1:n_decision_points  # Probably should rescale
    
    
    for (decision_point_index in 1:n_decision_points) {
      
      this_t <- all_decision_points[decision_point_index]
      these_ids <- all_subject_ids[selected_subjects]
      
      current_rows <- which((data[,id] %in% these_ids) &
                              (data[,time]==this_t)) 
      
      Y <- as.vector(data[current_rows,outcome])
      A <- as.vector(data[current_rows,treatment])
      if (is.null(covariates)) {
        Z <- NULL
      } else {
        Z <- as.matrix(data[current_rows,
                            covariates,drop=FALSE])
      }
      
      prob <- data[current_rows,rand_prob]
      avail <- data[current_rows,availability]
      
      # Compute weights
      W <- ifelse(A==1, avail/prob, avail/(1-prob))
      
      # Perform weighted least squares regression
      if (is.null(Z)) {
        weighted_ols_fit <- lm(Y ~ A, weights=W)
      } else {
        weighted_ols_fit <- lm(Y ~ A + Z, weights=W)
      }
      
      
      # " At each decision point t: Compute b(t), an estimate of the difference in the
      #  population mean proximal outcome for the two treatment groups. Compute
      # s_pool(t), the pooled standard deviation for the proximal outcome
      
      bt[decision_point_index] <- coef(weighted_ols_fit)["A"]
      n0[decision_point_index] <- length(which(A==code_for_no_prompt))
      n1[decision_point_index] <- length(which(A==code_for_yes_prompt))
      s0[decision_point_index] <- ifelse(length(which(A==code_for_no_prompt)>1),
                                         sd(Y[which(A==code_for_no_prompt)]),
                                         NA)
      s1[decision_point_index] <- ifelse(length(which(A==code_for_yes_prompt)>1),
                                         sd(Y[which(A==code_for_yes_prompt)]),
                                         NA)
      
      # "Even though we use an adjusted mean difference to
      # estimate ^b(t), we use the pooled standard deviation
      # at decision point t to estimate ^sigma(t) rather than
      # the residual standard error. Standardizing by the
      # residual standard error conflates treatment efficacy
      # and the correlation between pretreatment and
      # post-treatment outcome measurements (Feingold,
      # 2009, 2013). "
      
      if (standardize) {
      spoolt[decision_point_index] <- sqrt( ((n0[decision_point_index]-1)*(s0[decision_point_index]^2) +
                                               (n1[decision_point_index]-1)*(s1[decision_point_index]^2)) /
                                              ( n0[decision_point_index]+n1[decision_point_index]-2) )
      } else {
        spoolt[decision_point_index] <- 1  # treat pooled standard deviation as if it was 1 if you want a raw effect size
      }
      }
    
    # "3. Apply a linear smoothing procedure to the pairs
    # {(t1, b(t1)), (t2, b(t2)), . . . , (tT , b(tT ))}.
    # Evaluate the estimated function at
    # each of the points u1, . . . , uM"
    
    
     bt_loess <- predict(loess(bt ~ t_axis,
                               degree=loess_degree,
                               span=loess_span))
     
    # # 4. Apply the same smoothing procedure to the pairs
    # # {(t1, spool(t1)), (t2, spool(t2)), . . . , (tT ,
    # # spool(tT ))}. Evaluate the estimated
    # # function at the same points u1, . . . , uM
    # 
     spoolt_loess <- predict(loess(spoolt ~ t_axis,
                                   degree=loess_degree,
                                   span=loess_span))
     
     
    # # 5. Calculate the ratios  ^d(u) = ˆb(u) / ^s(u)
    # # for u1, . . . , uM.
    # 
     dt <- bt_loess / spoolt_loess
     
    return(dt)
  }
  
  if (do_bootstrap) {
    
    out_boot_1 <- boot(data=all_subject_ids,
                       statistic=estimate_mrt_effect_size_inner,
                       R=boot_replications)
    
    d_t_smoothed <- out_boot_1$t0
    
    ci_d_upper <- numeric(n_decision_points)
    ci_d_lower <- numeric(n_decision_points)
    for (this_t in 1:n_decision_points) {
      out_boot_ci_1 <- boot.ci(out_boot_1, index=this_t, type="perc")
      ci_d_lower[this_t] <- out_boot_ci_1$percent[4]
      ci_d_upper[this_t] <- out_boot_ci_1$percent[5]
    }
    return(data.frame(time = all_decision_points,
                      estimate = d_t_smoothed,
                      lower = ci_d_lower, 
                      upper = ci_d_upper))
  } else {
    fit1 <- estimate_mrt_effect_size_inner(  all_ids = all_subject_ids,
                                             selected_subjects = 1:n_participants)
    return(data.frame(time = all_decision_points,
                      estimate = fit1))
  }
   
} 
