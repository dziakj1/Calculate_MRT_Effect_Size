# Adapted by John Dziak from Tianchen Qian

# This MRT data has the following features that are commonly encountered in practice:
# 1. time-varying randomization probability (prob_treatment) that depend on history
# 2. availability indicator that is sometimes 0

simulate_mrt_data <- function(n_participants,
                              n_decision_points) {
  
  # Function to generate data for a single participant
  simulate_one_participant <- function(id, n_decision_points) {
    # Initialize variables
    decision_point <- 1:n_decision_points
    covariate1 <- numeric(n_decision_points)
    covariate2 <- numeric(n_decision_points)
    prob_treatment <- numeric(n_decision_points)
    treatment <- numeric(n_decision_points)
    availability <- numeric(n_decision_points)
    outcome <- numeric(n_decision_points)
    error <- numeric(n_decision_points) # AR1 error
    
    # AR1 parameters for error
    ar1_rho <- 0.6  # Correlation coefficient for AR1 process
    sigma <- 1      # Standard deviation of error
    
    
    # Auto-regressive coefficient for covariates
    ar_coef <- 0.8
    
    # Generate the time-varying treatment effect curve using a Beta distribution density
    treatment_effect_curve <- dbeta(seq(0, 1, length.out = n_decision_points), shape1 = 2, shape2 = 5)
    treatment_effect_curve <- treatment_effect_curve / max(treatment_effect_curve) * 2  # Scale to max effect = 2
    
    for (t in decision_point) {
      
      # Covariates evolve over time with an AR process and treatment effect
      if (t == 1) {
        # Set initial covariate values
        covariate1[1] <- rnorm(1, mean = 0, sd = 1)
        covariate2[1] <- rnorm(1, mean = 0, sd = 1)
      } else {
        covariate1[t] <- ar_coef * covariate1[t-1] + 0.1 * treatment[t-1] + rnorm(1, mean = 0, sd = 0.5)
        covariate2[t] <- ar_coef * covariate2[t-1] - 0.1 * treatment[t-1] + rnorm(1, mean = 0, sd = 0.5)
      }
      
      # Availability constraint
      availability[t] <- rbinom(1, size = 1, prob = 0.8)
      
      if (availability[t] == 1) {
        # Randomization probability depends on the most recent covariate
        prob_treatment[t] <- plogis(0.5 * covariate1[t] + 0.3 * covariate2[t])
        treatment[t] <- rbinom(1, size = 1, prob = prob_treatment[t])
      } else {
        # If unavailable, no treatment
        treatment[t] <- 0
      }
      
      # Generate AR1 error
      if (t == 1) {
        # Initialize AR1 error
        error[1] <- rnorm(1, mean = 0, sd = sigma)
      } else {
        error[t] <- ar1_rho * error[t-1] + rnorm(1, mean = 0, sd = sqrt(1 - ar1_rho^2) * sigma)
      }
      
      # Outcome is influenced by time-varying treatment effect, covariates, and error
      outcome[t] <- treatment_effect_curve[t] * treatment[t] + 
        0.5 * covariate1[t] + 
        0.3 * covariate2[t] + 
        error[t]
    }
    
    # Return data frame for this participant
    return(data.frame(
      id = id,
      decision_point = decision_point,
      availability = availability,
      prob_treatment = prob_treatment,
      treatment = treatment,
      covariate1 = covariate1,
      covariate2 = covariate2,
      treatment_effect = treatment_effect_curve,  # Add treatment effect for inspection
      sigma = sigma,
      outcome = outcome
    ))
    
  }
  
  mrt_data <- bind_rows(lapply(1:n_participants, 
                               simulate_one_participant, 
                               n_decision_points = n_decision_points))
  mrt_data$id <- 1000 + mrt_data$id  # just to make sure we're not assuming ID's start from 1
  
  return(mrt_data)
}
