# based on draft code by Tianchen (2024.11.13) 
rm(list = ls())

set.seed(16802)

library(tidyverse)
library(stats)
library(boot)
source("simulate_mrt_data.R")
source("calculate_mrt_effect_size.R")

loess_degree <- 1
loess_span <- 1/4

# Some features that are not yet implemented below but would be nice to have:
# (it would be great if John can help with this once the basics of the algorithm
# is implemented)
# - MRT where different individuals have different length of follow-ups
#.   (i.e., individuals have different numbers of decision points)
#  - decision point index may not start at 1 (sometimes start at 0).
#.   (Maybe we just force the software user to have everyone's decision point to 
#.   start at 1?)

# Generate data for all participants
#set.seed(123) # For reproducibility
n_participants <- 100
n_decision_points <- 50

mrt_data <- simulate_mrt_data(n_participants=n_participants,
                              n_decision_points=n_decision_points)
ans1 <- calculate_mrt_effect_size(mrt_data,
                                  id="id",
                                  outcome="outcome",
                                  treatment="treatment",
                                  time="decision_point",
                                  rand_prob="prob_treatment",
                                  availability = "availability",
                                  do_bootstrap = TRUE,
                                  boot_replications = 100,
                                  covariates=c("covariate1","covariate2"))

print(ans1)


