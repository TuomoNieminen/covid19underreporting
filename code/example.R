# Example on how to estimate seroprevalence and underrerporting as described in the article
# This example uses week-level data whereas the article was based on daily data
# The daily data is not public due to privacy considerations
# This example reproduces the results of Table 2 APROXIMATELY, not exactly
# Tuomo A. Nieminen 2022

# requires installation of r package rstan and STAN
library(rstan)

###########################
## Data and hyperparameters
# -------------------------

# Data from Table 2
df <- read.table("data/table2.csv", sep = "\t")
df$wednesday <- as.Date(df$StartDate) + 2
last_date <- as.Date(max(df$EndDate))

## Projection model data (Registered infections and time lag to antibodies)
# -----------

# COVID-19 weekly counts from cumulative counts
df$Study <- ifelse(is.na(df$Study), 0, df$Study)
df$StudyInfectionsWeek <- c(0, diff(df$Study))

# as the < 10 obs were not shown, assign the cumulative count from week 4 evenly to first 4 weeks
df$StudyInfectionsWeek[1:4] <- round(df$StudyInfectionsWeek[4]/4)

# compute individual times since infection
# assume all infections occurred on wednesdays
t_since_infection <- as.numeric(last_date - df$wednesday)
indv_t_since_infection <- rep(t_since_infection, times = df$StudyInfectionsWeek)


# Tan et al 2020 data
# ----

data.time2antibodies_distribution <- function() {
  
  data <-  list(
    t = c(7,10,14,21,28,35,42,49),
    n = c(58,62,61,54,35,22,15,5),
    y = c(2,12,31,32,26,17,13,4)
  )
  names(data) <- paste0("immunity_", names(data))
  
  data$immunity_N <- length(data$immunity_t)
  data
}

immu_data <- data.time2antibodies_distribution()

projection_model_data <- c(list(t = indv_t_since_infection, 
                                n_ttr = length(indv_t_since_infection), 
                                popsize = 1.00e6), 
                           immu_data)


## Estimation model data (Serosurvey and hyperparameters)
# ------

start_date <- as.Date("2020-04-09")
dates <- seq.Date(from=start_date, to = last_date, by = "day")
D <- length(dates) # 87 days
df_sero <- df[df$wednesday >= start_date, ]

n <- df_sero$Samples
y <- df_sero$Confirmationpos

# need to expand the weekly data into daily data for the stan code
# We simply input 6 days of no samples (sample size 0) for each week
n <- array(0,  D)
y <- array(0,  D)

# place serosurvey sample data to wednesdays 
y[(0:(nrow(df_sero)-1))*7 + 3] <- df_sero$Confirmationpos
n[(0:(nrow(df_sero)-1))*7 + 3] <- df_sero$Samples


# estimation model hyperparameters
# ----
spec <- 1
sigma_t_alpha <- 2
sigma_t_beta <- 40
prior_sigma <- 2 # logit scale
prior_mu <- 0.05 # probability scale
 

estimation_model_data <- list(y = y, 
                              n = n,
                              N = D,
                              prior_mu = prior_mu, 
                              prior_sigma = prior_sigma,
                              sigma_t_alpha = sigma_t_alpha,
                              sigma_t_beta = sigma_t_beta,
                              spec = spec)

# all input data combined
standata <- c(projection_model_data, 
              estimation_model_data, 
              list(SMOOTHING_WINDOW = 7,
                   positive_trend = 1)) 


###########
# Estimation
###########

# some sensible initial values for parameters

logit <- function(x) log(x / (1-x))

init <- function(chain_id) {
  
  init_sigma_t <- sigma_t_alpha / sigma_t_beta
  init_p_logit_inf <- rnorm(D, logit(sum(y)/sum(n)), init_sigma_t)
  
  list(sigma_t = init_sigma_t, 
       p_logit_inf_w = init_p_logit_inf,
       meanlog = log(14),  # in Tan et al. 50% had serovonverted by day 14 so use that as initial guess
       sdlog = 1)
}

# use multiple cores when sampling
ncore <- max(1, floor(parallel::detectCores() / 2))
options(mc.cores = ncore)

# read the stan model, may give a parser diagnostics message for integer division
model <- rstan::stan_model(file = file.path("code/model.stan"))

# do sampling
# may give 'Bayesian Fraction of Missing Information was low' warning.
fit <- rstan::sampling(model,
                       data = standata,
                       pars = c("p_sero_w", "p_test_w", "p_sero", "meanlog","sigma_t", "sdlog", "ttr_seroconverted", "ttr_p_sero", "unobserved_ratio"),
                       control  = list(adapt_delta = 0.98, max_treedepth = 20),
                       iter = 30000,
                       init = init)


# optinally, save the model
# saveRDS(fit, file = "data/fitted_model.Rds")


##########
# Results
##########

# Reproduce Table 3 (excluding screening test results)
# ----

# read the saved model
# fit <- readRDS("data/fitted_model.Rds")

# extract samples and compute 95% credible intervals
# ----
sero_projection <- summary(fit, 
             pars = c("ttr_p_sero"),
             probs = c(0.025, 0.975))$summary
projres <- as.data.frame((round(100*sero_projection,2)))
projres  <- with(projres, paste0(mean, " (",`2.5%`,"-", `97.5%`, ")"))

sero_estimate <- summary(fit, 
                           pars = c("p_sero"),
                           probs = c(0.025, 0.975))$summary
serores <- as.data.frame((round(100*sero_estimate,2)))
serores  <- with(serores, paste0(mean, " (",`2.5%`,"-", `97.5%`, ")"))

underreporting <- summary(fit, 
                           pars = c("unobserved_ratio"),
                           probs = c(0.025, 0.975))$summary
uderres <- as.data.frame((round(underreporting,2)))
uderres  <- with(uderres, paste0(mean, " (",`2.5%`,"-", `97.5%`, ")"))


# Aproximate reproduction of Table 3 main results, including:
# - COVID-19 -based seroprevalence projection
# - Confirmation test -based seroprevalence estimate
# - Confirmation test -based underreporting estimate
results_df <- cbind.data.frame(projres, serores, uderres)
results_df$Date <- seq.Date(from=start_date, to = last_date, by = "day")
results_df <- results_df[results_df$Date %in% df_sero$wednesday, ]
results_df <- results_df[, c(4, 1:3)]
colnames(results_df)  <- c("date", "Proj.Prevalence", "Est.Prevalence", "Underreporting")
View(results_df)
