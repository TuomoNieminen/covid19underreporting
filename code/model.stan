// Bayesian model for estimating underreporting of SARS-CoV-2 infections 
// during the first epidemic wave in 2020
// Implementation in STAN
// Tuomo Nieminen 2021
// tuomo.nieminen (at) thl.fi


data {
  
  // sero sample data for t = 1, .., N time points
  int<lower=1> N; // number of time points (days). This is D + 1 in the article.
  
    // length in days for which the seroprevalence is constant (e.g. 7)
  int<lower=1> SMOOTHING_WINDOW; 
  
    // indicator for using a positive trend in seroprevalence estimation
  int<lower=0, upper=1> positive_trend; 

  int<lower=0> n[N]; // number of serology samples for each day
  int<lower=0> y[N]; // number of positive samples for each day
  
    // hyperparameters for prevalence at t=1
  real<lower=0,upper=1> prior_mu; // seroprevalence mean, probability scale
  real<lower=0> prior_sigma; // seroprevalence stardard deviation, logit scale
  
    // strength of dependency in seroprevalence over time
    // hyperparameters for prior Gamma distribution of sigma
  real<lower=0>sigma_t_alpha;
  real<lower=0>sigma_t_beta;
  
    // antibody test specificity
  real<lower=0, upper=1> spec;
  
  
    // data to estimate the time to seroconversion distribution (Tan et al.)
  int<lower=1>immunity_N; // number of time points (8 in Tan et al. data)
  real<lower=0> immunity_t[immunity_N]; // days since symptom onset (Q_j)
  int immunity_n[immunity_N]; // number of samples tested at each day
  int immunity_y[immunity_N]; // number of positive samples at each day
  
  // FNIDR data
  int<lower=1> n_ttr; // number of COVID-19 cases
  real<lower=0> t[n_ttr]; // for each case, time since symptom onset to the 
                          // last day of interest, i.e. t_i = N - r_i - C
  real<lower=1> popsize;  // population size
  
}


transformed data {

  int PERIODS = (N - 1) / SMOOTHING_WINDOW + 1;  // number of periods (weeks)
  int<lower=0> n_period[PERIODS]; // number of serosamples during periods
  int<lower=0> y_period[PERIODS]; // number of positive samples during periods

    // helper index, used later to expand weekly seroprevalence for each day
  int<lower=1> ii[N]; 

  for(w in 1:PERIODS) {
    n_period[w] = 0;
    y_period[w] = 0;
  }
  
  for(i in 1:N) {
    int period = (i - 1) / SMOOTHING_WINDOW + 1;
    n_period[period] = n_period[period] + n[i];
    y_period[period] = y_period[period] + y[i];
    ii[i] = period; // values of ii are period indicators.
  }
  
}


parameters {
  
  real p_logit_sero_w[PERIODS]; // logit seroprevalence at w
  real<lower=0> sigma_t; //strength of dependency 
    
  // time to antibodies
  real meanlog;
  real<lower=0> sdlog;
  
}

transformed parameters {
  
  // SARS-CoV-2 seroprevalence
  real<lower=0, upper=1> p_sero_w[PERIODS] = inv_logit(p_logit_sero_w);
  
  // seropositivity prevalence
  real<lower=0, upper=1> p_test_w[PERIODS];

    //  seropositivity prevalence is a function of 
    //  seroprevalence using test specificity
  for(i in 1:PERIODS) {
    p_test_w[i] = p_sero_w[i] + (1 - p_sero_w[i])*(1-spec);
  }
  
}


model {
  
  // seroprevalence strength of dependency 
  sigma_t ~ gamma(sigma_t_alpha, sigma_t_beta);

  // prior for (logit) seroprevalence at t = 1
  p_logit_sero_w[1] ~ normal(logit(prior_mu), prior_sigma);

  // prior for (logit) seroprevalence at t = 2
  p_logit_sero_w[2] ~ normal(p_logit_sero_w[1], sigma_t);

  // prior for (logit) seroprevalence at t > 2
  for(w in 3:PERIODS) {
    real trend;
    
    if(positive_trend == 1)
        trend = max([0, p_logit_sero_w[w-1] - p_logit_sero_w[w-2]]);

    if(positive_trend == 0)
        trend = 0;
      
    p_logit_sero_w[w] ~ normal(p_logit_sero_w[w-1] + trend, sigma_t);
    
  }

  // likelihood for seropositive observations
  y_period ~ binomial(n_period, p_test_w);
  
  // likelihood for time to seroconversion
  for(i in 1:immunity_N) {
    real p_seroconversion;
    p_seroconversion = lognormal_cdf(immunity_t[i], meanlog, sdlog);
    immunity_y[i] ~ binomial(immunity_n[i], p_seroconversion);
  }
  
}

generated quantities {
  
  // FNIDR seroprevalence projection and corresponding ratio estimate
  
  real<lower=0> ttr_seroconverted[N]; // projected number seroconverted for each day
  real<lower=0, upper=1> ttr_p_sero[N]; // projected seroprevalence for each day
  real<lower=0> unobserved_ratio[N]; // The ratio (Delta in the article)
  real p_sero[N] = p_sero_w[ii]; //  estimated seroprevalence for each day
  
  // for each time index and each COVID-19 case, 
  // predict the probability of seroconversion (projection)
  for(time_index in 1:N) {
    ttr_seroconverted[time_index] = 0;
    
    for(i in 1:n_ttr) {
      real time_since_symptoms = max([0, t[i] - (N - time_index)]);
      real immu_prob = lognormal_cdf(time_since_symptoms, meanlog, sdlog);
      ttr_seroconverted[time_index] += immu_prob;
    }
    
    // seroprevalence projection for the time index
    ttr_p_sero[time_index] = ttr_seroconverted[time_index] / popsize;
    
    // ratio for the time index
    unobserved_ratio[time_index] = p_sero[time_index] / ttr_p_sero[time_index];
  }
  
}
