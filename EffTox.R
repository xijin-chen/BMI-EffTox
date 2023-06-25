library(rstan)
library(doParallel)
library(doRNG)
library(purrr)

EffTox_XC <- "
functions {
  real log_joint_pdf(int num_patients, int[] eff_miss,
                    int[] tox, int[] eff, int[] doses, real[] coded_doses,
                    real alpha, real beta, real gamma, real zeta, real eta, real psi) {
    real loglikeli;
    loglikeli = 0;
    for(j in 1:num_patients) {
      real prob_tox;
      real prob_eff;
      real likeli;
      prob_tox = inv_logit(alpha + beta * coded_doses[doses[j]]);
      prob_eff = inv_logit(gamma + zeta * coded_doses[doses[j]] + 
                 eta * (coded_doses[doses[j]])^2);
      if (eff_miss[j] == 1)
        likeli = (prob_tox)^tox[j] * (1 - prob_tox)^(1 - tox[j]);
      else
        likeli = prob_eff^eff[j] * (1 - prob_eff)^(1 - eff[j]) * prob_tox^tox[j] *
              (1 - prob_tox)^(1 - tox[j]) + (-1)^(eff[j] + tox[j]) * prob_eff *
              prob_tox * (1 - prob_eff) * (1 - prob_tox) *
              (exp(psi) - 1) / (exp(psi) + 1);

               
      loglikeli += log(likeli);
    }
    return loglikeli;
  }
}
          
data {
  // mean values matters considering the shape of dose-tox or dose effs curves
  real alpha_mean;
  real<lower=0> alpha_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real gamma_mean;
  real<lower=0> gamma_sd;
  real zeta_mean;
  real<lower=0> zeta_sd;
  real eta_mean;
  real<lower=0> eta_sd;
  real psi_mean;
  real<lower=0> psi_sd;
  int<lower=1> num_doses;
  real<lower=0> real_doses[num_doses];
  
  int <lower=0> num_patients;
  int<lower=1, upper=num_doses> doses[num_patients];
  int<lower=0, upper=1> tox[num_patients];
  int<lower=0, upper=2> eff[num_patients];
  // int eff[num_patients];
  int<lower=0, upper=1> eff_miss[num_patients];
  real Lp;
  real p1e;
  real p2t;
}

transformed data {
  real coded_doses[num_doses];
  real sum_term;
  sum_term = 0;
  for(i in 1:num_doses) {
    sum_term += log(real_doses[i]);
  }
  sum_term = sum_term/num_doses;
  for(i in 1:num_doses) {
    coded_doses[i] = log(real_doses[i]) - sum_term;
  }
}

parameters {
  real alpha;
  real beta;
  real gamma;
  real zeta;
  real eta;
  real psi;
}

transformed parameters {
  real<lower=0, upper=1> prob_tox[num_doses];
  real<lower=0, upper=1> prob_eff[num_doses];
  real utility[num_doses];
  for(i in 1:num_doses) {
    prob_tox[i] = inv_logit(alpha + beta * coded_doses[i]);
    prob_eff[i] = inv_logit(gamma + zeta * coded_doses[i] + eta * coded_doses[i]^2);
    utility[i] = 1-(((1-prob_eff[i])/(1-p1e))^Lp+((prob_tox[i])/(p2t))^Lp)^(1/Lp);
  }
}


model {
  target += normal_lpdf(alpha | alpha_mean, alpha_sd);
  if (1 == 1)
    target += normal_lpdf(beta | beta_mean, beta_sd);
  target += normal_lpdf(gamma | gamma_mean, gamma_sd);
  target += normal_lpdf(zeta | zeta_mean, zeta_sd);
  target += normal_lpdf(eta | eta_mean, eta_sd);
  target += normal_lpdf(psi | psi_mean, psi_sd);
  target += log_joint_pdf(num_patients, eff_miss,
            tox, eff, doses, coded_doses,
            alpha, beta, gamma, zeta, eta, psi);
}
"
model_efftox_xc <- stan_model(model_code = EffTox_XC, verbose = TRUE)


sim_EffTox <- function(decision_t = 20, num_doses = 5, eff_t = 12, tox_t = 4,
                       num_patients_ini = 0,  zeta_mean=2,
                       true_prob_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       true_prob_eff = c(0.00045, 0.00326, 0.03593, 0.89251, 0.98928),
                       real_doses = c(1.0, 2.0, 4.0, 6.6, 10.0),
                       Lp=2.07, p1e=0.5, p2t=0.65,
                       pe_ad=0.1, pt_ad=0.1, pitg=0.3, piel=0.5,
                       tox_ini = integer(length = 0), eff_ini = integer(length = 0), 
                       eff_miss_ini = integer(length = 0), doses_ini = integer(length = 0),
                       start_dose = 1, num_sims=5000, cohort_size = 3, chains=4){
  
  out <- foreach(sim = 1:num_sims, .export = "model_efftox_xc", .packages="rstan") %dorng% {
    # -----------------------------------------------------------
    stop <- decision_t+4
    cohorts <- seq(4, stop, by=4)
    
    eff_outcome = array(NA, dim = c(length(cohorts), length(cohorts), cohort_size))
    eff_miss_outcome = array(NA, dim = c(length(cohorts), length(cohorts), cohort_size))
    tox_outcome = matrix(NA, ncol=cohort_size, nrow=length(cohorts))
    
    eff_outcome_list <- array(NA)
    loc_num = array(NA, dim = c(length(cohorts), stop, cohort_size))
    
    administerdose = NA
    recommenddose = NA
    administerdose[1] <- start_dose 
    arrived_t <- NA
    
    tox = tox_ini
    eff = eff_ini
    eff_miss = eff_miss_ini
    num_patients = num_patients_ini
    doses = doses_ini
    
    efftox_ptox_mean <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    efftox_ptox_se <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    efftox_peff_mean <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    efftox_peff_se <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    efftox_utility_mean <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    efftox_utility_se <- matrix(NA, nrow=length(cohorts), ncol=num_doses)
    
    stop_efficacy <- array(0, dim = length(cohorts))
    stop_toxicity <- array(0, dim = length(cohorts))
    stop_other <- array(0, dim = length(cohorts))
    
    for(i in cohorts){
      cohort_num <- i/4
      arrived_t[cohort_num] <- cohorts[cohort_num]-4
      tox_outcome[cohort_num, 1:cohort_size] <- rbinom(cohort_size, 1, true_prob_tox[administerdose[cohort_num]])
      
      for(patient in 1:cohort_num){
        loc_num[patient, (arrived_t[patient]+1):i, ] <- (1:(i-arrived_t[patient]))
      }
      for(cohort in 1:cohort_size){
        for(j in 1:cohort_num){
          if(loc_num[j, cohort_num*4, cohort] < (eff_t-tox_t+4)){
            eff_outcome[j, cohort_num, cohort] <- 2
            eff_miss_outcome[j,cohort_num,cohort] <- 1
          }else{
            if(loc_num[j,cohort_num*4, cohort] == (eff_t-tox_t+4)){
              eff_outcome[j, cohort_num, cohort] <- rbinom(1, 1, true_prob_eff[administerdose[j]]) 
            }else{
              # 
              eff_outcome[j,cohort_num,cohort] <- eff_outcome[j,cohort_num-1,cohort]
            }
            eff_miss_outcome[j,cohort_num,cohort] <- 0
          }
        }
      }  
      
      num_patients = (num_patients + cohort_size)
      tox = array(t(tox_outcome[1:cohort_num,1:cohort_size]))
      eff = array(t(eff_outcome[1:cohort_num, cohort_num,(1:cohort_size)]))
      eff_miss = array(t(eff_miss_outcome[1:cohort_num, cohort_num,(1:cohort_size)]))
      doses = array(c(doses, rep(administerdose[cohort_num],cohort_size)))
      # ------- adjustments for cohorts consideration ---------
      dat_now <- list(alpha_mean = -5, alpha_sd = 3, 
                      beta_mean = 3, beta_sd = 3,
                      gamma_mean = 0, gamma_sd = 2, 
                      zeta_mean = zeta_mean, zeta_sd = 2,
                      eta_mean = 0, eta_sd = 0.2, psi_mean = 0, psi_sd = 1,
                      real_doses = array(real_doses),
                      num_doses = num_doses, num_patients = num_patients,
                      tox = tox, eff = eff, eff_miss = eff_miss, doses = doses,
                      Lp=Lp, p1e=p1e, p2t=p2t, 
                      pe_ad=pe_ad, pt_ad=pt_ad, pitg=pitg, piel=piel)
      fit <- sampling(model_efftox_xc, data = dat_now, chains=chains)
      sum_list <- summary(fit)$summary
      for (num_d in 1:num_doses){
        exps_utility <- paste("utility","[",num_d,"]",sep="")
        exps_toxicity <- paste("prob_tox","[",num_d,"]",sep="")
        exps_efficacy <- paste("prob_eff","[",num_d,"]",sep="")
        
        efftox_ptox_mean[cohort_num,num_d] <- sum_list[exps_toxicity,'mean']
        efftox_ptox_se[cohort_num,num_d] <- sum_list[exps_toxicity,'se_mean']
        efftox_peff_mean[cohort_num,num_d] <- sum_list[exps_efficacy,'mean']
        efftox_peff_se[cohort_num,num_d] <- sum_list[exps_efficacy,'se_mean']
        efftox_utility_mean[cohort_num,num_d] <- sum_list[exps_utility,'mean']
        efftox_utility_se[cohort_num,num_d] <- sum_list[exps_utility,'se_mean']
      }
      # ----- dose admissibility -------
      prob_tox_samp <- (rstan::extract(fit, 'prob_tox')[[1]])
      prob_acc_tox <- colMeans(prob_tox_samp < dat_now$pitg)
      prob_eff_samp <- (rstan::extract(fit, 'prob_eff')[[1]])
      prob_acc_eff <- colMeans(prob_eff_samp > dat_now$piel)
      acceptable <- (prob_acc_eff > dat_now$pe_ad) & 
        (prob_acc_tox > dat_now$pt_ad) 
      if(sum(acceptable)){
        acceptable_doses <- (1:num_doses)[acceptable]
        administerdose_current <- acceptable_doses[which.max(efftox_utility_mean[cohort_num,acceptable_doses])]
        if(administerdose_current > administerdose[cohort_num]+1){
          administerdose[cohort_num+1] <- administerdose[cohort_num]+ 1
        }else{
          administerdose[cohort_num+1] <- administerdose_current
        }
        recommenddose[cohort_num] <- administerdose[cohort_num+1]
      }else{
        if((sum(prob_acc_eff > dat_now$pe_ad)==0)){
          stop_efficacy[cohort_num:length(cohorts)] <- 1 
          break
        }else{
          if((sum(prob_acc_tox > dat_now$pt_ad)==0)){
            stop_toxicity[cohort_num:length(cohorts)] <- 1 
            break
          }else{
            stop_other[cohort_num:length(cohorts)] <- 1
            break
          }
        }
      }
    }
    # ------------------------------- output format --------------------------------
    eff_outcome <- replace(eff_outcome, eff_outcome==2, '-')
    length(recommenddose) <- length(cohorts)
    length(doses) <- cohort_size*length(cohorts)
    # ------------------------------- output --------------------------------
    iter_out <- list('recommended_dose' = recommenddose,
                     'doses_given' = doses,
                     'toxicities' = as.vector(t(tox_outcome)),
                     'efficacies' = eff_outcome,
                     'efftox_ptox_mean' = efftox_ptox_mean,
                     'efftox_ptox_se' = efftox_ptox_se,
                     'efftox_peff_mean' = efftox_peff_mean,
                     'efftox_peff_se' = efftox_peff_se,
                     'efftox_utility_mean' = efftox_utility_mean,
                     'efftox_utility_se' = efftox_utility_se,
                     'stop_for_efficacy' = stop_efficacy,
                     'stop_for_toxicity' = stop_toxicity,
                     'stop_for_other' = stop_other)
    
    return(iter_out)
  }
  return(transpose(out))
}