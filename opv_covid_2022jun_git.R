######## setup ####
rm(list=ls())

setwd("~/Documents/04 OPV/SIR")

library(foreign)
library(data.table)
library(ggplot2)
library(deSolve)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(ggrepel)
options(scipen=999) 

######## model parameters ####
init_inf = 0.00175
seroprev = 0.007 - init_inf # assume May 11 (when the seroprev started), 0.7% May-June 
N = 1
pop_size <-  1384660352*0.68 

# target: seroprevalence 20-23%. overall prevalence 21.5 (27 cases to each confirmed case) 
# https://www.thehindu.com/news/national/coronavirus-one-in-five-indians-have-been-exposed-to-coronavirus-icmr-survey-finds/article33751083.ece 
# to early January. by end of Feb prob another few%. assume 25% 

time_step = 365 # 

R0 = 1.15
day_exposed <- 3
day_asymp <- 3
day_symp <- 4 
day_hosp <- 4 
day_icu <- 10 

perc_exp_asymp <- 1

perc_asymp_symp <- 0.074
perc_asymp_recover <- 1 - perc_asymp_symp

perc_symp_hosp <- 0.075 
perc_symp_dead <- 0.008
perc_symp_recover <- 1 - perc_symp_hosp - perc_symp_dead
perc_symp_recover

perc_hosp_icu <- 0.1 
perc_hosp_dead <- 0.07
perc_hosp_recover <- 1 - perc_hosp_icu - perc_hosp_dead
perc_hosp_recover

perc_icu_recover <- 0.85
perc_icu_dead <- 1 - perc_icu_recover

asymp = perc_exp_asymp/day_exposed # assume 100% becomes asymptomatic
symp = perc_asymp_symp/day_asymp
hosp = perc_symp_hosp/day_symp
icu = perc_hosp_icu/day_hosp

g1 = perc_asymp_recover/day_asymp
g2 = perc_symp_recover/day_symp
g3 = perc_hosp_recover/day_hosp
g4 = perc_icu_recover/day_icu

m_symp = perc_symp_dead/day_symp
m_hosp = perc_hosp_dead/day_hosp #
m_icu = perc_icu_dead/day_icu

b12_ratio = 0.75
b2 = R0 * (symp + g1) /(b12_ratio + (symp / (hosp + g2 + m_symp))) 
b1 = b2 * b12_ratio
b2
b1 

######## vaccine parameters #### 
vx_cov <- 0.3

# covid vaccine 
vx_covid_v1_day <- 28 
vx_covid_v2_day <- 0
vx_covid_v1_eff_inf <- 0.74
vx_covid_v1_eff_sev <- 0.95

vx_covid_v2_eff <- 0

vx_covid_delay <- 250

# opv 
vx_opv_delay_start <- 100
vx_opv_delay_rate <- 1 
# vx_opv_v2_day <- 1
vx_opv_v1_eff <- 0.2 # 0.47
# vx_opv_v2_eff <- 0 

######## Cost parameters #### 

## COVID
vx_covid_cost = 10  # 10 per dose https://www.reuters.com/article/uk-health-coronavirus-usa-johnsonandjohn-idUKKCN2511V7?edition-redirect=uk
vx_covid_cost_l = vx_covid_cost * (1-0.25)
vx_covid_cost_u = vx_covid_cost * (1+0.25)

vx_covid_cost_delivery = 1.07 + 0.41 * 258.811 / 251.107 # 1.49 
vx_covid_cost_delivery_l = 1.07 * (1-0.25 ) + 0.10 * 258.811 / 251.107 # 0.90
# or vx_covid_cost_delivery * (1-0.5)
vx_covid_cost_delivery_u = 1.07 * (1+0.25)+ 1.14 * 258.811 / 251.107  # 2.51 
# or vx_covid_cost_delivery * (1+0.5)


dose_covid <- 1
# 1.07 per dose supplied (Gavi/COVID) in India (table 12), 10% wastage, 2020 values 
# not including purchase 
# vaccine J&J 734 rs per dose = $10 USD 
# https://www.indiatoday.in/coronavirus-outbreak/vaccine-updates/story/here-is-what-one-dose-of-various-covid-19-vaccines-will-cost-in-india-1758527-2021-01-13

vx_covid_wastage <- 0.1 # WHO AMC 2021 wastage assumption 10% 
vx_covid_wastage_l <- vx_covid_wastage * (1-0.5) 
vx_covid_wastage_u <- vx_covid_wastage * (1+0.5) 

## OPV
vx_opv_cost = 0.15 # * 	258.811 / 255.657 # 2019 to 2020 USD, 0.1518 
vx_opv_cost_l = vx_opv_cost * (1-0.25) # 0.12
vx_opv_cost_u = vx_opv_cost * (1+0.25) # 0.19
vx_opv_cost_delivery = 0.96 # (0.56 + 0.77 + 1.56) / 3
  # 0.56 # 0.49 in 2011  258.811 /224.939 Kar 2014 India 
vx_opv_cost_delivery_l = 0.56 # vx_opv_cost_delivery  * (1-0.5)
vx_opv_cost_delivery_u = 1.56 # vx_opv_cost_delivery  * (1+0.5)
dose_opv <- 1

vx_opv_wastage <- 0.1 # WHO AMC 2021 wastage assumption 10% 
vx_opv_wastage_l <- vx_opv_wastage * (1-0.5) 
vx_opv_wastage_u <- vx_opv_wastage * (1+0.5)


######## BCA parameters #### 
gni <- 6920 # 2019 PPP 
gni_us <- 57900 # 66080
vsl_floor <- gni * 20 # 44600

vsl_1 <- 9400000 * (gni/gni_us) ^1.5 
vsl_2 <- gni * 100 
vsl_3 <- gni* 160 

vsl_floor > vsl_1 # if TRUE use vsl_floor instead
# vsl_1 <- vsl_floor
# 

######## SEIR model #### 

SEIR <- function(pars){
  
  # print(pars)
  times <- seq(from = 0, to = time_step, by = 1)    
  
  yinit <- c( S = N * (1 - seroprev - init_inf - vx_coverage), 
              S_vx_temp = N * vx_coverage,
              S_vb = 0, 
              S_v1 = 0, 
              S_v2 = 0, 
              E =  N*init_inf,
              E_vb = 0,
              E_v1 = 0, 
              E_v2 = 0,
              I_asymp = 0,
              I_asymp_vb = 0,
              I_asymp_v1 = 0, 
              I_asymp_v2 = 0,
              I_symp = 0, 
              I_symp_vb = 0,
              I_symp_v1 = 0, 
              I_symp_v2 = 0,
              I_hosp = 0, 
              I_hosp_vb = 0,
              I_hosp_v1 = 0, 
              I_hosp_v2 = 0,
              I_icu = 0, 
              I_icu_vb = 0,
              I_icu_v1 = 0, 
              I_icu_v2 = 0,
              R = seroprev,
              D = 0, 
              
              I_ever_inf = 0,
              I_ever_symp = 0, 
              I_ever_hosp = 0, 
              I_ever_icu = 0,
              I_ever_inf_novx = 0,
              I_ever_inf_vxtemp = 0, 
              I_ever_inf_vb = 0,
              I_ever_inf_v1 = 0,
              I_ever_inf_v2 = 0)
  
  SIR_model <- function(times, yinit, pars){
    
    with(as.list(c(yinit,pars)), {
      
      lambda <- beta1 * I_asymp / N + beta2 * I_symp / N +
        (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N +
        (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N +
        (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N
      
      # # vaccine impact only incidence but not transmission 
      # lambda <- beta1 * I_asymp / N + beta2 * I_symp / N + 
      #   beta1 * I_asymp_vb / N + beta2 * I_symp_vb / N + 
      #   beta1 * I_asymp_v1 / N + beta2 * I_symp_v1 / N + 
      #    beta1 * I_asymp_v2 / N + beta2 * I_symp_v2 / N  
      
      dS  <- - lambda * S 
      dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
      dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
      dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
      dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
      
      dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
      dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
      dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
      dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
      
      dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
      dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
      dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
      dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
      
      dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
      dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
      dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
      dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
      
      dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
      dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
      dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
      dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
      
      dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
      dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
      dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
      dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
      
      dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
        gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
        gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
        gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
      
      dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
        death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
        death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
      
      dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
      dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
      dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
      dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
      
      # where infection occurs 
      dI_ever_inf_novx <- lambda * S 
      dI_ever_inf_vxtemp <-  lambda * S_vx_temp
      dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
      dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
      dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
      
      return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                    dE, dE_vb, dE_v1, dE_v2, 
                    dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                    dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                    dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                    dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                    dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                    dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
    
  }
  
  results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
  results <- as.data.frame(results)
  
  return(results)
  
}

#### scenario 0: no vaccine ####
time_step <- 365
vx_coverage <- 0

s0_pars <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
             death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
             asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
             vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
             vx_susp_vb = 0, vx_susp_v1 = 0, vx_susp_v2 = 0,
             vx_delay_rate_vb = 0, vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, 
             vx_sev_vb = 0, vx_sev_v1 = 0, vx_sev_v2 = 0)

results_s0 <- SEIR(s0_pars)
results_s0 <- data.table(results_s0)
results_s0[time == time_step, I_ever_inf/ N]
results_s0[time == 242, I_ever_inf/ N] # 20-23 by Jan 8 
results_s0[time == 100, I_ever_inf/ N] # Aug 17 (98) -Sep 22 (134), should be about 7.1% 
results_s0[time == time_step, I_ever_symp/ N]/results_s0[time == time_step, I_ever_inf/ N]  # should be 7.4

results_s0[, I_daily := I_ever_inf - shift(I_ever_inf, type = "lag")]
results_s0[, I_symp_daily := I_ever_symp - shift(I_ever_symp, type = "lag")]
results_s0[, D_daily := D - shift(D, type = "lag")]

  s0_inf_max <- results_s0[, max(I_daily, na.rm=TRUE)]
  s0_inf_max * pop_size # 1,425,212 infections per day (1.43 million)
  s0_inf_max * 1000000 # 1513 
  
  s0_symp_max <- results_s0[, max(I_symp_daily, na.rm=TRUE)]
  s0_symp_max * pop_size # 105211 per day 
  s0_symp_max * 1000000 # 112
 
results_s0[D_daily == s0_inf_max, min(time)] #   should be around 130 days (May 11 - Sep 16)
results_s0[time == time_step, D/ N] / results_s0[time == time_step, I_ever_inf/ N] # IFR (among all infected) should be 0.09% (0.0009) but at the end of second round? 
results_s0[time == time_step, D/ N] / results_s0[time == time_step, I_ever_symp/ N] # CFR (among symptomatic):  1.34

results_s0[time == 130, D/ N] / results_s0[time == 130, I_ever_inf/ N] # IFR (among all infected) should be 0.09% (0.0009) but at the end of second round? 
results_s0[time == 130, D/ N] / results_s0[time == 130, I_ever_symp/ N] # CFR (among symptomatic):  1.34
results_s0[time == 242, D/ N] / results_s0[time == 242, I_ever_inf/ N] # IFR (among all infected) should be 0.09% (0.0009) but at the end of second round? 
results_s0[time == 242, D/ N] / results_s0[time == 242, I_ever_symp/ N] # CFR (among symptomatic):  1.34

results_s0[, I_per_inf := I_daily / (I_asymp + I_symp)]
results_s0[, S_total := S + S_vx_temp + S_vb + S_v1 + S_v2]
results_s0[, S_except_v := S + S_vx_temp]
results_s0[, E_total := E + E_vb + E_v1 + E_v2]
results_s0[, I_asymp_total := I_asymp + I_asymp_vb + I_asymp_v1 + I_asymp_v2]
results_s0[, I_symp_total := I_symp + I_symp_vb + I_symp_v1 + I_symp_v2]
results_s0[, I_hosp_total := I_hosp + I_hosp_vb + I_hosp_v1 + I_hosp_v2]
results_s0[, I_icu_total := I_icu + I_icu_vb + I_icu_v1 + I_icu_v2]

results_s0[, Re := R0 * S_except_v] # should be S_total or S + S_temp_vx? (exclude the effective ones?)
# results_s0[, S_total - S_except_vxeff]

#### scenario 1c: simultaneous administration #### 
# s1_delay <- c(25,50,100,200,300)
s1_delay <- c(5,10,seq(25,350,25))
vx_cov_v <- c(0.3,0.5)
opv_eff_v <- c(0.2,0.6)

data_s1c <- data.table()
for (v in vx_cov_v){
  for (t in s1_delay){
    print(t)
    
    vx_covid_delay_vb <- 1
    vx_covid_imple_delay <- t
    vx_covid_eff_delay <- vx_covid_v1_day
    
    vx_covid_inf_v1 <- vx_covid_v1_eff_inf
    vx_covid_susp_v1 <- vx_covid_v1_eff_inf
    vx_covid_sev_v1 <- vx_covid_v1_eff_sev 
    
    # step 1 
    time_step <- vx_covid_imple_delay
    init_inf = 0.00175
    seroprev = 0.007 - init_inf
    N = 1
    vx_coverage <- 0 # must be zero here 
    s1c_pars_step1 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                        death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                        asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                        vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
                        vx_susp_vb = 0, vx_susp_v1 = 0, vx_susp_v2 = 0,
                        vx_delay_rate_vb = 0, vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, 
                        vx_sev_vb = 0, vx_sev_v1 = 0, vx_sev_v2 = 0) # same as s0_pars
    
    results_s1c_step1 <- SEIR(s1c_pars_step1)
    results_s1c_step1 <- data.table(results_s1c_step1)
    # results_s1c_step1[time == vx_covid_imple_delay, .(I_ever_inf)] # 7.9
    ## check with results_s0
    # results_s0[time == vx_covid_imple_delay, .(I_ever_inf)] # 7.9
    
    # step 2
    time_step <- 365 - vx_covid_imple_delay
    vx_coverage <- v # V here 
    # vaccinate 30% of S, NOT people who have been infected 
    
    vx_covid_delay_v1 <- 1/vx_covid_eff_delay # 28 days 
    results_s1c_step1[time == vx_covid_imple_delay, ]
    
    SEIR_s1c_step2 <- function(pars){
      
      # print(pars)
      times <- seq(from = 0, to = time_step, by = 1)    
      
      yinit <- c( S = max(results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(S)] - vx_coverage,0), 
                  S_vx_temp = min(vx_coverage, results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(S)]), 
                  S_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(S_vb)], 
                  S_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(S_v1)], 
                  S_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(S_v2)], 
                  E =  results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(E)], 
                  E_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(E_vb)], 
                  E_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(E_v1)], 
                  E_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(E_v2)], 
                  I_asymp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_asymp)], 
                  I_asymp_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_asymp_vb)], 
                  I_asymp_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_asymp_v1)],  
                  I_asymp_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_asymp_v2)], 
                  I_symp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_symp)], 
                  I_symp_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_symp_vb)], 
                  I_symp_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_symp_v1)],  
                  I_symp_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_symp_v2)], 
                  I_hosp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_hosp)], 
                  I_hosp_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_hosp_vb)], 
                  I_hosp_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_hosp_v1)],  
                  I_hosp_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_hosp_v2)], 
                  I_icu = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_icu)], 
                  I_icu_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_icu_vb)], 
                  I_icu_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_icu_v1)],  
                  I_icu_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_icu_v2)], 
                  R = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(R)], 
                  D = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(D)], 
                  
                  I_ever_inf = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf)], 
                  I_ever_symp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_symp)], 
                  I_ever_hosp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_hosp)], 
                  I_ever_icu = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_icu)], 
                  I_ever_inf_novx = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf_novx)], 
                  I_ever_inf_vxtemp = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf_vxtemp)], 
                  I_ever_inf_vb = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf_vb)], 
                  I_ever_inf_v1 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf_v1)], 
                  I_ever_inf_v2 = results_s1c_step1[time == nrow(results_s1c_step1)-1, as.numeric(I_ever_inf_v2)])
      
      SIR_model <- function(times, yinit, pars){
        
        with(as.list(c(yinit,pars)), {
          
          lambda <- beta1 * I_asymp / N + beta2 * I_symp / N +
            (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N +
            (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N +
            (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N
          
          # # vaccine impact only incidence but not transmission 
          # lambda <- beta1 * I_asymp / N + beta2 * I_symp / N + 
          #   beta1 * I_asymp_vb / N + beta2 * I_symp_vb / N + 
          #   beta1 * I_asymp_v1 / N + beta2 * I_symp_v1 / N + 
          #   beta1 * I_asymp_v2 / N + beta2 * I_symp_v2 / N  
          # 
          
          dS  <- - lambda * S 
          dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
          dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
          dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
          dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
          
          dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
          dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
          dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
          dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
          
          dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
          dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
          dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
          dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
          
          dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
          dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
          dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
          dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
          
          dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
          dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
          dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
          dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
          
          dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
          dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
          dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
          dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
          
          dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
            gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
            gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
            gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
          
          dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
            death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
            death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
          
          dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
          dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
          dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
          dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
          
          # where infection occurs 
          dI_ever_inf_novx <- lambda * S 
          dI_ever_inf_vxtemp <-  lambda * S_vx_temp
          dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
          dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
          dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
          
          return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                        dE, dE_vb, dE_v1, dE_v2, 
                        dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                        dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                        dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                        dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                        dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                        dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
        
      }
      
      results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
      results <- as.data.frame(results)
      
      return(results)
      
    }
    s1c_pars_step2 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                        death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                        asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                        vx_inf_vb = 0, vx_inf_v1 = vx_covid_inf_v1, vx_inf_v2 = 0,
                        vx_susp_vb = 0, vx_susp_v1 = vx_covid_susp_v1, vx_susp_v2 = 0,
                        vx_delay_rate_vb = 1, vx_delay_rate_v1 = vx_covid_delay_v1, vx_delay_rate_v2 = 0, 
                        vx_sev_vb = 0, vx_sev_v1 = vx_covid_sev_v1 , vx_sev_v2 = 0)
    
    results_s1c_step2 <- SEIR_s1c_step2(s1c_pars_step2)
    results_s1c_step2 <- data.table(results_s1c_step2)
    # results_s1c_step2[time == time_step, .(I_ever_inf)] 
    # results_s0[time == 365, .(I_ever_inf)]  - results_s1c_step2[time == time_step, .(I_ever_inf)] # 100 days averted 11.6% 
    
    # merge two results 
    results_s1c <- rbind(results_s1c_step1, results_s1c_step2[time != 0, ])
    results_s1c[, .(time)]
    results_s1c[, time := c(0: (nrow(results_s1c) - 1))]
    
    time_step <- 365
    # results_s1c[time == time_step, .(I_ever_inf)]
    # results_s0[time == time_step, .(I_ever_inf)]  - results_s1c[time == time_step, .(I_ever_inf)]
    
    s1c_inf <- results_s1c[time == time_step, round(I_ever_inf/ N,2)] 
    # results_s1c[round(I_ever_inf,2) == s1c_inf, min(time)] # 155
    # results_s1c[time == time_step, I_ever_inf/ N]
    # results_s1c[time == time_step, I_ever_hosp/ N] 
    # results_s1c[time == time_step, D/ N] 
    # 
    # results_s1c[time == time_step, I_ever_symp/ N] / results_s1c[time == time_step, I_ever_inf/ N] 
    # results_s1c[time == time_step, D/ N] / results_s1c[time == time_step, I_ever_inf/ N] 
    # results_s1c[time == time_step, D/ N] / results_s1c[time == time_step, I_ever_symp/ N] 
    
    results_s1c[, D_daily := D - shift(D, type = "lag")]
    s1c_inf_max <- results_s1c[, max(D_daily, na.rm=TRUE)]
    # results_s1c[D_daily == s1c_inf_max, min(time)]  # peak at 112 days 
    
    results_s1c[, I_daily := I_ever_inf - shift(I_ever_inf, type = "lag")]
    results_s1c[, I_per_inf := I_daily / (I_asymp + I_symp)]
    
    results_s1c[, S_total := S + S_vx_temp + S_vb + S_v1 + S_v2]
    results_s1c[, S_except_v := S + S_vx_temp]
    results_s1c[, E_total := E + E_vb + E_v1 + E_v2]
    results_s1c[, I_asymp_total := I_asymp + I_asymp_vb + I_asymp_v1 + I_asymp_v2]
    results_s1c[, I_symp_total := I_symp + I_symp_vb + I_symp_v1 + I_symp_v2]
    results_s1c[, I_hosp_total := I_hosp + I_hosp_vb + I_hosp_v1 + I_hosp_v2]
    results_s1c[, I_icu_total := I_icu + I_icu_vb + I_icu_v1 + I_icu_v2]
    # results_s1c[, Re := R0 * S_except_v] # should be S_total or S + S_temp_vx? (exclude the effective ones?)
    # results_s1c[, S_total - S_except_vxeff]
    
    results_s1c[,delay := t]
    results_s1c[,vx_cov := v]
    
    # save outputs 
    data_s1c <- rbind(data_s1c, results_s1c)
  }
}

#### scenario 1: simultaneous administration #### 
data_s1 <- data.table()
for (v in vx_cov_v){
  for (t in s1_delay){
    for (o_eff in opv_eff_v){
      print(c(t,o_eff,v))
      
      vx_covid_delay_vb <- 1
      vx_covid_imple_delay <-t
      vx_covid_eff_delay <- vx_covid_v1_day
      
      vx_covid_inf_v1 <- vx_covid_v1_eff_inf
      vx_covid_susp_v1 <- vx_covid_v1_eff_inf
      vx_covid_sev_v1 <- vx_covid_v1_eff_sev 
      
      vx_opv_delay_vb <- 1 # immediate effect
      vx_opv_v1_eff <- o_eff
      
      # step 1 
      time_step <- vx_covid_imple_delay
      init_inf = 0.00175
      seroprev = 0.007 - init_inf
      N = 1
      vx_coverage <- 0 # must be zero here 
      s1_pars_step1 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                         death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                         asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                         vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
                         vx_susp_vb = 0, vx_susp_v1 = 0, vx_susp_v2 = 0,
                         vx_delay_rate_vb = 0, vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, 
                         vx_sev_vb = 0, vx_sev_v1 = 0, vx_sev_v2 = 0) # same as s0_pars
      
      results_s1_step1 <- SEIR(s1_pars_step1)
      results_s1_step1 <- data.table(results_s1_step1)

      # step 2
      time_step <- 365 - vx_covid_imple_delay
      vx_coverage <- v

      vx_covid_delay_v1 <- 1/vx_covid_eff_delay # 28 days 

      SEIR_s1_step2 <- function(pars){
        
        # print(pars)
        times <- seq(from = 0, to = time_step, by = 1)    
        
        yinit <- c( S = max(results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(S)] - vx_coverage,0), 
                    S_vx_temp = min(vx_coverage, results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(S)]), 
                    S_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(S_vb)], 
                    S_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(S_v1)], 
                    S_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(S_v2)], 
                    E =  results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(E)], 
                    E_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(E_vb)], 
                    E_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(E_v1)], 
                    E_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(E_v2)], 
                    I_asymp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_asymp)], 
                    I_asymp_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_asymp_vb)], 
                    I_asymp_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_asymp_v1)],  
                    I_asymp_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_asymp_v2)], 
                    I_symp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_symp)], 
                    I_symp_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_symp_vb)], 
                    I_symp_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_symp_v1)],  
                    I_symp_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_symp_v2)], 
                    I_hosp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_hosp)], 
                    I_hosp_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_hosp_vb)], 
                    I_hosp_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_hosp_v1)],  
                    I_hosp_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_hosp_v2)], 
                    I_icu = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_icu)], 
                    I_icu_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_icu_vb)], 
                    I_icu_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_icu_v1)],  
                    I_icu_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_icu_v2)], 
                    R = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(R)], 
                    D = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(D)], 
                    
                    I_ever_inf = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf)], 
                    I_ever_symp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_symp)], 
                    I_ever_hosp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_hosp)], 
                    I_ever_icu = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_icu)], 
                    I_ever_inf_novx = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf_novx)], 
                    I_ever_inf_vxtemp = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf_vxtemp)], 
                    I_ever_inf_vb = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf_vb)], 
                    I_ever_inf_v1 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf_v1)], 
                    I_ever_inf_v2 = results_s1_step1[time == nrow(results_s1_step1)-1, as.numeric(I_ever_inf_v2)])
        
        SIR_model <- function(times, yinit, pars){
          
          with(as.list(c(yinit,pars)), {
            
            lambda <- beta1 * I_asymp / N + beta2 * I_symp / N +
              (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N +
              (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N +
              (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N
            
            dS  <- - lambda * S 
            dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
            dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
            dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
            dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
            
            dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
            dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
            dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
            dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
            
            dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
            dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
            dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
            dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
            
            dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
            dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
            dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
            dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
            
            dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
            dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
            dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
            dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
            
            dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
            dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
            dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
            dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
            
            dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
              gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
              gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
              gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
            
            dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
              death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
              death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
            
            dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
            dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
            dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
            dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
            
            # where infection occurs 
            dI_ever_inf_novx <- lambda * S 
            dI_ever_inf_vxtemp <-  lambda * S_vx_temp
            dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
            dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
            dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
            
            return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                          dE, dE_vb, dE_v1, dE_v2, 
                          dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                          dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                          dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                          dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                          dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                          dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
          
        }
        
        results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
        results <- as.data.frame(results)
        
        return(results)
        
      }
      s1_pars_step2 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                         death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                         asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                         vx_delay_rate_vb = vx_opv_delay_vb, vx_delay_rate_v1 = vx_covid_delay_v1, vx_delay_rate_v2 = 0, 
                         vx_inf_vb = 0, vx_inf_v1 = vx_covid_inf_v1, vx_inf_v2 = 0,
                         vx_susp_vb = vx_opv_v1_eff, vx_susp_v1 = vx_covid_susp_v1, vx_susp_v2 = 0,
                         vx_sev_vb = vx_opv_v1_eff, vx_sev_v1 = vx_covid_sev_v1 , vx_sev_v2 = 0)
      results_s1_step2 <- SEIR_s1_step2(s1_pars_step2)
      results_s1_step2 <- data.table(results_s1_step2)
      # results_s1_step2[time == time_step, .(I_ever_inf)] # 11.8
      # results_s0[time == 365, .(I_ever_inf)]  - results_s1_step2[time == time_step, .(I_ever_inf)]
      
      # merge two results 
      results_s1 <- rbind(results_s1_step1, results_s1_step2[time != 0, ])
      results_s1[, .(time)]
      results_s1[, time := c(0: (nrow(results_s1) - 1))]
      
      time_step <- 365

      s1_inf <- results_s1[time == time_step, round(I_ever_inf/ N,2)] 

      results_s1[, D_daily := D - shift(D, type = "lag")]
      s1_inf_max <- results_s1[, max(D_daily, na.rm=TRUE)]
      # results_s1[D_daily == s1_inf_max, min(time)]  # 
      # 
      results_s1[, I_daily := I_ever_inf - shift(I_ever_inf, type = "lag")]
      results_s1[, I_per_inf := I_daily / (I_asymp + I_symp)]
      
      results_s1[, S_total := S + S_vx_temp + S_vb + S_v1 + S_v2]
      results_s1[, S_except_v := S + S_vx_temp]
      results_s1[, E_total := E + E_vb + E_v1 + E_v2]
      results_s1[, I_asymp_total := I_asymp + I_asymp_vb + I_asymp_v1 + I_asymp_v2]
      results_s1[, I_symp_total := I_symp + I_symp_vb + I_symp_v1 + I_symp_v2]
      results_s1[, I_hosp_total := I_hosp + I_hosp_vb + I_hosp_v1 + I_hosp_v2]
      results_s1[, I_icu_total := I_icu + I_icu_vb + I_icu_v1 + I_icu_v2]

      results_s1[, delay := t]
      results_s1[, vx_cov := v]
      results_s1[, opv_eff := o_eff]
      
      # save outputs 
      data_s1 <- rbind(data_s1, results_s1)
      
    }
  }
}


#### scenario 2c: different day administration #### 
s2_opv_delay <- c(25,50,100)
s2_covid_delay <- seq(25,300,25) # c(50,100,150)
vx_cov_v <- c(0.3,0.5)
opv_eff_v <- c(0.2,0.6)

data_s2c <- data.table()

for (v in vx_cov_v){
  for (t in s2_opv_delay){
    for (d in s2_covid_delay){
      if(t + d <= 365) {
        print(c(d, t, v))
        
        vx_covid_delay_vb <- 1
        vx_covid_imple_delay <- (d + t)
        vx_covid_eff_delay <- vx_covid_v1_day
        
        vx_covid_inf_v1 <- vx_covid_v1_eff_inf
        vx_covid_susp_v1 <- vx_covid_v1_eff_inf
        vx_covid_sev_v1 <- vx_covid_v1_eff_sev 
        
        # step 1 
        time_step <- vx_covid_imple_delay 
        init_inf = 0.00175
        seroprev = 0.007 - init_inf
        N = 1
        vx_coverage <- 0 # must be zero here 
        s2c_pars_step1 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                            death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                            asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                            vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
                            vx_susp_vb = 0, vx_susp_v1 = 0, vx_susp_v2 = 0,
                            vx_delay_rate_vb = 0, vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, 
                            vx_sev_vb = 0, vx_sev_v1 = 0, vx_sev_v2 = 0) # same as s0_pars
        results_s2c_step1 <- SEIR(s2c_pars_step1)
        results_s2c_step1 <- data.table(results_s2c_step1)

        # step 2
        time_step <- 365 - vx_covid_imple_delay
        vx_coverage <- v
        vx_covid_delay_v1 <- 1/vx_covid_eff_delay # 28 days 

        SEIR_s2c_step2 <- function(pars){
          
          # print(pars)
          times <- seq(from = 0, to = time_step, by = 1)    
          
          yinit <- c( S = max(results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(S)] - vx_coverage,0), 
                      S_vx_temp = min(vx_coverage, results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(S)]), 
                      S_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(S_vb)], 
                      S_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(S_v1)], 
                      S_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(S_v2)], 
                      E =  results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(E)], 
                      E_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(E_vb)], 
                      E_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(E_v1)], 
                      E_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(E_v2)], 
                      I_asymp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_asymp)], 
                      I_asymp_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_asymp_vb)], 
                      I_asymp_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_asymp_v1)],  
                      I_asymp_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_asymp_v2)], 
                      I_symp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_symp)], 
                      I_symp_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_symp_vb)], 
                      I_symp_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_symp_v1)],  
                      I_symp_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_symp_v2)], 
                      I_hosp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_hosp)], 
                      I_hosp_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_hosp_vb)], 
                      I_hosp_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_hosp_v1)],  
                      I_hosp_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_hosp_v2)], 
                      I_icu = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_icu)], 
                      I_icu_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_icu_vb)], 
                      I_icu_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_icu_v1)],  
                      I_icu_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_icu_v2)], 
                      R = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(R)], 
                      D = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(D)], 
                      
                      I_ever_inf = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf)], 
                      I_ever_symp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_symp)], 
                      I_ever_hosp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_hosp)], 
                      I_ever_icu = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_icu)], 
                      I_ever_inf_novx = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf_novx)], 
                      I_ever_inf_vxtemp = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf_vxtemp)], 
                      I_ever_inf_vb = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf_vb)], 
                      I_ever_inf_v1 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf_v1)], 
                      I_ever_inf_v2 = results_s2c_step1[time == nrow(results_s2c_step1)-1, as.numeric(I_ever_inf_v2)])
          
          SIR_model <- function(times, yinit, pars){
            
            with(as.list(c(yinit,pars)), {
              
              lambda <- beta1 * I_asymp / N + beta2 * I_symp / N +
                (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N +
                (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N +
                (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N
  
              dS  <- - lambda * S 
              dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
              dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
              dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
              dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
              
              dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
              dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
              dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
              dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
              
              dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
              dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
              dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
              dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
              
              dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
              dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
              dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
              dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
              
              dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
              dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
              dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
              dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
              
              dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
              dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
              dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
              dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
              
              dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
                gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
                gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
                gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
              
              dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
                death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
                death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
              
              dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
              dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
              dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
              dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
              
              # where infection occurs 
              dI_ever_inf_novx <- lambda * S 
              dI_ever_inf_vxtemp <-  lambda * S_vx_temp
              dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
              dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
              dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
              
              return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                            dE, dE_vb, dE_v1, dE_v2, 
                            dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                            dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                            dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                            dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                            dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                            dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
            
          }
          
          results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
          results <- as.data.frame(results)
          
          return(results)
          
        }
        s2c_pars_step2 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                            death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                            asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                            vx_inf_vb = 0, vx_inf_v1 = vx_covid_inf_v1, vx_inf_v2 = 0,
                            vx_susp_vb = 0, vx_susp_v1 = vx_covid_susp_v1, vx_susp_v2 = 0,
                            vx_delay_rate_vb = 1, vx_delay_rate_v1 = vx_covid_delay_v1, vx_delay_rate_v2 = 0, 
                            vx_sev_vb = 0, vx_sev_v1 = vx_covid_sev_v1 , vx_sev_v2 = 0)
        
        results_s2c_step2 <- SEIR_s2c_step2(s2c_pars_step2)
        results_s2c_step2 <- data.table(results_s2c_step2)
        # results_s2c_step2[time == time_step, .(I_ever_inf)] 
        # results_s0[time == 365, .(I_ever_inf)]  - results_s2c_step2[time == time_step, .(I_ever_inf)]
        
        # merge two results 
        results_s2c <- rbind(results_s2c_step1, results_s2c_step2[time != 0, ])
        results_s2c[, .(time)]
        results_s2c[, time := c(0: (nrow(results_s2c) - 1))]
        
        time_step <- 365

        s2c_inf <- results_s2c[time == time_step, round(I_ever_inf/ N,2)] 

        results_s2c[, D_daily := D - shift(D, type = "lag")]
        s2c_inf_max <- results_s2c[, max(D_daily, na.rm=TRUE)]
        # results_s2c[D_daily == s2c_inf_max, min(time)]  # peak at 59 days 
        
        results_s2c[, I_daily := I_ever_inf - shift(I_ever_inf, type = "lag")]
        results_s2c[, I_per_inf := I_daily / (I_asymp + I_symp)]
        
        results_s2c[, S_total := S + S_vx_temp + S_vb + S_v1 + S_v2]
        results_s2c[, S_except_v := S + S_vx_temp]
        results_s2c[, E_total := E + E_vb + E_v1 + E_v2]
        results_s2c[, I_asymp_total := I_asymp + I_asymp_vb + I_asymp_v1 + I_asymp_v2]
        results_s2c[, I_symp_total := I_symp + I_symp_vb + I_symp_v1 + I_symp_v2]
        results_s2c[, I_hosp_total := I_hosp + I_hosp_vb + I_hosp_v1 + I_hosp_v2]
        results_s2c[, I_icu_total := I_icu + I_icu_vb + I_icu_v1 + I_icu_v2]

        results_s2c[, opv_delay := t]
        results_s2c[, covid_delay := d]
        results_s2c[, vx_cov := v]
        
        # save outputs 
        data_s2c <- rbind(data_s2c, results_s2c)
        
      }
      }
    }
}

#### scenario 2: different day administration #### 
data_s2 <- data.table()

for (v in vx_cov_v){
  for (o_eff in opv_eff_v){
    for (t in s2_opv_delay){
      for (d in s2_covid_delay){
        print(c(d, t, v, o_eff))
        
        if(t + d <= 365) {
          vx_covid_delay_vb <- 1
          vx_covid_imple_delay <- d
          vx_covid_eff_delay <- vx_covid_v1_day
          vx_covid_inf_v1 <- vx_covid_v1_eff_inf
          vx_covid_susp_v1 <- vx_covid_v1_eff_inf
          vx_covid_sev_v1 <- vx_covid_v1_eff_sev 
          
          vx_opv_day_start <- t
          vx_opv_delay_vb <- vx_opv_delay_rate # 1/35
          vx_opv_v1_eff <- o_eff
          
          # step 1: no vaccine
          time_step <- vx_opv_day_start
          init_inf = 0.00175
          seroprev = 0.007 - init_inf
          N = 1
          vx_coverage <- 0 # must be zero here 
          s2_pars_step1 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                             death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                             asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                             vx_delay_rate_vb = 0, vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, 
                             vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
                             vx_susp_vb = 0, vx_susp_v1 = 0, vx_susp_v2 = 0,
                             vx_sev_vb = 0, vx_sev_v1 = 0, vx_sev_v2 = 0) # same as s0_pars
          results_s2_step1 <- SEIR(s2_pars_step1)
          results_s2_step1 <- data.table(results_s2_step1)

          # step 2: OPV starting day 100
          time_step <- vx_covid_imple_delay
          vx_coverage <- v
          
          s2_pars_step2 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                             death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                             asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                             vx_delay_rate_vb = vx_opv_delay_vb , vx_delay_rate_v1 = 0, vx_delay_rate_v2 = 0, # 
                             vx_inf_vb = 0, vx_inf_v1 = 0, vx_inf_v2 = 0,
                             vx_susp_vb = vx_opv_v1_eff, vx_susp_v1 = 0, vx_susp_v2 = 0,
                             vx_sev_vb = vx_opv_v1_eff, vx_sev_v1 = 0, vx_sev_v2 = 0)
          SEIR_s2_step2 <- function(pars){
            
            # print(pars)
            times <- seq(from = 0, to = time_step, by = 1)    
            
            yinit <- c(S = max(results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(S)] - vx_coverage,0), 
                       S_vx_temp = min(vx_coverage, results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(S)]), 
                       S_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(S_vb)], 
                       S_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(S_v1)], 
                       S_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(S_v2)], 
                       E =  results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(E)], 
                       E_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(E_vb)], 
                       E_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(E_v1)], 
                       E_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(E_v2)], 
                       I_asymp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_asymp)], 
                       I_asymp_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_asymp_vb)], 
                       I_asymp_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_asymp_v1)],  
                       I_asymp_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_asymp_v2)], 
                       I_symp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_symp)], 
                       I_symp_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_symp_vb)], 
                       I_symp_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_symp_v1)],  
                       I_symp_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_symp_v2)], 
                       I_hosp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_hosp)], 
                       I_hosp_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_hosp_vb)], 
                       I_hosp_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_hosp_v1)],  
                       I_hosp_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_hosp_v2)], 
                       I_icu = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_icu)], 
                       I_icu_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_icu_vb)], 
                       I_icu_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_icu_v1)],  
                       I_icu_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_icu_v2)], 
                       R = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(R)], 
                       D = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(D)], 
                       
                       I_ever_inf = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf)], 
                       I_ever_symp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_symp)], 
                       I_ever_hosp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_hosp)], 
                       I_ever_icu = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_icu)], 
                       I_ever_inf_novx = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf_novx)], 
                       I_ever_inf_vxtemp = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf_vxtemp)], 
                       I_ever_inf_vb = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf_vb)], 
                       I_ever_inf_v1 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf_v1)], 
                       I_ever_inf_v2 = results_s2_step1[time == nrow(results_s2_step1)-1, as.numeric(I_ever_inf_v2)])
            
            SIR_model <- function(times, yinit, pars){
              
              with(as.list(c(yinit,pars)), {
                
                lambda <- beta1 * I_asymp / N + beta2 * I_symp / N + 
                  (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N + 
                  (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N + 
                  (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N  
                
                dS  <- - lambda * S 
                dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
                dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
                dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
                dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
                
                dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
                dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
                dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
                dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
                
                dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
                dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
                dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
                dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
                
                dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
                dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
                dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
                dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
                
                dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
                dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
                dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
                dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
                
                dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
                dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
                dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
                dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
                
                dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
                  gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
                  gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
                  gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
                
                dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
                  death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
                  death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
                
                dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
                dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
                dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
                dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
                
                # where infection occurs 
                dI_ever_inf_novx <- lambda * S 
                dI_ever_inf_vxtemp <-  lambda * S_vx_temp
                dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
                dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
                dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
                
                return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                              dE, dE_vb, dE_v1, dE_v2, 
                              dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                              dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                              dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                              dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                              dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                              dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
              
            }
            
            results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
            results <- as.data.frame(results)
            
            return(results)
            
          }
          results_s2_step2 <- SEIR_s2_step2(s2_pars_step2)
          results_s2_step2 <- data.table(results_s2_step2)

          # step 3: COVID starting later
          time_step <- 365 - vx_covid_imple_delay - vx_opv_day_start
          vx_covid_delay_v1 <- 1/vx_covid_eff_delay # 28 days 
          s2_pars_step3 <- c(beta1 = b1, beta2 = b2,  gamma1 = g1, gamma2 = g2, gamma3 = g3, gamma4 = g4, 
                             death_symp = m_symp, death_hosp = m_hosp, death_icu = m_icu, 
                             asymp_rate = asymp, symp_rate = symp, hosp_rate = hosp, icu_rate = icu,
                             vx_delay_rate_vb = vx_opv_delay_vb, vx_delay_rate_v1 = vx_covid_delay_v1, vx_delay_rate_v2 = 0, 
                             vx_inf_vb = 0, vx_inf_v1 = vx_covid_inf_v1, vx_inf_v2 = 0,
                             vx_susp_vb = vx_opv_v1_eff, vx_susp_v1 = vx_covid_susp_v1, vx_susp_v2 = 0,
                             vx_sev_vb = vx_opv_v1_eff, vx_sev_v1 = vx_covid_sev_v1, vx_sev_v2 = 0)
          SEIR_s2_step3 <- function(pars){
            
            # print(pars)
            times <- seq(from = 0, to = time_step, by = 1)    
            
            yinit <- c( S = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(S)], 
                        S_vx_temp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(S_vx_temp)], 
                        S_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(S_vb)], 
                        S_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(S_v1)], 
                        S_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(S_v2)], 
                        E =  results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(E)], 
                        E_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(E_vb)], 
                        E_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(E_v1)], 
                        E_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(E_v2)], 
                        I_asymp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_asymp)], 
                        I_asymp_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_asymp_vb)], 
                        I_asymp_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_asymp_v1)],  
                        I_asymp_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_asymp_v2)], 
                        I_symp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_symp)], 
                        I_symp_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_symp_vb)], 
                        I_symp_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_symp_v1)],  
                        I_symp_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_symp_v2)], 
                        I_hosp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_hosp)], 
                        I_hosp_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_hosp_vb)], 
                        I_hosp_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_hosp_v1)],  
                        I_hosp_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_hosp_v2)], 
                        I_icu = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_icu)], 
                        I_icu_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_icu_vb)], 
                        I_icu_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_icu_v1)],  
                        I_icu_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_icu_v2)], 
                        R = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(R)], 
                        D = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(D)], 
                        
                        I_ever_inf = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf)], 
                        I_ever_symp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_symp)], 
                        I_ever_hosp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_hosp)], 
                        I_ever_icu = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_icu)], 
                        I_ever_inf_novx = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf_novx)], 
                        I_ever_inf_vxtemp = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf_vxtemp)], 
                        I_ever_inf_vb = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf_vb)], 
                        I_ever_inf_v1 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf_v1)], 
                        I_ever_inf_v2 = results_s2_step2[time == nrow(results_s2_step2)-1, as.numeric(I_ever_inf_v2)])
            
            SIR_model <- function(times, yinit, pars){
              
              with(as.list(c(yinit,pars)), {
                
                lambda <- beta1 * I_asymp / N + beta2 * I_symp / N + 
                  (1-vx_inf_vb) * beta1 * I_asymp_vb / N + (1-vx_inf_vb) * beta2 * I_symp_vb / N + 
                  (1-vx_inf_v1) * beta1 * I_asymp_v1 / N + (1-vx_inf_v1) * beta2 * I_symp_v1 / N + 
                  (1-vx_inf_v2) * beta1 * I_asymp_v2 / N + (1-vx_inf_v2) * beta2 * I_symp_v2 / N  
                
                dS  <- - lambda * S 
                dS_vx_temp <- - lambda * S_vx_temp - vx_delay_rate_vb * S_vx_temp 
                dS_vb <- vx_delay_rate_vb * S_vx_temp - (1-vx_susp_vb) * lambda * S_vb - vx_delay_rate_v1 * S_vb
                dS_v1 <- vx_delay_rate_v1 * S_vb - (1-vx_susp_v1) * lambda * S_v1 - vx_delay_rate_v2 * S_v1
                dS_v2 <- vx_delay_rate_v2 * S_v1 - (1-vx_susp_v2) * lambda * S_v2    
                
                dE <- lambda * S + lambda * S_vx_temp - asymp_rate * E  
                dE_vb <- (1-vx_susp_vb) * lambda * S_vb - (1-vx_sev_vb) * asymp_rate * E_vb
                dE_v1 <- (1-vx_susp_v1) * lambda * S_v1 - (1-vx_sev_v1) * asymp_rate * E_v1
                dE_v2 <- (1-vx_susp_v2) * lambda * S_v2 - (1-vx_sev_v2) * asymp_rate * E_v2
                
                dI_asymp <- asymp_rate * E - gamma1 * I_asymp - symp_rate * I_asymp
                dI_asymp_vb <- (1-vx_sev_vb) * asymp_rate * E_vb - (1+vx_sev_vb) * gamma1 * I_asymp_vb - (1-vx_sev_vb) * symp_rate * I_asymp_vb
                dI_asymp_v1 <- (1-vx_sev_v1) * asymp_rate * E_v1 - (1+vx_sev_v1) * gamma1 * I_asymp_v1 - (1-vx_sev_v1) * symp_rate * I_asymp_v1
                dI_asymp_v2 <- (1-vx_sev_v2) * asymp_rate * E_v2 - (1+vx_sev_v2) * gamma1 * I_asymp_v2 - (1-vx_sev_v2) * symp_rate * I_asymp_v2
                
                dI_symp <- symp_rate * I_asymp - gamma2 * I_symp - hosp_rate * I_symp - death_symp * I_symp
                dI_symp_vb <- (1-vx_sev_vb) * symp_rate * I_asymp_vb - (1+vx_sev_vb) * gamma2 * I_symp_vb - (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1-vx_sev_vb) * death_symp * I_symp_vb
                dI_symp_v1 <- (1-vx_sev_v1) * symp_rate * I_asymp_v1 - (1+vx_sev_v1) * gamma2 * I_symp_v1 - (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1-vx_sev_v1) * death_symp * I_symp_v1
                dI_symp_v2 <- (1-vx_sev_v2) * symp_rate * I_asymp_v2 - (1+vx_sev_v2) * gamma2 * I_symp_v2 - (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1-vx_sev_v2) * death_symp * I_symp_v2
                
                dI_hosp <- hosp_rate * I_symp - gamma3 * I_hosp - icu_rate * I_hosp - death_hosp * I_hosp
                dI_hosp_vb <- (1-vx_sev_vb) * hosp_rate * I_symp_vb - (1+vx_sev_vb) * gamma3 * I_hosp_vb - (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1-vx_sev_vb) * death_hosp * I_hosp_vb
                dI_hosp_v1 <- (1-vx_sev_v1) * hosp_rate * I_symp_v1 - (1+vx_sev_v1) * gamma3 * I_hosp_v1 - (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1-vx_sev_v1) * death_hosp * I_hosp_v1
                dI_hosp_v2 <- (1-vx_sev_v2) * hosp_rate * I_symp_v2 - (1+vx_sev_v2) * gamma3 * I_hosp_v2 - (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1-vx_sev_v2) * death_hosp * I_hosp_v2
                
                dI_icu <- icu_rate * I_hosp - gamma4 * I_icu - death_icu * I_icu
                dI_icu_vb <- (1-vx_sev_vb) * icu_rate * I_hosp_vb - (1+vx_sev_vb) * gamma4 * I_icu_vb - (1-vx_sev_vb) * death_icu * I_icu_vb
                dI_icu_v1 <- (1-vx_sev_v1) * icu_rate * I_hosp_v1 - (1+vx_sev_v1) * gamma4 * I_icu_v1 - (1-vx_sev_v1) * death_icu * I_icu_v1
                dI_icu_v2 <- (1-vx_sev_v2) * icu_rate * I_hosp_v2 - (1+vx_sev_v2) * gamma4 * I_icu_v2 - (1-vx_sev_v2) * death_icu * I_icu_v2
                
                dR <- gamma1 * I_asymp + (1+vx_sev_vb) * gamma1 * I_asymp_vb + (1+vx_sev_v1) * gamma1 * I_asymp_v1 + (1+vx_sev_v2) * gamma1 * I_asymp_v2 + 
                  gamma2 * I_symp  + (1+vx_sev_vb) * gamma2 * I_symp_vb  + (1+vx_sev_v1) * gamma2 * I_symp_v1  + (1+vx_sev_v2) * gamma2 * I_symp_v2  +
                  gamma3 * I_hosp  + (1+vx_sev_vb) * gamma3 * I_hosp_vb  + (1+vx_sev_v1) * gamma3 * I_hosp_v1  + (1+vx_sev_v2) * gamma3 * I_hosp_v2  +
                  gamma4 * I_icu   + (1+vx_sev_vb) * gamma4 * I_icu_vb   + (1+vx_sev_v1) * gamma4 * I_icu_v1   + (1+vx_sev_v2) * gamma4 * I_icu_v2 
                
                dD <- death_symp * I_symp + (1-vx_sev_vb) * death_symp * I_symp_vb + (1-vx_sev_v1) * death_symp * I_symp_v1 + (1-vx_sev_v2) * death_symp * I_symp_v2 +
                  death_hosp * I_hosp + (1-vx_sev_vb) * death_hosp * I_hosp_vb + (1-vx_sev_v1) * death_hosp * I_hosp_v1 + (1-vx_sev_v2) * death_hosp * I_hosp_v2 +
                  death_icu * I_icu   + (1-vx_sev_vb) * death_icu * I_icu_vb   + (1-vx_sev_v1) * death_icu * I_icu_v1   + (1-vx_sev_v2) * death_icu * I_icu_v2
                
                dI_ever_inf <- lambda * S +  lambda * S_vx_temp + (1-vx_susp_vb) * lambda * S_vb + (1-vx_susp_v1) * lambda * S_v1 + (1-vx_susp_v2) * lambda * S_v2
                dI_ever_symp <- symp_rate * I_asymp + (1-vx_sev_vb) * symp_rate * I_asymp_vb + (1-vx_sev_v1) * symp_rate * I_asymp_v1 + (1-vx_sev_v2) * symp_rate * I_asymp_v2
                dI_ever_hosp <- hosp_rate * I_symp  + (1-vx_sev_vb) * hosp_rate * I_symp_vb  + (1-vx_sev_v1) * hosp_rate * I_symp_v1  + (1-vx_sev_v2) * hosp_rate * I_symp_v2
                dI_ever_icu <-  icu_rate * I_hosp   + (1-vx_sev_vb) * icu_rate * I_hosp_vb   + (1-vx_sev_v1) * icu_rate * I_hosp_v1   + (1-vx_sev_v2) * icu_rate * I_hosp_v2
                
                # where infection occurs 
                dI_ever_inf_novx <- lambda * S 
                dI_ever_inf_vxtemp <-  lambda * S_vx_temp
                dI_ever_inf_vb <- (1-vx_susp_vb) * lambda * S_vb
                dI_ever_inf_v1 <- (1-vx_susp_v1) * lambda * S_v1
                dI_ever_inf_v2 <- (1-vx_susp_v2) * lambda * S_v2
                
                return(list(c(dS, dS_vx_temp, dS_vb, dS_v1, dS_v2, 
                              dE, dE_vb, dE_v1, dE_v2, 
                              dI_asymp, dI_asymp_vb, dI_asymp_v1, dI_asymp_v2, 
                              dI_symp, dI_symp_vb, dI_symp_v1, dI_symp_v2, 
                              dI_hosp, dI_hosp_vb, dI_hosp_v1, dI_hosp_v2, 
                              dI_icu, dI_icu_vb, dI_icu_v1, dI_icu_v2, 
                              dR, dD, dI_ever_inf, dI_ever_symp, dI_ever_hosp, dI_ever_icu,
                              dI_ever_inf_novx, dI_ever_inf_vxtemp, dI_ever_inf_vb, dI_ever_inf_v1, dI_ever_inf_v2)))})
              
            }
            
            results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
            results <- as.data.frame(results)
            
            return(results)
            
          }
          results_s2_step3 <- SEIR_s2_step3(s2_pars_step3)
          results_s2_step3 <- data.table(results_s2_step3)

          # merge all results 
          results_s2 <- rbind(results_s2_step1, results_s2_step2[time!= 0, ], results_s2_step3[time!= 0, ])
          results_s2[, .(time)]
          results_s2[, time := c(0: (nrow(results_s2) - 1))]
          
          time_step <- 365
          s2_inf <- results_s2[time == time_step, round(I_ever_inf/ N,2)] 

          results_s2[, D_daily := D - shift(D, type = "lag")]
          s2_inf_max <- results_s2[, max(D_daily, na.rm=TRUE)]
          # results_s2[D_daily == s2_inf_max, min(time)]  # peak at 91 days 
          
          results_s2[, I_daily := I_ever_inf - shift(I_ever_inf, type = "lag")]
          results_s2[, I_per_inf := I_daily / (I_asymp + I_symp)]
          
          results_s2[, S_total := S + S_vx_temp + S_vb + S_v1 + S_v2]
          results_s2[, S_except_v := S + S_vx_temp]
          results_s2[, E_total := E + E_vb + E_v1 + E_v2]
          results_s2[, I_asymp_total := I_asymp + I_asymp_vb + I_asymp_v1 + I_asymp_v2]
          results_s2[, I_symp_total := I_symp + I_symp_vb + I_symp_v1 + I_symp_v2]
          results_s2[, I_hosp_total := I_hosp + I_hosp_vb + I_hosp_v1 + I_hosp_v2]
          results_s2[, I_icu_total := I_icu + I_icu_vb + I_icu_v1 + I_icu_v2]

          results_s2[, opv_delay := t]
          results_s2[, covid_delay := d]
          results_s2[, vx_cov := v]
          results_s2[, opv_eff := o_eff]
          
          # save outputs 
          data_s2 <- rbind(data_s2, results_s2)
        }
  
      }
    }
  }
}

#### results ####
# S1 
summary_s1 <- merge(data_s1[time == 365,.(delay, vx_cov, opv_eff, I_ever_inf, I_ever_symp,I_ever_hosp,I_ever_icu, D)],
                    data_s1c[time == 365,.(delay,vx_cov, I_ever_inf, I_ever_symp,I_ever_hosp,I_ever_icu, D)], by=c("delay","vx_cov"), all.x = TRUE)
names(summary_s1) <- c("delay","vx_cov","opv_eff",
                       "s1_ever_inf","s1_ever_symp","s1_ever_hosp","s1_ever_icu","s1_D","s1c_ever_inf","s1c_ever_symp","s1c_ever_hosp","s1c_ever_icu","s1c_D")
summary_s1[vx_cov == 0.3, s1c_cost := pop_size * 0.3 * dose_covid * (vx_covid_cost + vx_covid_cost_delivery) * (1 + vx_covid_wastage)]
summary_s1[vx_cov == 0.3, s1_cost := pop_size * 0.3 * (dose_covid * (vx_covid_cost + vx_covid_cost_delivery) + dose_opv * (vx_opv_cost + vx_opv_cost_delivery)) * (1 + vx_covid_wastage)]
summary_s1[vx_cov == 0.5, s1c_cost := pop_size * 0.5 * dose_covid * (vx_covid_cost + vx_covid_cost_delivery) * (1 + vx_covid_wastage)]
summary_s1[vx_cov == 0.5, s1_cost := pop_size * 0.5 * (dose_covid * (vx_covid_cost + vx_covid_cost_delivery) + dose_opv * (vx_opv_cost + vx_opv_cost_delivery)) * (1 + vx_covid_wastage)]

summary_s1[,.(s1c_D - s1_D)]
summary_s1[, dati := (s1c_D - s1_D) * 1000 / vx_cov]
summary_s1[, ce := (s1_cost - s1c_cost)/ (pop_size * (s1c_D - s1_D))]
summary_s1[, bc := pop_size * (s1c_D - s1_D) * vsl_1 / (s1_cost - s1c_cost)]
summary_s1
write.csv(summary_s1,"summary_s1.csv")

# S2
summary_s2 <- merge(data_s2[time == 365,.(opv_delay, covid_delay, vx_cov, opv_eff, I_ever_inf, I_ever_symp,I_ever_hosp,I_ever_icu, D)],
                    data_s2c[time == 365,.(opv_delay, covid_delay, vx_cov,  I_ever_inf, I_ever_symp,I_ever_hosp,I_ever_icu, D)], by=c("opv_delay","covid_delay","vx_cov"), all.x = TRUE)
names(summary_s2) <- c("opv_delay","covid_delay","vx_cov","opv_eff",
                       "s2_ever_inf","s2_ever_symp","s2_ever_hosp","s2_ever_icu","s2_D","s2c_ever_inf","s2c_ever_symp","s2c_ever_hosp","s2c_ever_icu","s2c_D")
summary_s2[vx_cov == 0.3, s2c_cost := pop_size * 0.3 * dose_covid * (vx_covid_cost + vx_covid_cost_delivery) * (1 + vx_covid_wastage)]
summary_s2[vx_cov == 0.3, s2_cost := pop_size * 0.3 * (dose_covid * (vx_covid_cost + vx_covid_cost_delivery) + dose_opv * (vx_opv_cost + vx_opv_cost_delivery)) * (1 + vx_covid_wastage)]
summary_s2[vx_cov == 0.5, s2c_cost := pop_size * 0.5 * dose_covid * (vx_covid_cost + vx_covid_cost_delivery) * (1 + vx_covid_wastage)]
summary_s2[vx_cov == 0.5, s2_cost := pop_size * 0.5 * (dose_covid * (vx_covid_cost + vx_covid_cost_delivery) + dose_opv * (vx_opv_cost + vx_opv_cost_delivery)) * (1 + vx_covid_wastage)]

summary_s2[,.(s2c_D - s2_D)]
summary_s2[, dati := (s2c_D - s2_D) * 1000 / vx_cov]
summary_s2[, ce := (s2_cost - s2c_cost)/ (pop_size * (s2c_D - s2_D))]
summary_s2[, bc := pop_size * (s2c_D - s2_D) * vsl_1 / (s2_cost - s2c_cost)]
summary_s2
summary_s2[, unique(covid_delay)]
write.csv(summary_s2,"summary_s2.csv")