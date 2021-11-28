setwd("~/covid19_cfr_ve_tokyo")
source("./src/utils.R")

# 1. Convolution of the incidence on a confirmed date with the function representing the relative frequency of the delay for each vaccine status####
## The code was written by Yura Ko K

get_case_period <- function(vac_status){
  
  df_case_vac_ <- df_case_vac %>% 
    rename(total = vac_status, date = "date_confirm") %>% 
    dplyr::select(date,age,total) 
  df_case_vac_before_vac <- df_case_vac_ %>% 
    filter(date <= "2021-04-11") %>% 
    group_by(age) %>% 
    complete(date = seq.Date(as.Date("2020-10-01"), as.Date("2021-06-10"), by="day")) %>% 
    replace_na(., replace = list(total= 0)) %>% 
    ungroup()
  df_case_vac_after_vac <- df_case_vac_ %>% 
    filter(date >= "2021-04-12")
  return(list(df_case_vac_before_vac,df_case_vac_after_vac))
}
confirm_u_before_vac <- get_case_period("not_vaccinated")[[1]]
confirm_u_after_vac <- get_case_period("not_vaccinated")[[2]]
confirm_v1_before_vac <- get_case_period("partial_vaccinated")[[1]]
confirm_v1_after_vac <- get_case_period("partial_vaccinated")[[2]]
confirm_v2_before_vac <- get_case_period("fully_vaccinated")[[1]]
confirm_v2_after_vac <- get_case_period("fully_vaccinated")[[2]]

df_case_vac_adjusted_delay <- tibble(t=NA,date=as.Date("2020-01-01"),age=NA,u=NA,v_1=NA,v_2=NA,phase=NA)
for (age_ in list_age){
  #before
  shape_before <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "before"))[[3]][1]
  scale_before <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "before"))[[4]][1]
  case_vac_adjusted_delay_u_before <- get_casenum_adjusted_delay(confirm_u_before_vac,age_,shape_before,scale_before) %>% 
    rename(u = "sum_n") %>% 
    mutate(phase = "before")
  
  case_vac_adjusted_delay_v_1_before <- get_casenum_adjusted_delay(confirm_v1_before_vac,age_,shape_before,scale_before) %>% 
    rename(v_1 = "sum_n") %>% 
    mutate(phase = "before")
  
  case_vac_adjusted_delay_v_2_before <- get_casenum_adjusted_delay(confirm_v2_before_vac,age_,shape_before,scale_before) %>% 
    rename(v_2 = "sum_n") %>% 
    mutate(phase = "before")
  
  count_vac_adjusted_delay_before <- left_join(case_vac_adjusted_delay_u_before,case_vac_adjusted_delay_v_1_before,by=c("t","date","age","phase")) %>% 
    left_join(case_vac_adjusted_delay_v_2_before,by=c("t","date","age","phase")) %>% 
    filter(date >= "2020-12-29" & date <= "2021-09-03") 
  
  #after
  shape_after <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "after"))[[3]][1]
  scale_after <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "after"))[[4]][1]
  case_vac_adjusted_delay_u_after <- get_casenum_adjusted_delay(confirm_u_after_vac,age_,shape_after,scale_after) %>% 
    rename(u = "sum_n") %>% 
    mutate(phase = "after")
  
  case_vac_adjusted_delay_v_1_after <- get_casenum_adjusted_delay(confirm_v1_after_vac,age_,shape_after,scale_after) %>% 
    rename(v_1 = "sum_n") %>% 
    mutate(phase = "after")
  
  case_vac_adjusted_delay_v_2_after <- get_casenum_adjusted_delay(confirm_v2_after_vac,age_,shape_after,scale_after) %>% 
    rename(v_2 = "sum_n") %>% 
    mutate(phase = "after")
  
  count_vac_adjusted_delay_after <- left_join(case_vac_adjusted_delay_u_after,case_vac_adjusted_delay_v_1_after,by=c("t","date","age","phase")) %>% 
    left_join(case_vac_adjusted_delay_v_2_after,by=c("t","date","age","phase")) %>% 
    filter(date >= "2020-12-29" & date <= "2021-09-03") 
  
  df_case_vac_adjusted_delay <- rbind(df_case_vac_adjusted_delay,count_vac_adjusted_delay_before,count_vac_adjusted_delay_after)
}

df_case_vac_adjusted_delay <- df_case_vac_adjusted_delay %>% 
  rename(not_vaccinated = "u",partial_vaccinated = "v_1",fully_vaccinated = "v_2") %>% 
  filter(!is.na(t)) %>% 
  group_by(date,age) %>% 
  mutate(not_vaccinated = sum(not_vaccinated),
         partial_vaccinated = sum(partial_vaccinated),
         fully_vaccinated = sum(fully_vaccinated)) %>% 
  ungroup() %>% 
  distinct(date,age,.keep_all = T) %>% 
  dplyr::select(date,age,not_vaccinated,partial_vaccinated,fully_vaccinated)

# write_xlsx(df_case_vac_adjusted_delay,"./output/data_casenum_vac_adjusted_delay.xlsx")
# df_case_vac_adjusted_delay <- read_excel("./output/data_casenum_vac_adjusted_delay.xlsx") %>% mutate_at(vars(date),~as.Date(.))

# 2. Joint estimation of time-varying CFR and VE from death by age group ####
## The code was written by Hiroaki Murayama

case <- df_case_vac_adjusted_delay %>%  
  mutate(unvac = rollmean(not_vaccinated, 7,fill = NA),
         vac1 = rollmean(partial_vaccinated, 7,fill = NA),
         vac2 = rollmean(fully_vaccinated, 7,fill = NA)) %>% 
  filter(date >= "2021-01-01" & date <= "2021-08-31") %>% 
  dplyr::select(date,age,unvac,vac1,vac2)

death <- df_death_vac %>%
  mutate(unvac = rollmean(not_vaccinated, 7,fill = NA),
         vac1 = rollmean(partial_vaccinated, 7,fill = NA),
         vac2 = rollmean(fully_vaccinated, 7,fill = NA)) %>% 
  filter(date >= "2021-01-01" & date <= "2021-08-31") %>% 
  dplyr::select(date,age,unvac,vac1,vac2) 

case30_50 <- case %>% filter(age=="30_50s") %>% as.data.frame()
case60 <- case %>% filter(age=="60s") %>% as.data.frame()
case70 <- case %>% filter(age=="70s") %>% as.data.frame()
case80 <- case %>% filter(age=="80s") %>% as.data.frame()
case90_100 <- case %>% filter(age=="90_100s") %>% as.data.frame()

death30_50 <- death %>% filter(age=="30_50s") %>% as.data.frame()
death60 <- death %>% filter(age=="60s") %>% as.data.frame()
death70 <- death %>% filter(age=="70s") %>% as.data.frame()
death80 <- death %>% filter(age=="80s") %>% as.data.frame()
death90_100 <- death %>% filter(age=="90_100s") %>% as.data.frame()

# stan model
Model_vac2 = "
functions{
real continuous_binomial_lpdf(real d, real theta, real r) {
    real lp = lgamma(theta + 1) - lgamma(d+1) - lgamma(theta-d+1) + d*log(r) + (theta-d)*log(1-r);
    return lp;
  }
}

data {
int T;
real Cn[T];
real C1v[T];
real C2v[T];
real Dn[T];
real D1v[T];
real D2v[T];
real vei1;
real vei2;
}

parameters {
real<lower=0, upper=1> ved1;
real<lower=0, upper=1> ved2;
real<lower=0, upper=1> p[T];
}

model{  
ved1 ~ beta(1,1);
ved2 ~ beta(1,1);
p ~ beta(1,1);

for(tt in 1:T){
   target += continuous_binomial_lpdf(Dn[tt] | Cn[tt], p[tt]) + continuous_binomial_lpdf(D1v[tt] | C1v[tt], (1-ved1)*p[tt]) 
             + continuous_binomial_lpdf(D2v[tt] | C2v[tt], (1-ved2)*p[tt]);
   }
}
generated quantities{
real<lower=0> eps1;
real<lower=0> eps2;

eps1 = vei1 +(1-vei1) * ved1;
eps2 = vei2 +(1-vei2) * ved2;
}
"

get_ve_cfr_res <- function(case_,death_,age_,vei1_,vei2_){
  
  T = nrow(case_) # number of days observed 
  Cn = case_$unvac # daily number of un vaccinated cases 
  C1v = case_$vac1 # daily number of vaccinated cases with the 1st dose 
  C2v = case_$vac2 # daily number of vaccinated cases with the 2nd dose
  Dn = death_$unvac # daily number of unvaccinated deaths 
  D1v = death_$vac1 # daily number of vaccinated deaths
  D2v = death_$vac2 # daily number of vaccinated deaths
  vei1 = vei1_ # ve1 against documented infection (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577)
  vei2 = vei2_ # ve2 against documented infection (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577)
  
  data = list(T=T, Cn=Cn, C1v=C1v, C2v=C2v, Dn=Dn, D1v=D1v, D2v=D2v, vei1=vei1, vei2=vei2)
  # specify parameters to monitor
  parameters = c("ved1","ved2","eps1","eps2","p")
  nuts_fit = stan(model_code=Model_vac2, data=data, pars=parameters, iter=10000, thin=10, warmup=1000, chain=4)
  
  options(repr.plot.width=10,repr.plot.height=4)
  mcmc_trace(nuts_fit,  pars = c("ved1","ved2","eps1","eps2"), size=0.6)+ facet_text(size = 15)
  
  ### extract the mcmc results 
  result_median <- summary(nuts_fit)$summary[, "50%"]
  result_median <- as.data.frame(result_median)
  
  result_2.5 <- summary(nuts_fit)$summary[, "2.5%"]
  result_2.5 <- as.data.frame(result_2.5)
  
  result_97.5 <- summary(nuts_fit)$summary[, "97.5%"]
  result_97.5 <- as.data.frame(result_97.5) 
  
  result_all <- cbind(result_median,result_2.5,result_97.5) %>% 
    mutate(age = age_)
  
  var <- tibble(var = c("ve1","ve2","eps1","eps2"))
  res_ve <- result_all %>% 
      slice(1:4) %>% 
      cbind(var) %>% 
      rename(ve_med = "result_median", ve_2.5 = "result_2.5", ve_97.5 = "result_97.5") %>% 
      dplyr::select(age,var,ve_med,ve_2.5,ve_97.5)
  
  tail <- nrow(result_all)-1
  res_cfr <- result_all %>% 
    slice(5:tail) %>% 
    mutate(date = seq(as.Date("2021-01-01"), as.Date("2021-08-31"), by = "day"))
  
  return(list(res_ve,res_cfr))
}

ve_cfr_res_30_50 <- get_ve_cfr_res(case30_50,death30_50,"30_50s",0.38,0.90)
ve_cfr_res_60 <- get_ve_cfr_res(case60,death60,"60s",0.812,0.926)
ve_cfr_res_70 <- get_ve_cfr_res(case70,death70,"70s",0.764,0.954)
ve_cfr_res_80 <- get_ve_cfr_res(case80,death80,"80s",0.803,0.961)
ve_cfr_res_90_100 <- get_ve_cfr_res(case90_100,death90_100,"90_100s",0.757,0.926)

res_ve <- rbind(ve_cfr_res_30_50[[1]],ve_cfr_res_60[[1]],
                ve_cfr_res_70[[1]],ve_cfr_res_80[[1]],
                ve_cfr_res_90_100[[1]]) %>% 
  arrange(var)
write_xlsx(res_ve,"./output/res_ve.xlsx")

res_cfr <- rbind(ve_cfr_res_30_50[[2]],ve_cfr_res_60[[2]],
                ve_cfr_res_70[[2]],ve_cfr_res_80[[2]],
                ve_cfr_res_90_100[[2]]) %>% 
  dplyr::select(date,age,result_median,result_2.5,result_97.5)
write_xlsx(res_cfr,"./output/res_jointly-estimated-cfr.xlsx")

# 3. The impact of indicators of the level of healthcare burden on CFR by age group ####
## The code was written by Hiroaki Murayama

asymp30_50 <- df_asymp_prop %>% filter(age=="30_50s") 
asymp60 <- df_asymp_prop %>% filter(age=="60s") 
asymp70 <- df_asymp_prop %>% filter(age=="70s") 
asymp80 <- df_asymp_prop %>% filter(age=="80s") 
asymp90_100 <- df_asymp_prop %>% filter(age=="90_100s") 

# stan code 
Model = "
// modified binomial distr.
functions{
real continuous_binomial_lpdf(real d, real theta, real r) {
    real lp = lgamma(theta + 1) - lgamma(d+1) - lgamma(theta-d+1) + d*log(r) + (theta-d)*log(1-r);
    return lp;
  }


// data related process
real f(real C, real D, real p) {
    return(continuous_binomial_lpdf(D | C, p));
  }

real log_lik_Simpson_1(real C, real D, real a, real b, int M) {
    vector[M+1] lp;
    real h;
    h = (b-a)/M;
    lp[1] = f(C, D, a);
    for (m in 1:(M/2))
      lp[2*m] = log(4) + f(C, D, a + h*(2*m-1));
    for (m in 1:(M/2-1))
      lp[2*m+1] = log(2) + f(C, D, a + h*2*m);
    lp[M+1] = f(C, D, b);
    return(log(h/3) + log_sum_exp(lp));
  }


// hyperparameter related integration
real f_2(real x1, real x2, real a1, real a2, real b, real kappa, real p) {
    return(beta_proportion_lpdf(p | inv_logit(a1*x1 + a2*x2 + b), kappa));
  }

real log_lik_Simpson_2(real x1, real x2, real a1, real a2, real b, real kappa, real p_lower, real p_upper, int M) {
    vector[M+1] lp;
    real h;
    h = (p_upper-p_lower)/M;
    lp[1] = f_2(x1, x2, a1, a2, b, kappa, p_lower);
    for (m in 1:(M/2))
      lp[2*m] = log(4) + f_2(x1, x2, a1, a2, b, kappa, p_lower + h*(2*m-1));
    for (m in 1:(M/2-1))
      lp[2*m+1] = log(2) + f_2(x1, x2, a1, a2, b, kappa, p_lower + h*2*m);
    lp[M+1] = f_2(x1, x2, a1, a2, b, kappa, p_upper);
    return(log(h/3) + log_sum_exp(lp));
  }
}

data {
int T;
real C[T];
real D[T];
vector[T] x1;
vector[T] x2;
}

parameters {
real<lower=0, upper=1> p[T];
real a1;
real a2;
real b;
real<lower=0> kappa;
}

transformed parameters{
real<lower=0> mu[T];
for(i in 1:T){
mu[i] = inv_logit(a1*x1[i] + a2*x2[i] + b);
}
}
model{  
a1 ~ normal(0,100);
a2 ~ normal(0,100);
b ~ normal(0,100);
kappa ~ normal(0,10); 

for(n in 1:T){
  target += continuous_binomial_lpdf(D[n] | C[n], p[n]) + beta_proportion_lpdf(p[n] | inv_logit(a1*x1[n] + a2*x2[n] + b), kappa);
   }
}

generated quantities{
real log_lik[T];
for(n in 1:T)
log_lik[n] = log_lik_Simpson_1(C[n], D[n], 0.001, 1, 200) + log_lik_Simpson_2(x1[n], x2[n], a1, a2, b, kappa, 0.001, 1, 200);
}
"

get_cfr_regression <- function(data_c,data_d,data_asymp,var_name,age_){
  # modify data into a form suitable for Stan
  
  x <- df_hcb %>% 
    rename(var = var_name) 
  T = nrow(data_c)
  C = data_c$unvac 
  D = data_d$unvac
  x1 = x$var
  x2 = data_asymp$prop_asymp # asymptomatic rate
  data = list(T=T, C=C, D=D, x1=x1, x2=x2) 
  # specify parameters to monitor
  parameters = c("p","a1","a2","b","kappa","log_lik")
  nuts_fit = stan(model_code=Model,data=data,pars=parameters,iter=10000,thin=10,warmup=1000,chain=4)
  
  get_waic <- function(nuts_fit){
    # compute WAIC from the returned object from Stan
    # the log likelihood must be named as 'log_lik'
    waic <- function(log_lik) {
      Tn <- - mean(log(colMeans(exp(log_lik))))
      fV <- mean(colMeans(log_lik^2) - colMeans(log_lik)^2)
      waic <- Tn + fV
      waic
    }
    
    nuts_fit %>% rstan::extract() %>% .$log_lik %>% waic()
  }
  waic = get_waic(nuts_fit)
  
  nuts_fit
  #summary(nuts_fit)
  options(repr.plot.width=10,repr.plot.height=5)
  color_scheme_set("mix-blue-pink")
  mcmc_trace(nuts_fit,  pars = c("a1","a2","b","kappa"), size=0.6)+ facet_text(size = 15)
  
  ### extract the mcmc results 
  result_median <- summary(nuts_fit)$summary[, "50%"]
  result_median <- as.data.frame(result_median)
  
  median_a<-  list()
  for(i in 1:T){
    name_a <- paste("a_raw[", i, "]", sep = "")
    median_a$a_raw[i] <- result_median [name_a ,]
    remove(name_a)
  }
  
  result_2.5 <- summary(nuts_fit)$summary[, "2.5%"]
  result_2.5 <- as.data.frame(result_2.5)
  
  lower_a<-  list()
  for(i in 1:T){
    name_a <- paste("a_raw[", i, "]", sep = "")
    lower_a$a_raw[i] <- result_2.5 [name_a ,]
    remove(name_a)
  }
  
  result_97.5 <- summary(nuts_fit)$summary[, "97.5%"]
  result_97.5 <- as.data.frame(result_97.5) 
  
  upper_a<-  list()
  for(i in 1:T){
    name_a <- paste("a_raw[", i, "]", sep = "")
    upper_a$a_raw[i] <- result_2.5 [name_a ,]
    remove(name_a)
  }
  
  result <- cbind(result_median,result_2.5,result_97.5) 
  result_a1 <- result[T+1,]
  result <- result[1:T,] %>% 
    mutate(date = seq(as.Date("2021-01-01"), as.Date("2021-08-31"), by = "day"))
  
  colnames(result) <- c("median","lower","upper","date") 
  head(result);tail(result)
  
  result_table <- tibble(Age = age_, Model = var_name) %>% 
    cbind(result_a1) %>% 
    mutate(waic = waic)
  
  write_xlsx(result,paste("./output/result_",age_,"_",var_name,".xlsx",sep = ""))
  return(list(result,result_table))
}
get_cfr_regression_res <- function(data_c,data_d,data_asymp,age_){
  res_cfr_regression_severe_lag1 <- get_cfr_regression(data_c,data_d,data_asymp,"severe_lag1",age_)
  res_cfr_regression_severe_lag3 <- get_cfr_regression(data_c,data_d,data_asymp,"severe_lag3",age_)
  res_cfr_regression_severe_lag5 <- get_cfr_regression(data_c,data_d,data_asymp,"severe_lag5",age_)
  res_cfr_regression_tokyo_rule_lag1 <- get_cfr_regression(data_c,data_d,data_asymp,"tokyo_rule_lag1",age_)
  res_cfr_regression_tokyo_rule_lag3 <- get_cfr_regression(data_c,data_d,data_asymp,"tokyo_rule_lag3",age_)
  res_cfr_regression_tokyo_rule_lag5 <- get_cfr_regression(data_c,data_d,data_asymp,"tokyo_rule_lag5",age_)
  res_cfr_regression_non_hospi_prop_lag1 <- get_cfr_regression(data_c,data_d,data_asymp,"non_hospi_prop_lag1",age_)
  res_cfr_regression_non_hospi_prop_lag3 <- get_cfr_regression(data_c,data_d,data_asymp,"non_hospi_prop_lag3",age_)
  res_cfr_regression_non_hospi_prop_lag5 <- get_cfr_regression(data_c,data_d,data_asymp,"non_hospi_prop_lag5",age_)
  res_cfr_regression_adjus_prop_lag1 <- get_cfr_regression(data_c,data_d,data_asymp,"adjus_prop_lag1",age_)
  res_cfr_regression_adjus_prop_lag3 <- get_cfr_regression(data_c,data_d,data_asymp,"adjus_prop_lag3",age_)
  res_cfr_regression_adjus_prop_lag5 <- get_cfr_regression(data_c,data_d,data_asymp,"adjus_prop_lag5",age_)
  
  res_cfr_regression <- rbind(res_cfr_regression_severe_lag1[[2]],res_cfr_regression_severe_lag3[[2]],res_cfr_regression_severe_lag5[[2]],
                              res_cfr_regression_tokyo_rule_lag1[[2]],res_cfr_regression_tokyo_rule_lag3[[2]],res_cfr_regression_tokyo_rule_lag5[[2]],
                              res_cfr_regression_non_hospi_prop_lag1[[2]],res_cfr_regression_non_hospi_prop_lag3[[2]],res_cfr_regression_non_hospi_prop_lag5[[2]],
                              res_cfr_regression_adjus_prop_lag1[[2]],res_cfr_regression_adjus_prop_lag3[[2]],res_cfr_regression_adjus_prop_lag5[[2]])
  
  return(res_cfr_regression)
}

cfr_regression_all_30_50s <- get_cfr_regression_res(inf30_50,death30_50,asymp30_50,"30_50s")
write_xlsx(cfr_regression_all_30_50s,"./output/cfr_regression_all_30_50s.xlsx")
cfr_regression_all_30_50s <- read_xlsx("./output/cfr_regression_all_30_50s.xlsx")

cfr_regression_all_60s <- get_cfr_regression_res(inf60,death60,asymp60,"60s")
write_xlsx(cfr_regression_all_60s,"./output/cfr_regression_all_60s.xlsx")
cfr_regression_all_60s <- read_xlsx("./output/cfr_regression_all_60s.xlsx")

cfr_regression_all_70s <- get_cfr_regression_res(inf70,death70,asymp70,"70s")
write_xlsx(cfr_regression_all_70s,"./output/cfr_regression_all_70s.xlsx")
cfr_regression_all_70s <- read_xlsx("./output/cfr_regression_all_70s.xlsx")

cfr_regression_all_80s <- get_cfr_regression_res(inf80,death80,asymp80,"80s")
write_xlsx(cfr_regression_all_80s,"./output/cfr_regression_all_80s.xlsx")
cfr_regression_all_80s <- read_xlsx("./output/cfr_regression_all_80s.xlsx")

cfr_regression_all_90_100s <- get_cfr_regression_res(inf90_100,death90_100,asymp90_100,"90_100s")
write_xlsx(cfr_regression_all_90_100s,"./output/cfr_regression_all_90s.xlsx")
cfr_regression_all_90_100s <- read_xlsx("./output/cfr_regression_all_90s.xlsx")

cfr_regression_all <- rbind(cfr_regression_all_30_50s,cfr_regression_all_60s,
                            cfr_regression_all_70s,cfr_regression_all_80s,
                            cfr_regression_all_90_100s) %>% 
  mutate(coefficient = paste(sprintf(result_median,fmt= "%0.3f")," (",sprintf(result_2.5,fmt= "%0.3f"),"–",sprintf(result_97.5,fmt= "%0.3f"),")",sep = "")) %>% 
  mutate(group_2 = substr(Model,1,5)) %>% 
  group_by(Age,group_2) %>% 
  filter(waic == min(waic)) %>% 
  dplyr::select(Age,Model,coefficient,waic) 

write_xlsx(cfr_regression_all,"./output/cfr_regression_all.xlsx")


