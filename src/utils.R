libraries = c("dplyr","magrittr","tidyr","reshape2","ggplot2","openxlsx","RColorBrewer","zoo",
              "readxl","writexl","gridExtra","MASS","readr","stats","pracma","stringr","mixdist","corpcor","rstan","bayesplot")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE) }

setwd("~/covid-19_cfr_ve_tokyo")

# 1. Import data ####

## linelist of cases ####

df_case <- read_excel("./data/data_case.xlsx") %>% mutate_at(vars(date_confirm),~as.Date(.)) 

## linelist of deaths ####

df_death <- read_excel("./data/data_death.xlsx") %>% mutate_at(vars(date_confirm,date_death),~as.Date(.)) 

## daily number of cases by vac status by epi-week ####

df_case_vac <- read_excel("./data/data_casenum_byvac.xlsx") %>% mutate_at(vars(date_confirm),~as.Date(.)) 

## daily number of deaths by vac status by epi-week ####

df_death_vac <- read_excel("./data/data_deathnum_byvac.xlsx") %>% mutate_at(vars(date),~as.Date(.)) 

## daily number of cases and proportion of fully vaccinated (7 days moving average) ####

df_case_vac_prop <- read_excel("./data/data_casenum_vacprop.xlsx") %>% mutate_at(vars(date_confirm),~as.Date(.)) 

## daily number of deaths and proportion of fully vaccinated (7 days moving average) ####

df_death_vac_prop <- read_excel("./data/data_deathnum_vacprop.xlsx") %>% mutate_at(vars(date_death),~as.Date(.)) 

## health care burden(hcb) ####

df_hcb <- read_excel("./data/data_healthcare-burden.xlsx") %>% mutate_at(vars(date),~as.Date(.)) 

## daily proportion of asymptomatic cases by age (delay adjusted) ####

df_asymp_prop <- read_excel("./data/data_asymp_prop.xlsx") %>% mutate_at(vars(date),~as.Date(.)) 

# 2. Estimation of the distribution of the delay from the date of diagnosis and the date of death by time period and age ####
## The code was written by Yura Ko K

df_death_before_vac <- df_death %>% 
  filter(date_confirm >= "2021-01-01" & date_confirm <= "2021-04-11")

df_death_after_vac <- df_death %>% 
  filter(date_confirm >= "2021-04-12")

log_likelihood_Weibull <- function(parameters, observed_data){
  a <- parameters[1]
  b <- parameters[2]
  n <- length(x = observed_data)
  return((n * log(x = a)) + ((a - 1) * sum(log(x = observed_data))) - (a * n * log(x = b)) - (sum(observed_data ^ a) / (b ^ a)))
}

get_dist_diag_death <- function(data_,age_){
  dist_diag_death <- data_ %>% 
    filter(age == age_) %>% 
    filter(!is.na(date_confirm)&!is.na(date_death)) %>% 
    mutate(delay = as.numeric(date_death - date_confirm)) %>% 
    filter(delay >= 0) %>% 
    mutate(delay = ifelse(delay <= 0, 0.01, delay)) %>% 
    pull(delay)
  
  optim_diag_death <- optim(par = c(1, 1),
                            fn = log_likelihood_Weibull,
                            observed_data = dist_diag_death,
                            control = list(fnscale = -1))
  
  get_diag_death <- function(t){pweibull(t,  shape=optim_diag_death[[1]][1], scale=optim_diag_death[[1]][2]) - 
      pweibull(t-1, shape=optim_diag_death[[1]][1], scale=optim_diag_death[[1]][2])}
  
  par <- optim_diag_death$par
  
  plot <- ggplot(data.frame(x = seq(0, 60, length.out = 60)), aes(x))+
    scale_x_continuous(breaks = seq(0,60,5))+
    xlab("days from confirmation to death")+
    ylab("probability")+
    stat_function(fun = dweibull, args = list(shape=par[[1]], scale=par[[2]]))+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=12, family="sans",color="black"),
          legend.position="top",
          legend.title=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  return(list(par,plot))
}

dist_diag_deaath_par_by_age <- tibble(age = NA, phase = NA, shape = NA, scale = NA)
list_age <- unique(df_death$age)
for (age_ in list_age){
  dist_before <- get_dist_diag_death(df_death_before_vac,age_)
  dist_after <- get_dist_diag_death(df_death_after_vac,age_)
  dsit_res <- tibble(age = c(age_,age_),
                     phase = c("before","after"),
                     shape = c(dist_before[[1]][1],dist_after[[1]][1]),
                     scale = c(dist_before[[1]][2],dist_after[[1]][2]))
  dist_diag_deaath_par_by_age <- rbind(dist_diag_deaath_par_by_age,dsit_res) %>% 
    filter(!is.na(age))
}

# 3. Define function for convolution of the incidence on a confirmed date with the function representing the relative frequency of the delay ####
## The code was written by Yura Ko K

get_unbiased_case <- function(n,shape_,scale_){
  get_dist <- function(t){pweibull(t,  shape=shape_, scale=scale_) - 
      pweibull(t-1, shape=shape_, scale=scale_)}
  
  df <- tibble(t_ = NA, unbiased_case = NA)
  for (t in seq(1:60)){
    df_t <- tibble(t_ = t, unbiased_case = n*get_dist(t))
    df <- rbind(df,df_t) 
  }
  df_ <- df %>% slice(-1)
  return(df_)
}
get_casenum_adjusted_delay <- function(data,age_group,shape_,scale_){
  
  df_case_ <- data %>% 
    filter(age==age_group) %>% 
    mutate(t=row_number()) 
  
  list_date <- df_case_$total
  
  df_base <- tibble(t_ = NA,unbiased_case = NA)
  for (i in seq(1:nrow(df_case_))) {
    df_ <- get_unbiased_case(list_date[[i]],shape_,scale_) %>% 
      mutate(t_ = t_ + i - 1) 
    df_base <- rbind(df_base,df_)
  }
  
  df_unbiased_case <- df_base %>% 
    filter(!is.na(t_)) %>% 
    group_by(t_) %>% 
    mutate(sum_n = sum(unbiased_case)) %>% 
    distinct(t_,.keep_all = T) %>% 
    mutate(t = t_ +1) %>% 
    ungroup()
  
  df_case_adjusted <- left_join(df_unbiased_case,df_case_,by="t") %>% 
    dplyr::select(t,date,age,sum_n) 
  
  return(df_case_adjusted)
}
