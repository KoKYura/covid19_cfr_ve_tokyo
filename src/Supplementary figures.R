setwd("~/covid-19_cfr_ve_tokyo")
source("./src/utils.R")

#############################################################################################################                                                                          #
#                                                                                                           #
## The code was written by Yura Ko K                                                                        #
#############################################################################################################

# 1. Import and processing data ####

## daily number of delay adjusted cases ####

df_case_count <- df_case %>% 
  group_by(age,date_confirm) %>% 
  count() %>% 
  group_by(age) %>% 
  complete(date_confirm = seq.Date(as.Date("2020-10-01"), as.Date("2021-09-03"), by="day")) %>% 
  replace_na(., replace = list(n = 0)) %>% 
  filter(date_confirm <= "2021-09-03") %>% 
  rename(date = "date_confirm" ,total = "n")

df_case_count_before <- df_case_count %>% 
  filter(date <= "2021-04-11") %>% 
  group_by(age) %>% 
  complete(date = seq.Date(as.Date("2020-10-01"), as.Date("2021-06-10"), by="day")) %>% 
  replace_na(., replace = list(total= 0)) %>% 
  ungroup()
df_case_count_after <- df_case_count %>% 
  filter(date >= "2021-04-12")

df_case_adjusted_count <- tibble(t=NA,date=as.Date("2020-01-01"),age=NA,case=NA,phase=NA)
for (age_ in list_age){
  shape_before <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "before"))[[3]][1]
  scale_before <- (dist_diag_deaath_par_by_age %>% filter(age == age_) %>% filter(phase == "before"))[[4]][1]
  
  df_case_adjusted_count_before <- get_casenum_adjusted_delay(df_case_count_before,age_,shape_before,scale_before) %>% 
    rename(case = "sum_n") %>% 
    mutate(phase = "before") %>% 
    filter(date >= "2020-12-29" & date <= "2021-09-03") 
  df_case_adjusted_count_after <- get_casenum_adjusted_delay(df_case_count_after,age_,shape_before,scale_before) %>% 
    rename(case = "sum_n") %>% 
    mutate(phase = "after") %>% 
    filter(date >= "2020-12-29" & date <= "2021-09-03") 
  
  df_case_adjusted_count <- rbind(df_case_adjusted_count,df_case_adjusted_count_before,df_case_adjusted_count_after)
}

df_case_adjusted_count <- df_case_adjusted_count %>% 
  filter(!is.na(age)) %>% 
  group_by(date,age) %>% 
  mutate(sum_case = sum(case)) %>% 
  ungroup() %>% 
  distinct(date,age,.keep_all = T) %>% 
  filter(date >= "2020-12-29" & date <= "2021-09-03") %>% 
  dplyr::select(age,date,sum_case) 

## daily number of deaths ####

df_death_count <- df_death %>% 
  group_by(age,date_death) %>% 
  count() %>% 
  group_by(age) %>% 
  complete(date_death = seq.Date(as.Date("2020-12-29"), as.Date("2021-09-03"), by="day")) %>% 
  replace_na(., replace = list(n = 0)) %>% 
  filter(date_death >= "2020-12-29" & date_death <= "2021-09-03") %>% 
  rename(date = "date_death", sum_death = "n")

## Jointly estimated unvaccinated CFR ####

df_unvaccinated_cfr <- read_excel("./output/res_jointly-estimated-cfr.xlsx") %>% mutate(date = as.Date(date))

# 2. Supplementary Figure 1.####
## The distribution of days from confirmation to death fitted to Weibull distribution by age group and time.

get_dist_comp <- function(data_1,data_2,age_){
  
  dist_diag_death_1 <- data_1 %>% 
    filter(age == age_) %>% 
    filter(!is.na(date_confirm)&!is.na(date_death)) %>% 
    mutate(delay = as.numeric(date_death - date_confirm)) %>% 
    filter(delay >= 0) %>% 
    mutate(delay = ifelse(delay <= 0, 0.01, delay)) %>% 
    pull(delay)
  
  dist_diag_death_2 <- data_2 %>% 
    filter(age == age_) %>% 
    filter(!is.na(date_confirm)&!is.na(date_death)) %>% 
    mutate(delay = as.numeric(date_death - date_confirm)) %>% 
    filter(delay >= 0) %>% 
    mutate(delay = ifelse(delay <= 0, 0.01, delay)) %>% 
    pull(delay)
  
  optim_diag_death_1 <- optim(par = c(1, 1),
                              fn = log_likelihood_Weibull,
                              observed_data = dist_diag_death_1,
                              control = list(fnscale = -1))
  optim_diag_death_2 <- optim(par = c(1, 1),
                              fn = log_likelihood_Weibull,
                              observed_data = dist_diag_death_2,
                              control = list(fnscale = -1))
  
  par_1 <- optim_diag_death_1$par
  par_2 <- optim_diag_death_2$par
  
  iqr_l_1 <- qweibull(0.25,shape=par_1[[1]], scale=par_1[[2]])
  iqr_u_1 <- qweibull(0.75,shape=par_1[[1]], scale=par_1[[2]])
  
  dweibull_limit_1 <- function(x){
    y <- dweibull(x,shape=par_1[[1]], scale=par_1[[2]])
    y[x < iqr_l_1 | x > iqr_u_1] <- NA
    return(y)
  }
  
  iqr_l_2 <- qweibull(0.25,shape=par_2[[1]], scale=par_2[[2]])
  iqr_u_2 <- qweibull(0.75,shape=par_2[[1]], scale=par_2[[2]])
  
  dweibull_limit_2 <- function(x){
    y <- dweibull(x,shape=par_2[[1]], scale=par_2[[2]])
    y[x < iqr_l_2 | x > iqr_u_2] <- NA
    return(y)
  }
  
  plot <- ggplot(data.frame(x = seq(0, 60, length.out = 60)), aes(x))+
    scale_x_continuous(breaks = seq(0,60,5))+
    xlab("Days from confirmation to death")+
    ylab("Probability")+
    stat_function(fun = dweibull, args = list(shape=par_1[[1]], scale=par_1[[2]]),colour = "#d8b365",size = 1.5)+
    stat_function(fun = dweibull_limit_1, geom= "area",fill="#d8b365",alpha=0.2)+
    stat_function(fun = dweibull, args = list(shape=par_2[[1]], scale=par_2[[2]]),colour = "#5ab4ac",size = 1.5)+
    stat_function(fun = dweibull_limit_2, geom= "area",fill="#5ab4ac",alpha=0.2)+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(axis.text =element_text(size=18, family="sans",color="black"),
          title = element_text(size=20, family="sans",color="black"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(plot)
}

dist_comp_30_50 <- get_dist_comp(df_death_before_vac,df_death_after_vac,"30_50s")
dist_comp_60 <- get_dist_comp(df_death_before_vac,df_death_after_vac,"60s")
dist_comp_70 <- get_dist_comp(df_death_before_vac,df_death_after_vac,"70s")
dist_comp_80 <- get_dist_comp(df_death_before_vac,df_death_after_vac,"80s")
dist_comp_90_ <- get_dist_comp(df_death_before_vac,df_death_after_vac,"90_100s")

g_dist_comp <-  plot_grid(dist_comp_30_50+ggtitle("30–50s"),
                          dist_comp_60+ggtitle("60s"),
                          dist_comp_70+ggtitle("70s"),
                          dist_comp_80+ggtitle("80s"),
                          dist_comp_90_+ggtitle("90s+"),
                          ncol=1,align = "v")

ggsave(g_dist_comp, filename = "./output/SuppFigure 1.png",width = 15,height = 18)

# 3. Supplementary Figure 2.####
## The estimated CFR including both vaccinated and unvaccinated cases by age group. 

df_case_death_count <- left_join(df_case_adjusted_count,df_death_count,by=c("date","age")) %>% 
  mutate(c = rollmean(sum_case, 7,fill = NA),
         d = rollmean(sum_death, 7,fill = NA)) %>% 
  filter(date >= "2021-01-01" & date <= "2021-08-31") 

#stan
Model = "
functions{
real continuous_binomial_lpdf(real d, real theta, real r) {
    real lp = lgamma(theta + 1) - lgamma(d+1) - lgamma(theta-d+1) + d*log(r) + (theta-d)*log(1-r);
    return lp;
  }
}

data {
int T;
real C[T];
real D[T];
}

parameters {
real<lower=0, upper=1> p[T];
}

model{  
p ~ beta(1,1);

for(tt in 1:T){
   target += continuous_binomial_lpdf(D[tt] | C[tt], p[tt]) ;
   }
}
"

get_cfr_all <- function(age_){
  
  # modify data into a form suitable for Stan
  data_ <- df_case_death_count %>% filter(age == age_)
  T = nrow(data_) # number of days observed 
  C = data_$c # daily number of un vaccinated cases 
  D = data_$d # daily number of unvaccinated deaths 
  
  data = list(T=T,C=C,D=D)
  # specify parameters to monitor
  parameters = c("p")
  nuts_fit = stan(model_code=Model,data=data,pars=parameters,iter=10000,thin=10,warmup=1000,chain=4)
  
  ### extract the mcmc results 
  result_median <- summary(nuts_fit)$summary[, "50%"]
  result_median <- as.data.frame(result_median)
  
  result_2.5 <- summary(nuts_fit)$summary[, "2.5%"]
  result_2.5 <- as.data.frame(result_2.5)
  
  result_97.5 <- summary(nuts_fit)$summary[, "97.5%"]
  result_97.5 <- as.data.frame(result_97.5) 
  
  
  result_ <- cbind(result_median,result_2.5,result_97.5) %>% 
    slice(1:243) %>% 
    mutate(date = seq(as.Date("2021-01-01"), as.Date("2021-08-31"), by = "day"))
  return(result_)
}

cfr_all_30_50s <- get_cfr_all("30_50s")
cfr_all_60s <- get_cfr_all("60s")
cfr_all_70s <- get_cfr_all("70s")
cfr_all_80s <- get_cfr_all("80s")
cfr_all_90s_ <- get_cfr_all("90_100s")

get_plot_cfr_all <- function(data,age_){
  g_result <- data %>% ggplot() +
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d"
    ) +
    geom_ribbon(aes(x = date, ymin = result_2.5*100, ymax = result_97.5*100), fill = "#EDBA13", alpha = 0.3) +
    geom_line(aes(x = date, y = result_median*100), color="#ED9521", size = 1) +
    labs(x="\nDate of death",y="Time delay adjusted CFR [%]\n") +
    ggtitle(age_) +
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=12, family="sans",color="black"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(g_result)
}

g_all_cfr <- grid.arrange(get_plot_cfr_all(cfr_all_30_50s,"30-50s"),
                          get_plot_cfr_all(cfr_all_60s,"60s"),
                          get_plot_cfr_all(cfr_all_70s,"70s"),
                          get_plot_cfr_all(cfr_all_80s,"80s"),
                          get_plot_cfr_all(cfr_all_90s_,"90s+"))
ggsave(plot = g_all_cfr,filename = "./output/SuppFigure 2.png",width = 12,height = 9)

# 4. Supplementary Figure 3.####
## The estimated CFR of only unvaccinated cases by age group. 

get_plot_cfr_unvaccinated <- function(age_,age_group){
  
  data_ <- df_unvaccinated_cfr %>% filter(age == age_)
  
  g_result <- data_ %>% ggplot() +
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d"
    ) +
    geom_ribbon(aes(x = date, ymin = result_2.5*100, ymax = result_97.5*100), fill = "#2c7fb8", alpha = 0.3) +
    geom_line(aes(x = date, y = result_median*100), color="#2c7fb8", size = 1) +
    labs(x="\nDate of death",y="Time delay adjusted CFR [%]\n") +
    ggtitle(age_group) +
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=12, family="sans",color="black"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(g_result)
}

g_all_cfr <- grid.arrange(get_plot_cfr_unvaccinated("30_50s","30-50s"),
                          get_plot_cfr_unvaccinated("60s","60s"),
                          get_plot_cfr_unvaccinated("70s","70s"),
                          get_plot_cfr_unvaccinated("80s","80s"),
                          get_plot_cfr_unvaccinated("90_100s","90s+"))
ggsave(plot = g_all_cfr,filename = "./output/SuppFigure 3.png",width = 12,height = 9)

# 5. Supplementary Figure 4.####
## The trajectories of indicators for health care burden in Tokyo. 

get_plot_hcb <- function(var,title){
  data_ <- df_hcb %>% 
    rename(var_ = var) 
  
  plot_hcb <- 
    ggplot(NULL) +
    ggtitle(title)+
    xlab("")+
    ylab("")+
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d",
      expand = c(0,0)
    ) +
    geom_line(data=data_,
              aes(x= date, 
                  y=var_
              ),
              stat="identity", 
              width=1
    )+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(axis.text =element_text(size=18, family="sans",color="black"),
          title = element_text(size=20, family="sans",color="black"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(plot_hcb)
}

plot_hcb <- plot_grid(get_plot_hcb("severe", "Number of severe cases"),
                           get_plot_hcb("tokyo_rule", "Number of Tokyo rule applied"),
                           get_plot_hcb("non_hospi_prop","Proportion of non-hospiralized cases"),
                           get_plot_hcb("adjus_prop","Proportion of cases in the process of cordinating the care site"),
                           ncol=1,align = "v")

ggsave(plot = plot_hcb,filename = "./output/SuppFigure 4.png",width = 15,height = 18)

# 6. Supplementary Figure 5.####
## The proportion of asymptomatic cases among diagnosed cases by age group.

plot_asymp_prop <-
  ggplot(NULL) +
  scale_color_manual(values = c('#66c2a5','#a6d854','#8da0cb','#e78ac3','#fc8d62'),labels = c("30–50s","60s","70s","80s","90s+")) +
  xlab("Date")+
  ylab("Propotion of asymptomatic cases (%)")+
  scale_x_date(
    date_breaks = "1 months",
    date_labels = "%b-%d"
  ) +
  geom_line(data = df_asymp_prop %>% 
              filter(date < "2021-09-01"),
            aes(x = date, 
                y = 100*prop_asymp,
                colour = age),
            stat = "identity",
            size = 0.7
  )+
  theme_bw() +
  theme(axis.text.x = element_text()) +
  theme(text=element_text(size=18, family="sans",color="black"),
        legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=18, family="sans",color="black"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot = plot_asymp_prop,filename = "./output/SuppFigure 5.png",width = 12,height = 9)



