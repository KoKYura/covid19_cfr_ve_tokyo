setwd("~/covid19_cfr_ve_tokyo")
source("./src/utils.R")

#############################################################################################################
# Cumulative number of deaths in counterfactual in which all positive individuals were vaccinated,          #
# in which all positive individuals were not vaccinated, and actual reported cases.                         #
#                                                                                                           #
## The code was written by Yura Ko K                                                                        #
#############################################################################################################

# 1. Import and processing data ####

## daily number of delay adjusted cases ####

df_case_count <- df_case %>% 
  group_by(age,date_confirm) %>% 
  count() %>% 
  group_by(age) %>% 
  complete(date_confirm = seq.Date(as.Date("2020-10-01"), as.Date("2021-08-31"), by="day")) %>% 
  replace_na(., replace = list(n = 0)) %>% 
  filter(date_confirm <= "2021-08-31") %>% 
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
  filter(date >= "2021-01-01" & date <= "2021-08-31") %>% 
  dplyr::select(age,date,sum_case) 

## daily number of deaths ####

df_death_count <- df_death %>% 
  group_by(age,date_death) %>% 
  count() %>% 
  group_by(age) %>% 
  complete(date_death = seq.Date(as.Date("2021-01-01"), as.Date("2021-08-31"), by="day")) %>% 
  replace_na(., replace = list(n = 0)) %>% 
  filter(date_death >= "2021-01-01" & date_death <= "2021-08-31") %>% 
  rename(date = "date_death", sum_death = "n")

## Estimated VE against death ####

df_ve <- read_excel("./output/res_ve.xlsx") 

## Jointly estimated unvaccinated CFR ####

df_unvaccinated_cfr <- read_excel("./output/res_jointly-estimated-cfr.xlsx") %>% mutate(date = as.Date(date))

# 2. Plot Figure 2  ####

get_plot <- function(age_group){
  
  df_ve_ <- df_ve %>% filter(age == age_group & var == "ve2")
  ve_d <- df_ve_ [[3]][1]
  ve_d_l <- df_ve_[[4]][1]
  ve_d_u <- df_ve_[[5]][1]
  
  data_combined <- df_unvaccinated_cfr %>% 
    left_join(df_case_adjusted_count,by=c("date","age")) %>% 
    left_join(df_death_count,by=c("date","age")) %>% 
    filter(age == age_group) %>% 
    mutate(d_u = sum_case * result_median,
           d_l_m = sum_case * result_median * (1-ve_d),
           d_l_u = sum_case * result_median * (1-ve_d_l),
           d_l_l = sum_case * result_median * (1-ve_d_u)) %>% 
    rename(d = "sum_death") %>% 
    dplyr::select(date,d,d_u,d_l_m,d_l_u,d_l_l) %>% 
    gather(2:6,key = "group", value = "death") %>% 
    group_by(group) %>% 
    mutate(cum_death = cumsum(death))
  
  unvac_cfr_est <- data_combined %>% 
    filter(group != "d_l_u" &group != "d_l_l")
  unvac_cfr_est$group <- fct_relevel(unvac_cfr_est$group,c("d","d_u","d_l_m"))
  
  unvac_cfr_ci <- data_combined %>% 
    filter(group == "d_l_u" |group == "d_l_l") %>% 
    dplyr::select(-death) %>% 
    spread(2:3,key = group,value = cum_death)
  
  g_death_comp <- 
    ggplot(NULL) +
    ggtitle("")+
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d"
    ) +
    scale_colour_manual(values = c("#d7191c","#2c7fb8","#f4a582"),labels = c("Actual","Counterfactual (all unvaccinated)","Counterfactual (all vaccinated)"))+
    xlab("Date of death")+
    ylab("Number of death")+
    geom_ribbon(data=unvac_cfr_ci, aes(x = date, ymin = d_l_l, ymax = d_l_u), fill = "#f4a582", alpha = 0.2) +
    geom_line(data=unvac_cfr_est,
              aes(x= date, 
                  y= cum_death,
                  colour = group
              ),
              stat="identity", 
    )+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=18, family="sans",color="black"),
          legend.position = "none",
          legend.title = element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(g_death_comp)
}

plot_all <- plot_grid(get_plot("30_50s")+ggtitle("30-50s"),
                            get_plot("60s")+ggtitle("60s"),
                            get_plot("70s")+ggtitle("70s"),
                            get_plot("80s")+ggtitle("80s"),
                            get_plot("90_100s")+ggtitle("90s+"),ncol = 2,align = "v")

ggsave(plot = g_averteddeath,filename = "./output/Figure 2.png",width = 15,height = 9)


