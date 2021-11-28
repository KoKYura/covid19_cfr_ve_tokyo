setwd("~/covid19_cfr_ve_tokyo")
source("./src/utils.R")

#############################################################################################################
# Comparison of CFR estimated from the indicators of the level of health care burden and asymptomatic rate, #
# with CFR jointly estimated with VE by Bayesian estimation from the daily number of positive cases and     #
# deaths by vaccination status.                                                                             #
#                                                                                                           #
## The code was written by Yura Ko K                                                                        #
#############################################################################################################

# 1. Import and processing data ####

## Jointly estimated unvaccinated CFR ####

df_joint_estimated_cfr <- read_excel("./output/res_jointly-estimated-cfr.xlsx") %>% mutate(date = as.Date(date)) %>% 
  rename(je_cfr_m = "result_median", je_cfr_l = "result_2.5", je_cfr_u = "result_97.5")

## Predicted CFR by the regression of the health care burden indicator selected by WAIC####

df_regression_cfr_70s <- read_xlsx("./output/result_70s_severe_lag3.xlsx") %>% mutate(date = as.Date(date)) %>% mutate(age = "70s")
df_regression_cfr_80s <- read_xlsx("./output/result_80s_severe_lag1.xlsx") %>% mutate(date = as.Date(date)) %>% mutate(age = "80s")
df_regression_cfr_90s_ <- read_xlsx("./output/result_90_100s_severe_lag1.xlsx") %>% mutate(date = as.Date(date)) %>% mutate(age = "90_100s")

df_regression_cfr <- rbind(df_regression_cfr_70s,df_regression_cfr_80s,df_regression_cfr_90s_) %>% 
  rename(reg_cfr_m = "median", reg_cfr_l = "lower", reg_cfr_u = "upper")

## Combine both the results of CFR####

df_cfr <- left_join(df_joint_estimated_cfr,df_regression_cfr, by = c("date","age")) %>% 
  filter(age == "70s"|age == "80s"|age == "90_100s")

# 2. Plot Figure 3  ####

get_plot <- function(age_){
  data_ <- df_cfr %>% filter(age == age_)
  g_comp <- 
    data_ %>% 
    ggplot() +
    geom_ribbon(aes(x = date, ymin = je_cfr_l*100, ymax = je_cfr_u*100), fill = "#2c7fb8", alpha = 0.3) +
    geom_ribbon(aes(x = date, ymin = reg_cfr_l*100, ymax = reg_cfr_u*100), fill = "white", alpha = 0.3) +
    geom_ribbon(aes(x = date, ymin = reg_cfr_l*100, ymax = reg_cfr_u*100), fill = "#dd1c77", alpha = 0.3) +
    geom_line(aes(x = date, y = je_cfr_m*100), color="#2c7fb8", size = 1, alpha=1) +
    geom_line(aes(x = date, y = reg_cfr_m*100), color="#dd1c77", size = 1, alpha=1) +
    labs(x="",y="Time delay adjusted CFR [%]\n") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12)) +
    theme(text=element_text(size=12, family="sans",color="black"),
          axis.text=element_text(angle = 0),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d"
    ) 

  return(g_comp)
}

plot_all <- plot_grid(get_plot("70s")+ggtitle("70s")+coord_cartesian(ylim =c(0,15)),
                      get_plot("80s")+ggtitle("80s")+coord_cartesian(ylim =c(5,30)),
                      get_plot("90_100s")+ggtitle("90s+")+coord_cartesian(ylim =c(5,40)),
                      ncol = 1,align = "v")
ggsave(plot = plot_all,filename = "./output/Figure 3.png",width = 12,height = 9)

