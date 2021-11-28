setwd("~/covid19_cfr_ve_tokyo")
source("./src/utils.R")

#############################################################################################################
# Epi curves and the proportion of fully vaccinated                                                         #
#                                                                                                           #
## The code was written by Yura Ko K                                                                        #
#############################################################################################################

get_plot <- function(age_,age_name){
  
  data_case <- df_case_vac_prop %>% filter(age == age_)
  sec_case <- max(data_case$n)/100
  
  data_death <- df_death_vac_prop %>% filter(age == age_)
  sec_death <- max(data_death$n)/100
  
  g_case <- 
    ggplot(NULL) +
    ggtitle(age_name)+
    xlab("Date of confirmation")+
    ylab("Number of \ncases/deaths")+
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d",
      expand = c(0,0)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~ ./sec_case , name = "")) +
    geom_bar(data=data_case,
             aes(x= date_confirm, 
                 y=n
             ),
             stat="identity", 
             position="stack",
             fill="#fdae61",
             width=1,
             alpha=0.3
    )+
    geom_line(data=data_case,
              aes(x= date_confirm, 
                  y=prop_2vac*sec_case
              ),
              colour="#fdae61",
              size = 1
    )+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=18, family="sans",color="black"),
          axis.text.x =element_text(size=10),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  g_death <- 
    ggplot(NULL) +
    ggtitle(age_name)+
    xlab("Date of death")+
    ylab("")+
    scale_x_date(
      date_breaks = "1 months",
      date_labels = "%b-%d",
      expand = c(0,0)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~ ./sec_death , name = "Proportion of \nfully vaccinated (%)")) +
    geom_bar(data=data_death,
             aes(x= date_death, 
                 y=n
             ),
             stat="identity", 
             position="stack",
             fill="#d7191c",
             width=1,
             alpha=0.3
    )+
    geom_line(data=data_death,
              aes(x= date_death, 
                  y=prop_2vac*sec_death
              ),
              colour="#d7191c",
              size = 1
    )+
    theme_bw() +
    theme(axis.text.x = element_text()) +
    theme(text=element_text(size=18, family="sans",color="black"),
          axis.text.x =element_text(size=10),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(list(g_case,g_death))
}

plot_30_50s <- get_plot("30_50s","30-50s")
plot_60s <- get_plot("60s","60s")
plot_70s <- get_plot("70s","70s")
plot_80s <- get_plot("80s","80s")
plot_90_100s <- get_plot("90_100s","90s+")

plot_all <- plot_grid(plot_30_50s[[1]],plot_30_50s[[2]],
                       plot_60s[[1]],plot_60s[[2]],
                       plot_70s[[1]],plot_70s[[2]],
                       plot_80s[[1]],plot_80s[[2]],
                       plot_90_100s[[1]],plot_90_100s[[2]],ncol = 2,align = "v")

ggsave(plot=plot_all,filename = "./output/Figure 1.png",width = 12,height = 15)

