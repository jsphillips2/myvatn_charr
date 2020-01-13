#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# import data
data = read_csv("data/myvatn_char_clean.csv")

# set theme
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))





#==========
#========== Prepare Data
#==========

data = data %>%
  #Fill gaps in "Est.age"
  mutate(est_age = ifelse(is.na(age)==F, age, est_age)
  ) %>%
  #Remove June data & samples without ages
  filter(month != 6, is.na(est_age)==F) %>%
  #Define Stage
  mutate(stage = factor(ifelse(est_age == 1, "first", 
                          ifelse(est_age==2, "second", 
                                 ifelse(est_age==3, "third", "adult"))))
         )


#count number of fish per stage per year
dd = data %>%
  group_by(year, stage) %>%
  #Count number of fish
  summarize(
    count = length(length) 
  ) %>%
  #Merge with stage_year to keep track of 0's
  full_join(data %>% expand(stage, year)) %>% 
  #Repalce NA's with 1's
  #This only affects data for 1st year fish, which aren't used in the model
  mutate(
    count = ifelse(is.na(count)==T, 1L, count) 
  ) %>%
  ungroup()

#plot
dd %>% 
ggplot(aes(year, count, color=stage))+ 
  geom_line(size=0.8, alpha=0.8)+
  scale_color_manual(values=c("black","gray50","dodgerblue2","firebrick2"))+
  scale_y_continuous("Char Count",trans="log", limits=c(1,600), breaks=c(3,12,48,192))

dd_wide = dd %>%
  # log transformed
  mutate(count = count) %>%
  # convert to wide format
  spread(stage, count) 





#==========
#========== Package and Export
#==========

# define variables in environment
N_obs = dd_wide[,c("first","second","third","adult")] %>% t()
n_obs = log(N_obs)
nYears = ncol(n_obs)
nStages = nrow(n_obs)

# export
# stan_rdump(c("nYears", "nStages", "N_obs", "n_obs"), file = "model/data_list.R")
# write_csv(dd, "data/myvatn_char_counts.csv")












