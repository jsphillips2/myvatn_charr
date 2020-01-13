#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import data
survey = read_csv("data/raw_data/myvatn_char.csv")
catch = read_csv("data/raw_data/myvatn_catch.csv")





#==========
#========== Survey data
#==========

# clean data
# replace year of 0 with 2002 (probably a mistake in the data)
# only keep fall surveys
survey_clean = survey %>%
  rename(number = NR, year = AR, month = TIMI, area_1 = ST, area_2 = STOD, length = LENGD, mass = TYNGD, 
         sex = KYN, maturity = KTR, age = ALD, est_age = ALDUR) %>%
  select(number, year, month, area_1, area_2, length, mass, sex, maturity, age, est_age) %>%
  mutate(sex = ifelse(sex == 0, NA, sex),
         age = ifelse(age == 0, NA, age),
         est_age = ifelse(est_age == 0, NA, est_age),
         year = ifelse(year == 0, 2002, year)) 

# export
# write_csv(survey_clean, "data/myvatn_char_clean.csv")





#==========
#========== Survey data
#==========

# rename variables
# replace numeric indeces with strings
catch_clean = catch %>%
  rename(farm = NUMER, method = VEID, season = TIMI, year = ÁR, month = MAN, 
         day = DAG, week = VIKA, location = STADUR, basin = FLOI, net_nights = LAGNIR, 
         char = BLEIKJA, trout = URRIÐI, filter = FILTER__) %>%
  mutate(method = ifelse(method == 1, "net", "hook"),
         season = ifelse(season == 1, "winter", "summer")) %>%
  select(-VAR00001,-X15)

# export
# write_csv(catch_clean, "data/myvatn_catch_clean.csv")









