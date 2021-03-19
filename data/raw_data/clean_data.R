#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)

# import data
survey_8617 <- read_csv("data/raw_data/myvatn_char_1986_2017.csv")
survey_1820 <- read_csv("data/raw_data/myvatn_char_2018_2020.csv")




#==========
#========== Survey data
#==========

# clean data
# replace year of 0 with 2002 (probably a mistake in the data)
# only keep fall surveys
clean_8617 <- survey_8617 %>%
  rename(number = NR, year = AR, month = TIMI, area_1 = ST, area_2 = STOD, 
         length = LENGD, mass = TYNGD, 
         sex = KYN, maturity = KTR, age = ALD, est_age = ALDUR) %>%
  select(number, year, month, area_1, area_2, length, mass, sex, maturity, age, est_age) %>%
  mutate(sex = ifelse(sex == 0, NA, sex),
         age = ifelse(age == 0, NA, age),
         est_age = ifelse(est_age == 0, NA, est_age),
         year = ifelse(year == 0, 2002, year)) 

clean_1820 <- survey_1820 %>%
  mutate(year = year(DAGS),
         month = month(DAGS)) %>%
  rename(number = NR, area_1 = ST, 
         length = LENGD, mass = TYNGD, 
         sex = KYN, maturity = KTR, age = ALD,) %>%
  select(number, year, month, area_1,length, mass, sex, maturity, age) %>%
  mutate(sex = ifelse(sex == 0, NA, sex),
         age = ifelse(age == 0, NA, age),
         year = ifelse(year == 0, 2002, year)) 

clean_8620 <- bind_rows(clean_8617, clean_1820) 

# export
# write_csv(clean_8620, "data/myvatn_char_clean_1986_2020.csv")
