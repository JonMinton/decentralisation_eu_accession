# 27/7/2016

rm(list = ls())

# Script to compare distribution of EU  and EU-Accession populations 
# from 2001 to 2011 at English & Welsh LSOAs

require(pacman)


pacman::p_load(
  stringr, car, tidyr, 
  readr, readxl,
  dplyr, purrr,
  ggplot2, tmap
)



# Tidy 2001 data  ---------------------------------------------------------

dta_01 <- read_csv("data/2001_census/cob_2001_sex.csv")

value_codes <- read_excel("data/2001_census/value_code_labels.xlsx")
eu_groupings <- read_excel("data/2001_census/eu_2001_categories.xlsx", na = "NA")
geo_lookup <- read_csv("data/2001_geographical_lookup/OA01_LSOA01_MSOA01_EW_LU.csv")



value_codes %>% mutate(
  code = value + 150000,
  code = paste0("cs0", as.character(code))
) -> value_codes

dta_01 %>% 
  gather("code", "count", -`Zone Code`) %>% 
  left_join(value_codes) %>% 
  rename(OA01CD = `Zone Code`) %>% 
  left_join(geo_lookup) %>% 
  select(lsoa = LSOA01CD, geography, sex, count) %>% 
  group_by(lsoa, geography, sex) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>%   # .$geography %>% unique %>% write.csv(file = "clipboard", .)
  left_join(eu_groupings) %>% 
  filter(sex == "total") %>% 
  filter(!is.na(geography_short)) %>% 
  group_by(lsoa, geography_short) %>% 
  summarise(count = sum(count)) -> simplified_2001






# Tidy 2011 data  ---------------------------------------------------------



