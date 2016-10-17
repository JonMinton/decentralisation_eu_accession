# New attempt to combine and tidy data from 2001 and 2011 census 

rm(list = ls())

# Script to compare distribution of EU  and EU-Accession populations 
# from 2001 to 2011 at English & Welsh LSOAs

require(pacman)


pacman::p_load(
  stringr, forcats, tidyr, 
  readr, readxl,
  dplyr, purrr,
  ggplot2, tmap
)





# 2001 data - initial exploration -----------------------------------------


dta_01 <- read_csv("data/2001_census/20161017123252419_CTRBIR_UNIT/Data_CTRBIR_UNIT.csv")

dta_01_mta <- read_csv("data/2001_census/20161017123252419_CTRBIR_UNIT/Meta_CTRBIR_UNIT.csv")

dta_01

# Need to find those rows of dta_01_mta where topic == Country of birth

dta_01_mta %>% 
  filter(TOPIC == "Country of birth") %>% 
  select(CDU_FIELD_NAME, CATEGORY) -> key_value


# Now to slim down dta_01
dta_01 %>% 
  select(GEO_CODE, F15589:F96215) %>% 
  slice(-1) %>% 
  gather(code, count, -GEO_CODE) %>% 
  mutate(count = as.numeric(count)) %>% 
  inner_join(key_value, by = c("code" = "CDU_FIELD_NAME")) %>% 
  select(msoa = GEO_CODE, code, category = CATEGORY, count) -> census_2001_tidy

census_2001_tidy %>% group_by(category) %>% summarise(count = sum(count))




# Now to attempt something similar for 2011 

dta_11 <- read_csv("data/2011_census/20161017151535831_COBDET_UNIT/Data_COBDET_UNIT.csv")

dta_11_mta <- read_csv("data/2011_census/20161017151535831_COBDET_UNIT/Meta_COBDET_UNIT.csv")

dta_11_mta %>% 
  filter(TOPIC == "Country of birth (condensed for England and Wales) [E][W]") %>% 
  select(CDU_FIELD_NAME, CATEGORY) -> key_value

dta_11  %>% 
  select(GEO_CODE, F1219:F1290) %>% 
  slice(-1) %>% 
  gather(code, count, -GEO_CODE) %>% 
  mutate(count = as.numeric(count)) %>% 
  inner_join(key_value, by = c("code" = "CDU_FIELD_NAME")) %>% 
  select(msoa = GEO_CODE, code, category = CATEGORY, count) -> census_2011_tidy

census_2011_tidy %>% group_by(category) %>% summarise(count = sum(count, na.rm =T))


# Simplify 2001 data

lookup_01 <- read_excel("data/eu_countries_consistent_categories.xlsx", "c_2001")
lookup_11 <- read_excel("data/eu_countries_consistent_categories.xlsx", "c_2011")

census_2001_tidy %>% 
  inner_join(lookup_01, by = c("category" = "geography")) %>% 
  filter(geography_short != "NA") %>% 
  group_by(msoa, geography_short) %>% 
  summarise(count = sum(count)) -> census_2001_tidy

census_2011_tidy %>% 
  inner_join(lookup_11, by = c("category" = "country")) %>%
  .[,1:5] %>% 
  select(msoa, category, count, geography_short = country_short) %>% 
  filter(geography_short != "NA") %>% 
  group_by(msoa, geography_short) %>% 
  summarise(count = sum(count)) -> census_2011_tidy


census_2001_tidy %>% 
  rename(cob = geography_short) %>% 
  mutate(census = 2001) %>% 
  select(census, msoa, cob, count) -> census_2001_tidy

census_2011_tidy %>% 
  rename(cob = geography_short) %>% 
  mutate(census = 2011) %>% 
  ungroup() %>%
  select(census, msoa, cob, count) -> census_2011_tidy

census_both <- bind_rows(census_2001_tidy, census_2011_tidy)


write_csv(x = census_both, "data/derived/cob_both_censuses_simplified.csv")



