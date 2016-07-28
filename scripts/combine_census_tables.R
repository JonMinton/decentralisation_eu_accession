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
eu_groupings <- read_excel("data/eu_countries_consistent_categories.xlsx",sheet = "c_2001", na = "NA")
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
  summarise(count = sum(count)) %>% 
  mutate(year = 2001) %>% 
  select(lsoa, year, geography_short, count) -> simplified_2001


write_csv(simplified_2001, "data/derived/simplified_2001.csv")

rm(list = ls())
gc()

# Tidy 2011 data  ---------------------------------------------------------

dta_2011 <- read_csv("data/2011_census/201672713242535_AGE_COB_ECOACT_UNIT/Data_AGE_COB_ECOACT_UNIT.csv")
#meta_2011 <- read_csv("data/2011_census/201672713242535_AGE_COB_ECOACT_UNIT/Meta_AGE_COB_ECOACT_UNIT.csv")

eu_groupings <- read_excel("data/eu_countries_consistent_categories.xlsx",sheet = "c_2011", na = "NA")

ln_1 <- dta_2011 %>%  
  .[1,]  %>% 
  map_chr( ~.[[1]]) # convert to character vector


nms <- names(dta_2011)
nms2 <- ifelse(is.na(ln_1), nms, ln_1)
names(dta_2011) <- nms2

rm(nms, nms2, ln_1)
dta_2011 %>% 
  slice(-1) %>% 
  gather(key = "category", value = "count", -c(1:5)) %>% 
  filter(!is.na(count)) %>% 
  select(lsoa = GEO_CODE, category, count) %>% 
  separate(category, into = c("age", "country", "economic_activity", "unit"), sep = "-") %>% 
  filter(str_detect(age, "16 and")) %>%  # For 'Age : Age 16 and over'
  select(lsoa, country, count) %>% 
  mutate(country = str_replace(country, pattern = "Country of birth : ", replacement = "")) %>% 
  mutate(country = str_trim(country)) %>% 
  left_join(eu_groupings) %>% 
  select(lsoa, country_short, count) %>% 
  mutate(count = as.numeric(count)) %>%  
  group_by(lsoa, country_short) %>%
  summarise(count = sum(count)) %>% 
  filter(!is.na(country_short)) %>% 
  rename(geography_short = country_short) %>% 
  mutate(year = 2011) %>% 
  select(lsoa, year, geography_short, count) -> simplified_2011

write_csv(simplified_2011, path = "data/derived/simplified_2011.csv")


rm(list = ls())
gc()



# Join the two tables -----------------------------------------------------

d_01 <- read_csv("data/derived/simplified_2001.csv")
d_11 <- read_csv("data/derived/simplified_2011.csv")

d_both <- bind_rows(d_01, d_11)

d_both %>% 
  group_by(lsoa, year) %>% 
  mutate(proportion = count / sum(count)) %>% 
  ungroup() -> d_both

write_csv(d_both, "data/derived/simplified_both.csv")





         







         
