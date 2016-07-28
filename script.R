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


dta <- read_csv("data/derived/simplified_both.csv")

# 
# shp_eng <- read_shape(file = "shapefiles/England_low_soa_2001/england_low_soa_2001.shp")
# 
# dta  %>% filter(year == 2001, geography_short == "UK")  -> tmp
# 
# tmp %>% append_data(data = ., shp = shp_eng, key.shp = "zonecode", key.data = "lsoa") -> 
#   shp_eng_2001
# 
# shp_eng_2001 %>% 
#   tm_shape(.) + tm_fill(col = "proportion")
# 


# LSOAs seem too highly detailed. Looking to move to MSOAs instead

# I already have LSOAs to MSOA lookup for 2001

lookup <- read_csv("data/2001_geographical_lookup/OA01_LSOA01_MSOA01_EW_LU.csv")
         
lookup %>% 
  select(lsoa = LSOA01CD, msoa = MSOA01CD) %>% 
  distinct() -> simple_lookup


dta %>% 
  left_join(simple_lookup) %>% 
  select(msoa, year, geography_short, count) %>% 
  group_by(msoa, year, geography_short) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(msoa, year) %>% 
  mutate(proportion = count / sum(count)) -> dta_msoa

 
# 

shp_msoa <- read_shape(file = "shapefiles/Middle_layer_super_output_areas_(E+W)_2001_Boundaries_(Full_Extent)_V2/MSOA_2001_EW_BFE_V2.shp")

# # Create cartogram
# dta_msoa  %>% 
#   filter(year == 2001)  %>% 
#   select(msoa, count)  %>% 
#   group_by(msoa)  %>% 
#   summarise(count = sum(count))  %>% 
#   ungroup() -> pop_2001
# 
# append_data(shp = shp_msoa, data = pop_2001, key.shp = "MSOA01CD", key.data = "msoa") -> shp_msoa_count
# 
# shp_msoa_count[!is.na(shp_msoa_count$msoa),] -> shp_msoa_count
# 
# write_shape(shp = shp_msoa_count, file = "shapefiles/msoa_2001_pops/shp_2001_pop.shp")


# Append proportion UK born in different places 

dta_msoa %>% 
  filter(year == 2001, geography_short == "UK") %>% 
  ungroup() %>% 
  select(msoa, proportion) -> tmp 

append_data(shp = shp_msoa, data = tmp, key.shp = "MSOA01CD", key.data = "msoa") %>% 
  tm_shape(.) + 
  tm_polygons("proportion", border.col = NULL)




