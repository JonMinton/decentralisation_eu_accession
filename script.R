# 27/7/2016

rm(list = ls())

# Script to compare distribution of EU  and EU-Accession populations 
# from 2001 to 2011 at English & Welsh LSOAs

require(pacman)


pacman::p_load(
  stringr, car, tidyr, 
  readr, readxl,
  dplyr, purrr,
  ggplot2, tmap,
  RColorBrewer
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

# Remove the numbers from the end of MSOA01NM and we have place names. We can use those 
# 
# lookup %>% 
#   mutate(
#     place_name = str_replace(MSOA01NM, "[0-9]{1,3}$", ""),
#     place_name = str_trim(place_name)
#          ) %>% 
#  group_by(place_name) %>% tally %>% arrange(desc(n)) %>% View


# TTWA lookup 

ttwa_lookup <- read_csv("data/LSOA01_TTWA01_UK_LU.csv")

ttwa_lookup %>% 
  select(LSOA01CD, TTWA01CD, TTWA01NM) %>% 
  inner_join(lookup) %>% 
  select(lsoa = LSOA01CD, msoa = MSOA01CD, ttwa = TTWA01CD, ttwa_name = TTWA01NM) -> ttwa_msoa_lookup

dta %>% 
  left_join(ttwa_msoa_lookup) %>% 
  group_by(ttwa, ttwa_name, msoa, year, geography_short) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(msoa, year) %>% 
  mutate(proportion = count / sum(count)) -> dta_msoa_ttwa

dta %>% 
  left_join(ttwa_msoa_lookup) %>% 
  group_by(ttwa, ttwa_name, lsoa, year, geography_short) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  group_by(lsoa, year) %>% 
  mutate(proportion = count / sum(count)) -> dta_lsoa_ttwa


# shapefile, MSOA

shp_msoa <- read_shape(file = "shapefiles/Middle_layer_super_output_areas_(E+W)_2001_Boundaries_(Full_Extent)_V2/MSOA_2001_EW_BFE_V2.shp")


# Shapefile, LSOA 

shp_lsoa <- read_shape(file = "shapefiles/England_low_soa_2001/england_low_soa_2001.shp")

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

# Do just for Sheffield, MSOA

dta_msoa_ttwa %>% 
  filter(year == 2011) %>% 
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion")

# Do just for Sheffield, LSOA

dta_lsoa_ttwa %>% 
  filter(year == 2001) %>%
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_lsoa, data = ., key.shp = "zonecode", key.data = "lsoa", ignore.na = T) %>% 
  .[!is.na(.$ttwa),] %>% 
  tm_shape(.) + 
  tm_fill("proportion")

dta_lsoa_ttwa %>% 
  filter(year == 2011) %>%
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_lsoa, data = ., key.shp = "zonecode", key.data = "lsoa", ignore.na = T) %>% 
  .[!is.na(.$ttwa),] %>% 
  tm_shape(.) + 
  tm_fill("proportion", breaks = seq(0, 1, by = 0.1))

dta_msoa_ttwa %>% 
  filter(year == 2001) %>% 
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "York") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion")


# For individual cities LSOAs look OK 

dta_msoa_ttwa %>% 
  filter(year == 2011) %>% 
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion")



# Plot MSOAs by proportion UK born in 2001
dta_msoa %>% 
  filter(geography_short == "UK") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    c("p_2001", "p_2011"), 
    palette = brewer.pal(10, "Spectral"),
    breaks = seq(0, 1.0, by = 0.1), showNA = F,
    border.col = NULL
  )


# % point change in proportion from 2001 to 2011
dta_msoa %>% 
  filter(geography_short == "UK") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = p_2011 - p_2001) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = seq(-0.25, 0.25, by = 0.05), showNA = F,
    border.col = NULL
  )

# % change in proportion from 2001 to 2011
dta_msoa %>% 
  filter(geography_short == "UK") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = p_2011 / p_2001) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = c(0.4, 0.6, 0.9, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20, 1.40, 1.8), showNA = F,
    border.col = NULL
  ) 

# Tried to start up leaflet. Did not work.... 

# The outlier - 
# This could be Lakenheath, near RAF Lakenheath
# RAF Lakenheath is a USAF base. 
# There could have been a large reduction in the size of US presence 
# between 2001 and 2011, leading to an increase in the UK born population here. 

#To test this, let's look at the Other proportion

# % point change in proportion from 2001 to 2011
dta_msoa %>% 
  filter(geography_short == "Other") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = p_2011 - p_2001) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = seq(-0.25, 0.25, by = 0.05), showNA = F,
    border.col = NULL
  )
# This clearly indicates the same region as having the largest fall in 'Other' population



# Now, back to looking at change in EU proportion 

dta_msoa %>% 
  filter(geography_short == "Europe - Not Old EU") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = p_2011 - p_2001) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = seq(-0.25, 0.25, by = 0.05), showNA = F,
    border.col = NULL
  )    

# % rather than % point change
dta_msoa %>% 
  filter(geography_short == "Europe - Not Old EU") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = (p_2011 + 0.001)  / (p_2001 + 0.001)) %>% # Continuity correction
  filter(!is.na(change)) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200), showNA = F,
    border.col = NULL
  )    

dta_msoa %>% 
  filter(geography_short == "Europe - Old EU or Western Europe") %>% 
  ungroup() %>% 
  select(msoa, year, proportion) %>% 
  spread(year, proportion) %>% 
  rename(p_2001 = `2001`, p_2011 = `2011`) %>% 
  mutate(change = (p_2011 + 0.001)  / (p_2001 + 0.001)) %>% # Continuity correction
  filter(!is.na(change)) %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>% 
  tm_shape(.) + 
  tm_polygons(
    "change", 
    palette = brewer.pal(10, "Spectral"),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1, 2, 5, 10), showNA = F,
    border.col = NULL
  )    


# This definitely points to the two population groups being 
# very different in their migration trends between 2001 and 2011


# Sheffield TTWA 

lookup %>% 
  inner_join(ttwa_lookup)
  
  
  
  dta_msoa
dta
# 