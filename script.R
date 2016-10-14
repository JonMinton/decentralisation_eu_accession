# 27/7/2016

rm(list = ls())

# Script to compare distribution of EU  and EU-Accession populations 
# from 2001 to 2011 at English & Welsh LSOAs

require(pacman)


pacman::p_load(
  stringr, car, tidyr, 
  readr, readxl,
  dplyr, purrr,
  shapefiles, sp,
  spdep, rgeos, maptools,
  tmap, classInt,
  ggplot2, 
  RColorBrewer,
  CARBayes, MCMCpack,
  truncdist
)

#### Source the functions for running the model
source("scripts/RCI.R")
source("scripts/binomial.MCARleroux.R")
Rcpp::sourceCpp("scripts/aqmen.cpp")


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
  dplyr::select(LSOA01CD, TTWA01CD, TTWA01NM) %>% 
  inner_join(lookup) %>% 
  dplyr::select(lsoa = LSOA01CD, msoa = MSOA01CD, ttwa = TTWA01CD, ttwa_name = TTWA01NM) -> ttwa_msoa_lookup

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



# Automate process 
# For each TTWA 

# Proportion UK born 
# Proportion Old EU born
# Proportion New EU Born
#   - in 2001 and 2011 
create_ttwa_shp <- function(dta, shp) {
  dta %>% 
    append_data(shp = shp, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na = T) %>% 
    .[!is.na(.$ttwa),] -> out
  out
}

create_ttwa_tmap <- function(shp, title){
   tm_shape(shp) + 
    tm_fill("proportion", palette = "Paired", style = "quantile", n = 10) + 
    tm_legend(legend.outside = T, legend.outside.position = "right") -> tm
  
  save_tmap(tm, filename = paste0("tmaps/proportions/",title, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  NULL
}

dta_msoa_ttwa %>% 
  filter(!is.na(year), !is.na(ttwa_name)) %>% 
  group_by(ttwa_name, year, geography_short) %>% 
  nest() %>% 
  mutate(ttwa_shp = map(data, create_ttwa_shp, shp = shp_msoa)) %>% 
  mutate(title_name = paste0(ttwa_name, "_" , year, "_" , geography_short)) %>%  
  mutate(tmp = walk2(ttwa_shp, title_name, create_ttwa_tmap))


create_change_ttwa_tmap <- function(shp, title){
  tm_shape(shp) + 
    tm_fill("change", palette = "Paired", style = "quantile", n = 10) + 
    tm_legend(legend.outside = T, legend.outside.position = "right") -> tm
  
  save_tmap(tm, filename = paste0("tmaps/change_over_censuses/",title, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  NULL
}

debug(create_change_ttwa_tmap)
dta_msoa_ttwa %>% 
  filter(!is.na(year), !is.na(ttwa_name)) %>% 
  dplyr::select(ttwa, ttwa_name, msoa, year, geography_short, proportion) %>% 
  spread(year, proportion) %>%
  mutate(change = `2011` - `2001`) %>% 
  filter(!is.na(change)) %>% 
  dplyr::select(-`2001`, -`2011`) %>% 
  ungroup() %>% 
  group_by(ttwa_name, geography_short) %>% 
  nest() %>% 
  mutate(title_name = paste0(ttwa_name, "_", geography_short)) %>% 
  mutate(ttwa_shp = map(data, create_ttwa_shp, shp = shp_msoa)) -> tmp

tmp %>% 
  mutate(tmp = walk2(ttwa_shp, title_name, create_change_ttwa_tmap))

 



# Append proportion UK born in different places 

# Do just for Sheffield, MSOA

dta_msoa_ttwa %>% 
  filter(year == 2011) %>% 
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Proportion UK Born, 2011")

dta_msoa_ttwa %>% 
  filter(year == 2001) %>% 
  filter(geography_short == "UK") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Proportion UK Born, 2001")

dta_msoa_ttwa %>% 
  filter(year == 2001) %>% 
  filter(geography_short == "Europe - Old EU or Western Europe") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Old EU, 2001", palette = "Paired", breaks = seq(0, 0.25, by = 0.0125)) + 
  tm_legend(legend.outside = T, legend.outside.position = "right")

dta_msoa_ttwa %>% 
  filter(year == 2011) %>% 
  filter(geography_short == "Europe - Old EU or Western Europe") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Old EU, 2011", palette = "Paired", breaks = seq(0, 0.25, by = 0.0125))


dta_msoa_ttwa %>% 
  filter(year == 2001) %>% 
  filter(geography_short == "Europe - Not Old EU") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Accession States? 2001", palette = "Paired", breaks = seq(0, 0.0125, by = 0.00125)) + 
  tm_legend()

dta_msoa_ttwa %>% 
  filter(year == 2011) %>% 
  filter(geography_short == "Europe - Not Old EU") %>% 
  filter(ttwa_name == "Sheffield and Rotherham") %>% 
  append_data(shp = shp_msoa, data = ., key.shp = "MSOA01CD", key.data = "msoa", ignore.na =T) %>%  
  .[!is.na(.$ttwa),]  %>% 
  tm_shape(.) + 
  tm_fill("proportion", title = "Accession States, 2011", palette = "Paired", breaks = seq(0, 0.00625, by = 0.0125))


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




# Attempt to replicate Gavin's results with his data and code -------------


#### Check the boundary changes for different TTWAs. if not much not unchanged LSOAs or DATA zones would be used
### Import the English COB data in 2001
cob2001 <- read_csv("data/from_gavin/England_LSOA_2001.csv")

### Import the English COB data in 2001
cob2011 <- read_csv("data/from_gavin/England_LSOA_2011.csv")

#### For the three cities (London, Manchester and Sheffield), check the LSOA boundary changes during the two census years
# link LSOAs to the TTWA defined in 2007
lsoa_to_ttwa <- read_csv("data/from_gavin/LSOA01_TTWA07_UK_LU.csv")

# changes of LSOAs boundaries
changes_boundary <- read_csv("data/from_gavin/LSOA01_LSOA11_EW_LUv2.csv") %>%
  select(LSOA01CD,LSOA11CD,CHGIND)
# 
# 
# 
# match.ind <- match(cob2001$code,lsoa_to_ttwa$LSOA01CD)
# sum(is.na(match.ind))
# tmp <- cob2001
# tmp$ttwa <- lsoa_to_ttwa$TTWA07NM[match.ind]
# data.frame(table(tmp$ttwa)) -> tmp1

# cob2001 %>% 
#   inner_join(lsoa_to_ttwa, by = c("code"= "LSOA01CD")) %>%
#   rename(ttwa = TTWA07NM) %>% 
#   group_by(ttwa) %>% tally() %>% data.frame -> tmp2

# Check all elements in tmp1 and tmp2 are identical:
# sapply(tmp1 == tmp2, function(x) {x} )  %>% all(. == T)

cob2001 %>% 
  inner_join(
    lsoa_to_ttwa, by = c("code" = "LSOA01CD")
  ) %>% 
  mutate(
    eu15_01    = Republic_of_Ireland + Other_EU_countries,
    uk_01      = England + Scotland + Wales + Northern_Ireland,
    total_01   = All_people,
    lsoa_01_cd = code,
    eu12_01    = NA
  ) %>%
  select(
    lsoa_01_cd, lsoa_01_nm = LSOA01NM,
    ttwa_07_cd = TTWA07CD, ttwa_07_nm = TTWA07NM, 
    total_01, eu15_01, eu12_01, uk_01
  ) -> cob_01_tdy

cob2011 %>% 
  mutate(
    eu15_11    = Ireland + EU_2001,
    eu12_11    = EU_2001_TO_2011,
    uk_11      = UK,
    total_11   = all,
    lsoa_11_cd = code,
    ttwa_07_cd = NA,
    ttwa_07_nm = NA
  ) %>%
  select(
    lsoa_11_cd, lsoa_11_nm = Name,
    ttwa_07_cd, ttwa_07_nm,
    total_11, eu15_11, eu12_11, uk_11
    ) -> cob_11_tdy



# link to the data of boundary changes
changes_boundary %>% 
  rename(lsoa_01_cd = LSOA01CD, lsoa_11_cd = LSOA11CD) %>% 
  left_join(cob_01_tdy, by = "lsoa_01_cd") %>% 
  left_join(cob_11_tdy, by = "lsoa_11_cd") -> boundary_changes_1


### we use LSOA boundary in 2001 as baseline geographical units
# case 1; split--one lsoa was split into two or more lsoa in 2011. in this case, aggregate the COB data in 2011 based on 2001 lsoas
boundary_changes_1 %>% filter(CHGIND == "S")
#grepl("*.2011",names(changes.boundary.2))

split_data  <- boundary_changes_1 %>% 
  filter(CHGIND == "S") %>%
  select(lsoa_01_cd, total_11, uk_11, eu15_11, eu12_11) %>%
  group_by(lsoa_01_cd) %>%
  summarise_each(funs(sum)) 


# # link back to the
# s.ind <- changes.boundary.1$CHGIND =="S"
# match.ind <- match(changes.boundary.1$LSOA01CD[s.ind],data.split$LSOA01CD)
# changes.boundary.1[s.ind,c("all_2011","UK_2011","EU15_2011","EU12_2011")] <- data.split[match.ind,2:5]

# case 2: two or more LSOAs were merged into a single LSOA in 2011. in this case we need aggregate the 2001 LSOAs accordingly.

merged_data <- boundary_changes_1 %>% 
  filter(CHGIND == "M") %>%
  select(lsoa_11_cd, total_01, uk_01,eu15_01, eu12_01) %>%
  group_by(lsoa_11_cd) %>%
  summarise_each(funs(sum))

unchanged_data <-  boundary_changes_1 %>% 
  filter(CHGIND == "U") %>% 
  select(
    lsoa_01_cd, lsoa_01_nm, CHGIND, lsoa_11_cd, lsoa_11_nm, 
    ttwa_07_cd = ttwa_07_cd.x, ttwa_07_nm = ttwa_07_nm.x,
    total_01, eu15_01, eu12_01, uk_01,
    total_11, eu15_11, eu12_11, uk_11
    )

boundary_changes_1 %>% 
  filter(CHGIND == "S") %>% 
  select(
    lsoa_01_cd, lsoa_01_nm, CHGIND, lsoa_11_cd, lsoa_11_nm, 
    ttwa_07_cd = ttwa_07_cd.x, ttwa_07_nm = ttwa_07_nm.x,
    total_01, eu15_01, eu12_01, uk_01
  ) %>% 
  left_join(split_data) %>% 
  select(
    lsoa_01_cd, lsoa_01_nm, CHGIND, lsoa_11_cd, lsoa_11_nm, 
    ttwa_07_cd,  ttwa_07_nm,
    total_01, eu15_01, eu12_01, uk_01,
    total_11, eu15_11, eu12_11, uk_11
  ) -> all_split_data

boundary_changes_1 %>% 
  filter(CHGIND == "M") %>% 
  select(
    lsoa_01_cd, lsoa_01_nm, CHGIND, lsoa_11_cd, lsoa_11_nm, 
    ttwa_07_cd = ttwa_07_cd.x, ttwa_07_nm = ttwa_07_nm.x,
    total_11, eu15_11, eu12_11, uk_11
  ) %>% 
  left_join(merged_data) %>% 
  select(
    lsoa_01_cd, lsoa_01_nm, CHGIND, lsoa_11_cd, lsoa_11_nm, 
    ttwa_07_cd,  ttwa_07_nm,
    total_01, eu15_01, eu12_01, uk_01,
    total_11, eu15_11, eu12_11, uk_11
    ) -> all_merged_data

unchanged_data %>%
bind_rows(all_merged_data) %>%
bind_rows(all_split_data) -> all_recombined_data

all_recombined_data %>% 
  unique()

# 
# # link back
# m.ind <- changes.boundary.1$CHGIND =="M"
# match.ind <- match(changes.boundary.1$LSOA11CD[m.ind],data.merge$LSOA11CD)
# changes.boundary.1[m.ind,names(data.merge)[-1]] <- data.merge[match.ind,2:4]
# 
# changes.boundary.1  %>%
#   filter(ttwa %in% c("London","Manchester","Sheffield & Rotherham")) %>%
#   group_by(ttwa,CHGIND) %>%
#   summarise(num=length(ttwa))

# select unique ids of LSOAs in 2001
# COB.data <- changes.boundary.1[!duplicated(changes.boundary.1$LSOA01CD),] %>% data.frame
all_recombined_data %>% 
  filter(!duplicated(lsoa_01_cd)) -> tidy_count_data

rm(all_recombined_data, all_merged_data, all_split_data, split_data, merged_data)


# drop complex boundary changes, they only account for a very small proportion of the data
#COB.data <- COB.data[!COB.data$CHGIND == "X",]
# COB.data %>% tbl_df() %>% filter(ttwa %in% c("London","Manchester","Sheffield & Rotherham")) %>%
#   filter(CHGIND == "X") %>%
#   group_by(ttwa) %>%
#   summarise(num=length(ttwa))

# COB.data <- COB.data[!COB.data$CHGIND == "X",]
################################################################################################

######################
###### Sheffield #####
######################
# 
# sheffield.map <- readOGR(dsn=".",layer="Sheffield_ttwa",stringsAsFactors = FALSE)
# sheffield.map$lsoa <- as.character(sheffield.map$LSOA04CD)

# JM : LSOA map 

lsoa_eng <- read_shape(file = "shapefiles/England_low_soa_2001/england_low_soa_2001.shp")

lsoa_eng %>% 
  append_data(
    key.shp = "zonecode", key.data = "lsoa_01_cd", ignore.na = T,
    data = tidy_count_data, 
    shp = . 
  ) %>% 
  split(f = .@data$ttwa_07_nm) -> eng_by_ttwa


# cob data
# sheffield.data <- COB.data[COB.data$ttwa == "Sheffield & Rotherham",]
# head(sheffield.data)
sheffield_data <- eng_by_ttwa[["Sheffield & Rotherham"]]@data %>% tbl_df
sheffield_data 

# sheffield.data[sheffield.data$CHGIND == "M",]
sheffield_data[sheffield_data$CHGIND == "M",]

# length(unique(sheffield.data$LSOA11CD))
# 
# # first link cob to the spatial data, then dissolve the polygons
# match.ind <- match(sheffield.map$lsoa,sheffield.data$LSOA01CD)
match.ind <- match(eng_by_ttwa[["Sheffield & Rotherham"]]$lsoa_01_cd, sheffield_data$lsoa_01_cd)
# # drop the two polygons
eng_by_ttwa[["Sheffield & Rotherham"]] <- eng_by_ttwa[["Sheffield & Rotherham"]][!is.na(match.ind),]
eng_by_ttwa[["Sheffield & Rotherham"]]@data <- data.frame(
  eng_by_ttwa[["Sheffield & Rotherham"]]@data,
  eng_by_ttwa[["Sheffield & Rotherham"]][match.ind[!is.na(match.ind)],]
  )
# head(sheffield.map@data)

#### Now dissolve the Sheffield ttwa polygons based on the ids of LSOAs in 2011
# library(rgeos)
# if (rgeosStatus()) {
  tmp_map <- unionSpatialPolygons(
    eng_by_ttwa[["Sheffield & Rotherham"]], 
    IDs = as.character(eng_by_ttwa[["Sheffield & Rotherham"]]$lsoa_11_cd)
  )
#  temp.map <- unionSpatialPolygons(sheffield.map, IDs = as.character(sheffield.map$LSOA11CD)) 
# }
# if(rgeosStatus()) {
  tmp_df <- as(
    eng_by_ttwa[["Sheffield & Rotherham"]]@data
    , "data.frame"
  )[!duplicated(eng_by_ttwa[["Sheffield & Rotherham"]]$lsoa_11_cd),]
#  temp.df <- as(sheffield.map@data,"data.frame")[!duplicated(sheffield.map$LSOA11CD),]
#  row.names(temp.df) <- temp.df$LSOA11CD
  row.names(tmp_df) <- tmp_df$lsoa_11_cd
#    sheffield.map.final <- SpatialPolygonsDataFrame(temp.map,temp.df)
  eng_by_ttwa[["Sheffield & Rotherham"]] <- SpatialPolygonsDataFrame(
    tmp_map, tmp_df   
  )
# }


# save.image("segregation14_06.RData")
# load("segregation14_06.RData")

##### Mapping the EU migrants distributions
# head(sheffield.map.final@data)

all_t1 <- eng_by_ttwa[["Sheffield & Rotherham"]]$total_01  
all_t2 <- eng_by_ttwa[["Sheffield & Rotherham"]]$total_11

uk_t1 <- eng_by_ttwa[["Sheffield & Rotherham"]]$uk_01
uk_t2 <- eng_by_ttwa[["Sheffield & Rotherham"]]$uk_11

eu15_t1 <- eng_by_ttwa[["Sheffield & Rotherham"]]$eu15_01
eu15_t2 <- eng_by_ttwa[["Sheffield & Rotherham"]]$eu15_11

eu12_t2 <- eng_by_ttwa[["Sheffield & Rotherham"]]$eu12_11
eng_by_ttwa[["Sheffield & Rotherham"]]$prop_eu15_t1 <- prop_eu15_t1 <- eu15_t1 / all_t1
eng_by_ttwa[["Sheffield & Rotherham"]]$prop_eu15_t2 <- prop_eu15_t2 <- eu15_t2 / all_t2
eng_by_ttwa[["Sheffield & Rotherham"]]$prop_eu12_t2 <- prop_eu12_t2 <- eu12_t2 / all_t2

f5 <- classIntervals(c(
  prop_eu15_t1, prop_eu15_t2, 
  prop_eu12_t2
  ),n=5,style="fisher")
f5_brks <- round(f5$brks,digits=3)
f5_brks[1] <- f5_brks[1] - 0.001
f5_brks[6] <- f5_brks[6] + 0.001

eng_by_ttwa[["Sheffield & Rotherham"]] %>% 
  tm_shape(.) + 
  tm_fill(
    col = c("prop_eu15_t1", "prop_eu15_t2", "prop_eu12_t2"), 
    breaks = f5_brks
  )




#### Create the spatial structure objects
W.nb.city <- poly2nb(eng_by_ttwa[["Sheffield & Rotherham"]])
W.list.city <- nb2listw(W.nb.city, style = "B")
W.city <- nb2mat(W.nb.city, style = "B")
n.city <- nrow(W.city)

#### Specify the city centre
#city.centre.hall <- c(435197, 387228)
#city.centre.shopping <- c(435327, 387229)
#city.centre.rail <- c(435834, 386954)
city.centre.middle <- c(435453, 387137)


#### Compute the spatial autocorrelation using Moran's I
moran.mc(x=eng_by_ttwa[["Sheffield & Rotherham"]]@data$prop_eu15_t1, listw=W.list.city, nsim=10000)
moran.mc(x=eng_by_ttwa[["Sheffield & Rotherham"]]@data$prop_eu15_t2, listw=W.list.city, nsim=10000)
moran.mc(x=eng_by_ttwa[["Sheffield & Rotherham"]]@data$prop_eu12_t2, listw=W.list.city, nsim=10000)

#### The temporal dependence
# plot(prop.eu15.2001, prop.eu15.2011, col="red", pch=19, xlab="2001", ylab="2011")
# abline(0,1, col="blue")
# cor.test(prop.eu15.2001,prop.eu15.2011)
# 
eng_by_ttwa[["Sheffield & Rotherham"]]@data %>% 
  ggplot(., aes(x = prop_eu15_t1, y = prop_eu15_t2)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

cor.test(prop_eu15_t2, prop_eu15_t1)

################################
#### Fit the model
################################
#### MCMC quantities
burnin <- 10000
n.sample <- 20000
thin <- 10
n.keep <- (n.sample - burnin)/thin

#### Format the data
Y.mat <- cbind(eu15_t1, eu15_t2)
Y <- as.numeric(t(Y.mat))
N.mat <- cbind(all_t1, all_t2)
N <- as.numeric(t(N.mat))

#### Run the model
model <- binomial.MCARleroux(formula=Y~1, trials=N, W=W.city, burnin=burnin, n.sample=n.sample, thin=thin)
model$summary.results

#### Compute the coordinates and ordering from the city centre
coords <- coordinates(eng_by_ttwa[["Sheffield & Rotherham"]])
dist.cc <- sqrt((coords[,1] - city.centre.middle[1])^2 + (coords[ ,2] - city.centre.middle[2])^2)
dist.order <- order(dist.cc)

#### Compute the global RCI and D
indicators.post <- array(NA, c(n.keep,4))
colnames(indicators.post) <- c("RCI2001", "RCI2011", "D2001", "D2011")

for(i in 1:n.keep)
{
  ## Compute the probability and fitted values for the ith posterior sample
  logit <- model$samples$beta[i, ] + model$samples$phi[i, ]
  prob <- exp(logit) / (1 + exp(logit))
  prob.mat <- matrix(prob, nrow=n.city, byrow=TRUE)
  fitted.mat <- N.mat * prob.mat
  
  ## Compute the RCI for both years
  indicators.post[i, 1] <-RCI(fitted.mat[ ,1], all_t1, dist.order)
  indicators.post[i, 2] <-RCI(fitted.mat[ ,2], all_t2, dist.order)
  
  ## Compute D for both years
  p.2001 <- prob.mat[ ,1]
  p.2001.av <- sum(p.2001 * all_t1) / sum(all_t2)
  indicators.post[i, 3] <- sum(all_t1 * abs(p.2001 - p.2001.av)) / (2 * sum(all_t1) * p.2001.av * (1-p.2001.av))   
  
  p.2011 <- prob.mat[ ,2]
  p.2011.av <- sum(p.2011 * all_t2) / sum(all_t2)
  indicators.post[i, 4] <- sum(all_t2 * abs(p.2011 - p.2011.av)) / (2 * sum(all_t2) * p.2011.av * (1-p.2011.av))   
}

## Summarise the results
## RCI and D in 2001 and 2011 - estimate and 95% Credible Interval
round(apply(indicators.post, 2, quantile, c(0.5, 0.025, 0.975)),3)

## Differences in RCI and D in 2011 - 2001
round(quantile(indicators.post[ ,2] - indicators.post[ ,1], c(0.5, 0.025, 0.975)),3)
round(quantile(indicators.post[ ,4] - indicators.post[ ,3], c(0.5, 0.025, 0.975)),3)


##################################################
#### Compute the local RCI for each area and plot
##################################################
fitted <- matrix(model$fitted.values, nrow=n.city, byrow=TRUE)

#use RCI path to calculate Local RCI for each k
K.range <- seq(10,n.city,10)
K.length <- length(K.range)
RCI.local.K <- array(NA, c(n.city,2,K.length))
rownames(RCI.local.K) <- eng_by_ttwa[["Sheffield & Rotherham"]]$lsoa_01_cd

for(k in 1:K.length) {
  for(i in 1:n.city)
  {
    ## Compute the ordering from the current DZ
    centre <- coords[i, ]
    dist.cc <- sqrt((coords[ ,1] - centre[1])^2 + (coords[ ,2] - centre[2])^2)
    dist.order <- order(dist.cc) 
    
    ## Compute the RCI
    RCI.local.K[i,1,k] <- RCI.path(fitted[ ,1], all_t1, dist.order,K=K.range[k])
    RCI.local.K[i,2,k] <- RCI.path(fitted[ ,2], all_t2, dist.order,K=K.range[k])
  }
  print(k)
}

## summarise the median and 95% interval of local RCI for each area at each year
dim(RCI.local.K)
RCI.local.2001 <- apply(RCI.local.K[,1,],1,quantile,c(0.5,0.025,0.975))
RCI.local.2011 <- apply(RCI.local.K[,2,],1,quantile,c(0.5,0.025,0.975))

RCI.local.2001 <- t(RCI.local.2001)
RCI.local.2011 <- t(RCI.local.2011)
# the number of significant ones
sum(RCI.local.2001[,2]*RCI.local.2001[,3] > 0)
sum(RCI.local.2011[,2]*RCI.local.2011[,3] > 0)

sum(rownames(RCI.local.2001)==rownames(RCI.local.2011))

RCI.local <- cbind(RCI.local.2001,RCI.local.2011)
colnames(RCI.local) <- c("RCI.2001","lower.2001","upper.2001","RCI.2011","lower.2011","upper.2011")
RCI.local <- data.frame(RCI.local)
RCI.local$centralised2001 <- RCI.local$RCI.2001 > 0 & RCI.local$lower.2001*RCI.local$upper.2001 > 0
RCI.local$decetra2001 <- RCI.local$RCI.2001 < 0 & RCI.local$lower.2001*RCI.local$upper.2001 > 0

RCI.local$centralised2011 <- RCI.local$RCI.2011 > 0 & RCI.local$lower.2011*RCI.local$upper.2011 > 0
RCI.local$decetra2011 <- RCI.local$RCI.2011 < 0 & RCI.local$lower.2011*RCI.local$upper.2011 > 0

## map the local RCI
sheffield.map@data <- data.frame(sheffield.map@data,RCI.local)

#year 2001
plot(sheffield.map)
plot(sheffield.map[sheffield.map$centralised2001==1,],col="red",add=TRUE,lwd=1.5)
plot(sheffield.map[sheffield.map$decetra2001==1,],col="green",add=TRUE,lwd=1.5)
#year 2011
plot(sheffield.map)
plot(sheffield.map[sheffield.map$centralised2011==1,],col="red",add=TRUE,lwd=1.5)
plot(sheffield.map[sheffield.map$decetra2011==1,],col="green",add=TRUE,lwd=1.5)

## Map the local RCI
range <- quantile(c(min(sheffield.map$RCI.2001-0.01),max(sheffield.map$RCI.2011+0.01)), seq(0, 1, 0.1))
n.range <- length(range) - 1
spplot(sheffield.map,c("RCI.2001","RCI.2011"), 
       at=range, scales=list(draw=TRUE), xlab="Easting", names.attr=c("Local RCI 2001", "Local RCI 2011"),
       ylab="Northing", col.regions=hsv(0.7, seq(0.2,1,length.out=n.range),1), col="transparent")


writeOGR(sheffield.map,dsn=".",layer="sheffield.map",,driver="ESRI Shapefile")
