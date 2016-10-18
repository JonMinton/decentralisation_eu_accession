# Link postcode to MSOA 01 
# ONly do this once as involves some very large files

rm(list = ls())

pacman::p_load(
  readr, dplyr, tidyr,
  stringr
)

postcode_to_msoa2011 <- read_csv("data/big_lookups/National_Statistics_Postcode_Lookup_August_2016_Centroids.csv")

msoa01_to_msoa11 <- read_csv("data/big_lookups/middle_layer_super_output_areas_(2001)_to_middle_layer_super_output_areas_(2011)_to_local_authority_districts_(2011)_e+w_lookup/MSOA01_MSOA11_LAD11_EW_LU.csv")


# Sheffield
# S1 2 HH - Town Hall
# Liverpool
# L1 1HF - St John's Gardens
# Manchester
# M2 5WR - near town hall
# Birmingham
# B3 2EW - near cathedral 
# London
# SQ1P 3JY - Near Palace of Westminster

# Derby
# DE1 3GP

# Bristol
# BS1 3DY - Near Castle Park

# York
# YO1 7JJ - Near York Minster

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^S1[ ]{1,}2HH"
  )) %>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_sheffield

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^L1[ ]{1,}1HF"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_liverpool

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^M2[ ]{1,}5WR"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_manchester

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^B3[ ]{1,}2EW"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_birmingham

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^SW1P3JY"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_london

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^DE1 3GP"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_derby

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^BS1 3DY"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_bristol

postcode_to_msoa2011 %>% 
  filter(str_detect(
    pcd, "^YO1 7JJ"
  ))%>% 
  dplyr::select(pcd, lat, long, lsoa11, msoa11) -> loc_york


# All locations 

loc_all <- bind_rows(
  loc_sheffield,
  loc_liverpool,
  loc_manchester,
  loc_birmingham,
  loc_london,
  loc_derby,
  loc_bristol,
  loc_york
)


msoa01_to_msoa11 %>% 
  dplyr::select(msoa01 = MSOA01CD, 
         msoa11 = MSOA11CD
         ) %>% 
  inner_join(loc_all) -> linked_data

write_csv(linked_data, "data/derived/msoa01_for_centroids.csv")
  
