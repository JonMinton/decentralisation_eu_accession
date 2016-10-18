# 27/7/2016

rm(list = ls())

# Script to compare distribution of EU  and EU-Accession populations 
# from 2001 to 2011 at English & Welsh LSOAs

require(pacman)


pacman::p_load(
  stringr, tidyr, 
  shapefiles, sp,
  spdep, rgeos, maptools,
  tmap, classInt,
  ggplot2, 
  RColorBrewer,
  CARBayes, MCMCpack,
  truncdist,
  readr, readxl,
  dplyr, purrr,
  forcats
  
)

#### Source the functions for running the model
source("scripts/RCI.R")
source("scripts/binomial.MCARleroux.R")
Rcpp::sourceCpp("scripts/aqmen.cpp")



dta <- read_csv("data/derived/cob_both_censuses_simplified.csv")

# Link msoa to ttwa 
ttwa_lookup <- read_csv("data/LSOA01_TTWA01_UK_LU.csv")
lookup <- read_csv("data/2001_geographical_lookup/OA01_LSOA01_MSOA01_EW_LU.csv")
ttwa_lookup %>% 
  dplyr::select(LSOA01CD, TTWA01CD, TTWA01NM) %>% 
  inner_join(lookup) %>% 
  dplyr::select(lsoa = LSOA01CD, msoa = MSOA01CD, ttwa = TTWA01CD, ttwa_name = TTWA01NM) -> ttwa_msoa_lookup


#Names of TTWAs to use 

ttwas_to_use <- c(
  "Sheffield and Rotherham",
  "Liverpool", 
  "Manchester",
  "Birmingham", 
  "London", 
  "Derby",  
  "Bristol",  
  "York"  
)

# TTWA centroid MSOA finder 

ttwa_msoa_centroids <- read_csv("data/derived/msoa01_for_centroids.csv") %>% 
  mutate(ttwa = c("London", "Manchester", "Liverpool", "Sheffield and Rotherham", "Birmingham", "York", "Derby", "Bristol")) %>% 
  dplyr::select(ttwa, msoa = msoa01)


# Tasks: 

# For each TTWA
# produce reduced map with only relevant TTWA MSOAs 


# Shapefile for MSOA 

shp_msoa <- read_shape(file = "shapefiles/Middle_layer_super_output_areas_(E+W)_2001_Boundaries_(Full_Extent)_V2/MSOA_2001_EW_BFE_V2.shp")


# Simple version of task: do just for Sheffield 

ttwa_msoa_lookup

get_shapefile_for_ttwa <- function(ttwa_nm, shp, dta_ttwa){
  dta_ttwa  %>% 
    filter(ttwa_name == ttwa_nm)  %>% 
    dplyr::select(msoa)  %>% 
    unique()   %>% 
    .$msoa -> msoas_of_interest
  shp2 <- shp[shp$MSOA01CD %in% msoas_of_interest,]
  shp2
}

get_data_for_ttwa <- function(ttwa, dta, dta_ttwa){
  dta_ttwa[dta_ttwa$ttwa_name == ttwa,] %>% 
    dplyr::select(msoa) %>% 
    unique() %>% 
    .$msoa -> msoas_of_interest
  dta %>% 
    filter(msoa %in% msoas_of_interest) -> out
  out  
}

calc_distance_to_centre <- function(shp, msoa_centre){
  
  gCentroid(shp, byid = T) %>% 
    as(., "data.frame") -> tmp
  
  data_frame(msoa = as.character(shp@data$MSOA01CD), x = tmp$x, y = tmp$y) %>% 
    mutate(centre = msoa == msoa_centre) %>% 
    mutate(distance_to_centre = ((x - x[centre])^2 + (y - y[centre])^2)^0.5) -> output
  output
}

calc_distance_prop_relationship <- function(dta, dist_cent){
  dta  %>% 
    group_by(census, msoa)   %>% 
    mutate(prop = count / sum(count)) -> dta_prop
  dist_cent %>% dplyr::select(msoa, distance_to_centre) %>% 
    right_join(dta_prop) -> out
  out
}

# distance and proportion data for all ttwas 

dta_shp_ttwa <- data_frame(
  ttwa = ttwas_to_use,
  msoa_centroid = c(
    "E02001641","E02001379","E02001060","E02001885","E02000979","E02002808","E02003043","E02002784"
  ),
  data = map(ttwa, get_data_for_ttwa, dta =dta , dta_ttwa = ttwa_msoa_lookup),
  shp_for_ttwa = map(.x = ttwa, .f = get_shapefile_for_ttwa, shp = shp_msoa, dta_ttwa = ttwa_msoa_lookup),
  distance_to_centroid = map2(shp_for_ttwa, msoa_centroid, calc_distance_to_centre),
  dta_prop_dist = map2(data, distance_to_centroid, calc_distance_prop_relationship)
)


# Flatten above to dataframe
ttwa_dist_props <- dta_shp_ttwa %>% 
  dplyr::select(ttwa, dta_prop_dist) %>% 
  unnest()

ttwa_dist_props %>% 
  group_by(ttwa) %>% 
  mutate(ttwa_size = sum(count)) %>% 
  ungroup() %>% 
  mutate(ttwa = ifelse(ttwa == "Sheffield and Rotherham", "Sheffield", ttwa)) %>% 
  mutate(ttwa = fct_reorder(ttwa, ttwa_size)) %>% 
  mutate(cob = fct_relevel(cob, "UK", "rWE", "EE", "Other")) %>% 
  mutate(census = factor(census)) %>% 
  mutate(distance_to_centre = distance_to_centre / 1000) %>% 
  ggplot(., aes(x = distance_to_centre, y = prop, group = census, colour = census)) +
  geom_point(shape = ".", alpha = 0.2) + stat_smooth(se = F) + 
  facet_grid(cob ~ ttwa, scales = "free") +
  labs(x = "Distance to centre (km)", y = "Proportion of population in this category")

ggsave("figures/proportion_distance_plots.png", height = 20, width = 30, units = "cm", dpi = 300)

# Decile plots 

ttwa_dist_props %>% 
  group_by(ttwa) %>% 
  mutate(ttwa_size = sum(count)) %>% 
  ungroup() %>% 
  mutate(ttwa = ifelse(ttwa == "Sheffield and Rotherham", "Sheffield", ttwa)) %>% 
  mutate(ttwa = fct_reorder(ttwa, ttwa_size)) %>% 
  mutate(cob = fct_relevel(cob, "UK", "rWE", "EE", "Other")) %>% 
  mutate(census = factor(census)) %>% 
  group_by(ttwa) %>% 
  mutate(dist_decile = factor(ntile(distance_to_centre, 10))) %>% 
  group_by(ttwa, dist_decile, cob, census) %>% 
  summarise(mean_prop = mean(prop), sd_prop = sd(prop), n = sum(count), se_prop = sd_prop / (n ^ 0.5)) %>% 
  ggplot(., aes(x = dist_decile, y = mean_prop, group = census, colour = census)) + 
  geom_line() + geom_point() + 
  geom_linerange(aes(ymax = mean_prop + 2 * se_prop, ymin = mean_prop - 2 * se_prop)) +
  facet_grid(cob ~ ttwa, scales = "free") +
  labs(x = "Decile from centre (1 = nearest)", y = "Proportion of group in decile")
ggsave("figures/proportion_decile_distance_plots.png", height = 20, width = 30, units = "cm", dpi = 300)





# Code for doing RCI calculation  - stuck due to areas with no neighbours 



################################
#### Fit the model
################################
#### MCMC quantities
burnin <- 10000
n.sample <- 20000
thin <- 10
n.keep <- (n.sample - burnin)/thin

# For each ttwa, 
# for each year 
# for each group 

# Start off with Sheffield & Rotherham, Eastern European, 2011 

dta_shp_ttwa %>% 
  filter(ttwa == "Sheffield and Rotherham") %>% 
  .$shp_for_ttwa %>% .[[1]] -> this_ttwa

dta_shp_ttwa %>% 
  filter(ttwa == "Sheffield and Rotherham") %>% 
  .$data %>% .[[1]] -> this_data

this_data %>% 
  filter(census == 2001) %>% 
  group_by(msoa) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  filter(cob == "EE") -> this_data_ss

this_ttwa %>% 
  append_data(
    shp = ., data = this_data_ss,
    key.shp = "MSOA01CD", key.data = "msoa",
    ignore.na = T
  ) %>% 
  .[!is.na(.$count),] -> this_dta_shp

W_nb <- poly2nb(this_dta_shp)
W_list <- nb2listw(W_nb, style = "B") 
W <- nb2mat(W_nb, style = "B")
n <- nrow(W)

Y <- as.integer(this_dta_shp@data$count)

N <- as.integer(this_dta_shp@data$total)

model_01 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  

# This doesn't work as some areas have no neighbours


dta_shp_ttwa %>% 
  filter(ttwa == "Derby") %>% 
  .$shp_for_ttwa %>% .[[1]] -> this_ttwa

dta_shp_ttwa %>% 
  filter(ttwa == "Derby") %>% 
  .$data %>% .[[1]] -> this_data

this_data %>% 
  filter(census == 2001) %>% 
  group_by(msoa) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  filter(cob == "EE") -> this_data_ss

this_ttwa %>% 
  append_data(
    shp = ., data = this_data_ss,
    key.shp = "MSOA01CD", key.data = "msoa",
    ignore.na = T
  ) %>% 
  .[!is.na(.$count),] -> this_dta_shp

W_nb <- poly2nb(this_dta_shp)
W_list <- nb2listw(W_nb, style = "B") 
W <- nb2mat(W_nb, style = "B")
n <- nrow(W)

Y <- as.integer(this_dta_shp@data$count)

N <- as.integer(this_dta_shp@data$total)

model_01 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  



# Now to revise do_model to run all available years 

do_model <- function(place, initial_dist = 18000){
  
  this_dist <- initial_dist
  has_islands <- T
  
  while(has_islands){
    dz_2011  %>% 
      append_data(
        shp = . , data = centre_distance, 
        key.shp = "DataZone", key.data = "datazone"
      ) %>% 
      .[.$nearest_centre == place,] %>% 
      .[.$distance_to_centre <= this_dist,] -> dz_city
    # plot(dz_city, main = this_dist)
    # browser()
    
    w <-  poly2nb(dz_city)   
    if (any(card(w) == 0)){
      this_dist <- this_dist + 1000
    } else {
      has_islands <- F
    }
  }
  
  simd_2011_reweighted %>% 
    dplyr::select(dz_2011, year, pop_total, pop_incomedeprived) -> tmp
  
  tmp %>% 
    dplyr::select(-pop_incomedeprived) %>% 
    mutate(year = paste0("pop_total_", year)) %>% 
    spread(year, pop_total) -> pops
  
  tmp %>% 
    dplyr::select(-pop_total) %>% 
    mutate(year = paste0("pop_id_", year)) %>% 
    spread(year, pop_incomedeprived) -> incdeps
  
  popinc <- pops %>% inner_join(incdeps)
  rm(tmp, pops, incdeps)
  
  
  dz_city %>% 
    append_data(
      shp = ., data = popinc,
      key.shp = "DataZone", key.data = "dz_2011",
      ignore.na = T
    ) %>% 
    .[!is.na(.$pop_total_2004),] -> dz_city
  
  
  W_nb <- poly2nb(dz_city)
  # distance to allow the 'islands' to be included 
  W_list <- nb2listw(W_nb, style = "B")
  W <- nb2mat(W_nb, style = "B")
  n <- nrow(W)
  
  # Need to remove unconnected datazones ('islands')
  
  
  
  #### Format the data
  Y.mat <- cbind(
    as.integer(dz_city@data$pop_id_2004),
    as.integer(dz_city@data$pop_id_2006), 
    as.integer(dz_city@data$pop_id_2009), 
    as.integer(dz_city@data$pop_id_2012),
    as.integer(dz_city@data$pop_id_2016) 
  )
  
  N.mat <- cbind(
    as.integer(dz_city@data$pop_total_2004),
    as.integer(dz_city@data$pop_total_2006), 
    as.integer(dz_city@data$pop_total_2009), 
    as.integer(dz_city@data$pop_total_2012),
    as.integer(dz_city@data$pop_total_2016) 
  )
  
  denom_too_small <- N.mat[,5] < 1
  N.mat[denom_too_small, 5] <- 1
  
  Y.mat[denom_too_small, 5] <- 0
  
  
  
  #### Run the model
  Y <- Y.mat[,1] ; N <- N.mat[,1]
  model_01 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  
  
  Y <- Y.mat[,2] ; N <- N.mat[,2]
  model_02 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  
  
  Y <- Y.mat[,3] ; N <- N.mat[,3]
  model_03 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  
  
  Y <- Y.mat[,4] ; N <- N.mat[,4]
  model_04 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  
  
  #placeholder for new data 
  Y <- Y.mat[,5] ; N <- N.mat[,5]
  model_05 <- S.CARleroux(Y ~ 1, trials = N, W = W, family = "binomial", burnin = burnin, n.sample = n.sample, thin = thin )  
  
  #  model$summary.results
  
  #### Compute the coordinates and ordering from the city centre
  dist.order <- order(dz_city@data$distance_to_centre)
  
  #### Compute the global RCI and D
  indicators.post <- array(NA, c(n.keep,10))
  colnames(indicators.post) <- c(
    "RCI_2004", "RCI_2006","RCI_2009","RCI_2012", "RCI_2016", 
    "D_2004", "D_2006", "D_2009", "D_2012", "D_2016"
  )
  
  all_2004 <- as.integer(N.mat[,1])
  all_2006 <- as.integer(N.mat[,2])
  all_2009 <- as.integer(N.mat[,3])
  all_2012 <- as.integer(N.mat[,4])
  all_2016 <- as.integer(N.mat[,5]) # placeholder for new data 
  
  
  for(i in 1:n.keep)
  {
    ## Compute the probability and fitted values for the ith posterior sample
    logit_01 <- model_01$samples$beta[i, ] + model_01$samples$phi[i, ]
    prob_01 <- exp(logit_01) / (1 + exp(logit_01))
    prob.mat_01 <- matrix(prob_01, nrow=n, byrow=TRUE)
    fitted.mat_01 <- N.mat[,1] * prob.mat_01
    # Explore from here 
    
    logit_02 <- model_02$samples$beta[i, ] + model_02$samples$phi[i, ]
    prob_02 <- exp(logit_02) / (1 + exp(logit_02))
    prob.mat_02 <- matrix(prob_02, nrow=n, byrow=TRUE)
    fitted.mat_02 <- N.mat[,2] * prob.mat_02
    
    
    logit_03 <- model_03$samples$beta[i, ] + model_03$samples$phi[i, ]
    prob_03 <- exp(logit_03) / (1 + exp(logit_03))
    prob.mat_03 <- matrix(prob_03, nrow=n, byrow=TRUE)
    fitted.mat_03 <- N.mat[,3] * prob.mat_03
    
    logit_04 <- model_04$samples$beta[i, ] + model_04$samples$phi[i, ]
    prob_04 <- exp(logit_04) / (1 + exp(logit_04))
    prob.mat_04 <- matrix(prob_04, nrow=n, byrow=TRUE)
    fitted.mat_04 <- N.mat[,4] * prob.mat_04
    
    #Placeholder for new data 
    logit_05 <- model_05$samples$beta[i, ] + model_05$samples$phi[i, ]
    prob_05 <- exp(logit_05) / (1 + exp(logit_05))
    prob.mat_05 <- matrix(prob_05, nrow=n, byrow=TRUE)
    fitted.mat_05 <- N.mat[,5] * prob.mat_05
    
    ## Compute the RCI for both years
    indicators.post[i, 1] <- RCI(fitted.mat_01, as.integer(N.mat[,1]), dist.order)
    indicators.post[i, 2] <- RCI(fitted.mat_02, as.integer(N.mat[,2]), dist.order)
    indicators.post[i, 3] <- RCI(fitted.mat_03, as.integer(N.mat[,3]), dist.order)
    indicators.post[i, 4] <- RCI(fitted.mat_04, as.integer(N.mat[,4]), dist.order)
    indicators.post[i, 5] <- RCI(fitted.mat_05, as.integer(N.mat[,5]), dist.order) # placeholder
    
    ## Compute D for both years
    p_2004 <- prob.mat_01
    p_2004_av <- sum(p_2004 * all_2004) / sum(all_2004)
    indicators.post[i, 6] <- sum(all_2004 * abs(p_2004 - p_2004_av)) / (2 * sum(all_2004) * p_2004_av * (1-p_2004_av))   
    
    p_2006 <- prob.mat_02
    p_2006_av <- sum(p_2006 * all_2006) / sum(all_2006)
    indicators.post[i, 7] <- sum(all_2006 * abs(p_2006 - p_2006_av)) / (2 * sum(all_2006) * p_2006_av * (1-p_2006_av))   
    
    p_2009 <- prob.mat_03
    p_2009_av <- sum(p_2009 * all_2009) / sum(all_2009)
    indicators.post[i, 8] <- sum(all_2009 * abs(p_2009 - p_2009_av)) / (2 * sum(all_2009) * p_2009_av * (1-p_2009_av))   
    
    p_2012 <- prob.mat_04
    p_2012_av <- sum(p_2012 * all_2012) / sum(all_2012)
    indicators.post[i, 9] <- sum(all_2012 * abs(p_2012 - p_2012_av)) / (2 * sum(all_2012) * p_2012_av * (1-p_2012_av))   
    
    p_2016 <- prob.mat_04
    p_2016_av <- sum(p_2016 * all_2016) / sum(all_2016)
    indicators.post[i, 10] <- sum(all_2016 * abs(p_2016 - p_2016_av)) / (2 * sum(all_2016) * p_2016_av * (1-p_2016_av))   
    
  }
  
  indicators.post  
}




indicators_aberdeen <- do_model("Aberdeen")
indicators_dundee <- do_model("Dundee")
indicators_edinburgh <- do_model("Edinburgh")
indicators_glasgow <- do_model("Glasgow")

# Tidy all indicator draws 

tidy_indicators <- function(mtrx, place){
  mtrx %>% 
    data.frame() %>% 
    tbl_df %>% 
    mutate(draw = 1:dim(mtrx)[1]) %>% 
    gather(mp, value, -draw) %>% 
    mutate(place = place) %>% 
    separate(mp, into = c("measure", "period")) %>% 
    dplyr::select(place, measure, period, draw, value)
}

tind_aberdeen <- tidy_indicators(indicators_aberdeen, "Aberdeen")
tind_dundee <- tidy_indicators(indicators_dundee, "Dundee")
tind_edinburgh <- tidy_indicators(indicators_edinburgh, "Edinburgh")
tind_glasgow <- tidy_indicators(indicators_glasgow, "Glasgow")

tidy_indicators_all <- reduce(
  list(tind_aberdeen, tind_dundee, tind_edinburgh, tind_glasgow),
  bind_rows
)


write_csv(tidy_indicators_all, path = "data/all_posterior_draws.csv")





# Older material ----------------------------------------------------------




# 2011 centroids
ttwa_centroids <- c(
  Aberdeen = "S01006646",
  Glasgow = "S01010265",
  Edinburgh = "S01008677",
  Dundee = "S01007705",
  Inverness = "S01010620",
  Perth = "S01011939",
  `Falkirk and Stirling` = "S01013067"
)

# using rgeos::gCentroid

calc_distance_to_centres <- function(shp, code_centre){
  gCentroid(shp, byid = T) %>% 
    as(., "data.frame") -> tmp
  
  data_frame(dz = as.character(shp@data$DataZone), x = tmp$x, y = tmp$y) %>% 
    mutate(centre = dz == code_centre) %>% 
    mutate(distance_to_centre = ((x - x[centre])^2 + (y - y[centre])^2)^0.5) -> output
  output
}


fn <- function(val, nm){
  dz_2011 %>% 
    calc_distance_to_centres(., val) %>% 
    .[c(1, 5)] -> out 
  names(out) <- c("datazone", nm)
  out
}

# Find nearest centre and distance to nearest centre
map2(ttwa_centroids, names(ttwa_centroids), fn) %>% 
  reduce(., inner_join) %>% 
  gather(place, distance, -datazone) %>%
  arrange(datazone) %>% 
  group_by(datazone) %>% 
  mutate(min_distance = min(distance)) %>%
  filter(distance == min_distance) %>% 
  ungroup() %>% 
  transmute(datazone, nearest_centre = place, distance_to_centre = distance) -> centre_distance
# TTWA lookup 



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
debug(create_ttwa_tmap)
ttwa_msoa_lookup %>% select(msoa, ttwa_name) %>% filter(!duplicated(.)) -> tm_short
dta %>% 
  left_join(tm_short) %>% 
  filter(!is.na(census), !is.na(ttwa_name)) %>% 
  group_by(msoa, census) %>% 
  mutate(proportion = count / sum(count)) %>%
  group_by(ttwa_name, census, cob) %>% 
  nest() %>% 
  slice(1:10) %>% 
  mutate(ttwa_shp = map(data, create_ttwa_shp, shp = shp_msoa)) %>% 
  mutate(title_name = paste0(ttwa_name, "_" , census, "_" , cob)) -> tmp 
tmp %>%   
  mutate(tmp = walk2(ttwa_shp, title_name, create_ttwa_tmap))


create_change_ttwa_tmap <- function(shp, title){
  tm_shape(shp) + 
    tm_fill("change", palette = "Paired", style = "quantile", n = 10) + 
    tm_legend(legend.outside = T, legend.outside.position = "right") -> tm
  
  save_tmap(tm, filename = paste0("tmaps/change_over_censuses/",title, ".png"), width = 20, height = 20, units = "cm", dpi = 300)
  NULL
}


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
