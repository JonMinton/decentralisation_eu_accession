############################### Inference on COB segregation file_ 14/06/2016


#### Sheffield RCI ANALYSES

#setwd("C:/Users/Guanpeng/Dropbox/Migrants in the city confernce papers/segregation measures/analysis")
setwd("/Users/gavin/Dropbox/Migrants in the city confernce papers/segregation measures/analysis")

#### Load the libraries needed
library(Rcpp)
library(shapefiles)
library(sp)
library(spdep)
library(CARBayes)
library(MCMCpack)
library(truncdist)

library(rgdal)
library(maptools)
require(RColorBrewer)
library(dplyr)
#### Source the functions for running the model
source('RCI.R')
source('binomial.MCARleroux.R')
Rcpp::sourceCpp('aqmen.cpp')


#### Check the boundary changes for different TTWAs. if not much not unchanged LSOAs or DATA zones would be used
### Import the English COB data in 2001
cob2001 <- read.csv("Data/England_LSOA_2001.csv",header=TRUE,stringsAsFactors = FALSE)
head(cob2001)
str(cob2001)
summary(cob2001)

### Import the English COB data in 2001
cob2011 <- read.csv("Data/England_LSOA_2011.csv",header=TRUE,stringsAsFactors = FALSE)
head(cob2011)
str(cob2011)
summary(cob2011)

#### For the three cities (London, Manchester and Sheffield), check the LSOA boundary changes during the two census years
# link LSOAs to the TTWA defined in 2007
lsoa2ttwa <- read.csv("Data/LSOA01_TTWA07_UK_LU.csv",header = TRUE, stringsAsFactors = FALSE)
match.ind <- match(cob2001$code,lsoa2ttwa$LSOA01CD)
sum(is.na(match.ind))
cob2001$ttwa <- lsoa2ttwa$TTWA07NM[match.ind]
data.frame(table(cob2001$ttwa))

# changes of LSOAs boundaries
changes.boundary <- read.csv("Data/LSOA01_LSOA11_EW_LUv2.csv",header = TRUE,stringsAsFactors = FALSE) %>% tbl_df() %>%
  select(LSOA01CD,LSOA11CD,CHGIND)

# then link the COB data
#2001
cob2001tbl <- cob2001 %>% tbl_df() %>%
  mutate(EU15_2001=Republic_of_Ireland + Other_EU_countries,
         UK_2001=England + Scotland + Wales + Northern_Ireland,
         all_2001=All_people,LSOA01CD=code) %>%
  select(Name,LSOA01CD,ttwa,all_2001,UK_2001,EU15_2001)

#2011
cob2011tbl <- cob2011 %>% tbl_df() %>%
  mutate(EU15_2011=Ireland + EU_2001,EU12_2011=EU_2001_TO_2011,
         UK_2011=UK,all_2011=all,LSOA11CD=code) %>%
  select(LSOA11CD,all_2011,UK_2011,EU15_2011,EU12_2011,EU12_2011)

# link to the data of boundary changes
changes.boundary.1 <- left_join(changes.boundary,cob2001tbl,by="LSOA01CD") 
changes.boundary.1 <- left_join(changes.boundary.1,cob2011tbl,by="LSOA11CD")
  
### we use LSOA boundary in 2001 as baseline geographical units
# case 1; split--one lsoa was split into two or more lsoa in 2011. in this case, aggregate the COB data in 2011 based on 2001 lsoas
changes.boundary.1 %>% filter(CHGIND == "S")
#grepl("*.2011",names(changes.boundary.2))

data.split <- changes.boundary.1 %>% filter(CHGIND == "S") %>%
  select(LSOA01CD,all_2011,UK_2011,EU15_2011,EU12_2011) %>%
  group_by(LSOA01CD) %>%
  summarise(all_2011=sum(all_2011),UK_2011=sum(UK_2011),EU15_2011=sum(EU15_2011),EU12_2011=sum(EU12_2011))

# link back to the
s.ind <- changes.boundary.1$CHGIND =="S"
match.ind <- match(changes.boundary.1$LSOA01CD[s.ind],data.split$LSOA01CD)
changes.boundary.1[s.ind,c("all_2011","UK_2011","EU15_2011","EU12_2011")] <- data.split[match.ind,2:5]

# case 2: two or more LSOAs were merged into a single LSOA in 2011. in this case we need aggregate the 2001 LSOAs accordingly.
changes.boundary.1 %>% filter(CHGIND == "M")

data.merge <- changes.boundary.1 %>% filter(CHGIND == "M") %>%
  select(LSOA11CD,all_2001,UK_2001,EU15_2001) %>%
  group_by(LSOA11CD) %>%
  summarise(all_2001=sum(all_2001),UK_2001=sum(UK_2001),EU15_2001=sum(EU15_2001))

# link back
m.ind <- changes.boundary.1$CHGIND =="M"
match.ind <- match(changes.boundary.1$LSOA11CD[m.ind],data.merge$LSOA11CD)
changes.boundary.1[m.ind,names(data.merge)[-1]] <- data.merge[match.ind,2:4]

changes.boundary.1  %>%
  filter(ttwa %in% c("London","Manchester","Sheffield & Rotherham")) %>%
  group_by(ttwa,CHGIND) %>%
  summarise(num=length(ttwa))

# select unique ids of LSOAs in 2001
COB.data <- changes.boundary.1[!duplicated(changes.boundary.1$LSOA01CD),] %>% data.frame

# drop complex boundary changes, they only account for a very small proportion of the data
#COB.data <- COB.data[!COB.data$CHGIND == "X",]
COB.data %>% tbl_df() %>% filter(ttwa %in% c("London","Manchester","Sheffield & Rotherham")) %>%
  filter(CHGIND == "X") %>%
  group_by(ttwa) %>%
  summarise(num=length(ttwa))

COB.data <- COB.data[!COB.data$CHGIND == "X",]
################################################################################################

######################
###### Sheffield #####
######################

sheffield.map <- readOGR(dsn=".",layer="Sheffield_ttwa",stringsAsFactors = FALSE)
sheffield.map$lsoa <- as.character(sheffield.map$LSOA04CD)
# cob data
sheffield.data <- COB.data[COB.data$ttwa == "Sheffield & Rotherham",]
head(sheffield.data)
sheffield.data[sheffield.data$CHGIND == "M",]
length(unique(sheffield.data$LSOA11CD))

# first link cob to the spatial data, then dissolve the polygons
match.ind <- match(sheffield.map$lsoa,sheffield.data$LSOA01CD)
# drop the two polygons
sheffield.map <- sheffield.map[!is.na(match.ind),]
sheffield.map@data <- data.frame(sheffield.map@data,sheffield.data[match.ind[!is.na(match.ind)],])
head(sheffield.map@data)

#### Now dissolve the Sheffield ttwa polygons based on the ids of LSOAs in 2011
library(rgeos)
if (rgeosStatus()) {
  temp.map <- unionSpatialPolygons(sheffield.map, IDs = as.character(sheffield.map$LSOA11CD)) 
}
if(rgeosStatus()) {
  temp.df <- as(sheffield.map@data,"data.frame")[!duplicated(sheffield.map$LSOA11CD),]
  row.names(temp.df) <- temp.df$LSOA11CD
  sheffield.map.final <- SpatialPolygonsDataFrame(temp.map,temp.df)   
}


save.image("segregation14_06.RData")
load("segregation14_06.RData")

##### Mapping the EU migrants distributions
head(sheffield.map.final@data)
all_2001 <- sheffield.map.final$all_2001
all_2011 <- sheffield.map.final$all_2011
UK_2001 <- sheffield.map.final$UK_2001
UK_2011 <- sheffield.map.final$UK_2011
EU15_2001 <- sheffield.map.final$EU15_2001
EU15_2011 <- sheffield.map.final$EU15_2011
EU12_2011 <- sheffield.map.final$EU12_2011
sheffield.map.final$prop.eu15.2001 <- prop.eu15.2001 <- EU15_2001 / all_2001
sheffield.map.final$prop.eu15.2011 <- prop.eu15.2011 <- EU15_2011 / all_2011
sheffield.map.final$prop.eu12.2011 <-prop.eu12.2011<- EU12_2011 / all_2011


library(classInt)
f5 <- classIntervals(c(prop.eu15.2001,prop.eu15.2011,prop.eu12.2011),n=5,style="fisher")
f5.brks <- round(f5$brks,digits=3)
f5.brks[1] <- f5.brks[1] - 0.001
f5.brks[6] <- f5.brks[6] + 0.001

xx <- sheffield.map.final@bbox
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(),
                 offset = c(min(xx[1,] + 1000),xx[2,1] + 1000), scale = 10000, fill=c("transparent","black"),which=3)
centre.sp <- SpatialPoints(matrix(city.centre.middle,nrow=1),proj4string = CRS(proj4string(sheffield.map.final)))
city.centre = list("sp.points", centre.sp, col="black", pch=3, cex=0.8)

text1 <- list("sp.text", c(min(xx[1,]) + 1000,xx[2,1] + 2000), "0",which=3,cex=0.8)
text2 <- list("sp.text", c(min(xx[1,])+13000,xx[2,1] + 2000), "10 km",which=3,cex=0.8)
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(),
                   offset = c(min(xx[1,])+2000,xx[2,1] + 3000), scale = 3000,which=3)


tiff(filename="results/sheffieldCOB.tif",units="in",width = 8, height = 5, res=300)

spplot(sheffield.map.final, c("prop.eu15.2001","prop.eu15.2011","prop.eu12.2011"), sp.layout=list(city.centre,scalebar,text1,text2,northarrow),
       at=f5.brks, col.regions=brewer.pal(5,"Reds"), 
       col="transparent",pretty=TRUE,as.table=TRUE,layout=c(2,2),skip=c(FALSE,FALSE,TRUE,FALSE),
       colorkey=list(labels=list(at=f5.brks)))

dev.off()

#### Create the spatial structure objects
W.nb.city <- poly2nb(sheffield.map.final)
W.list.city <- nb2listw(W.nb.city, style = "B")
W.city <- nb2mat(W.nb.city, style = "B")
n.city <- nrow(W.city)

#### Specify the city centre
city.centre.hall <- c(435197, 387228)
city.centre.shopping <- c(435327, 387229)
city.centre.rail <- c(435834, 386954)
city.centre.middle <- c(435453, 387137)


#### Compute the spatial autocorrelation using Moran's I
moran.mc(x=sheffield.map.final@data$prop.eu15.2001, listw=W.list.city, nsim=10000)
moran.mc(x=sheffield.map.final@data$prop.eu15.2011, listw=W.list.city, nsim=10000)
moran.mc(x=sheffield.map.final@data$prop.eu12.2011, listw=W.list.city, nsim=10000)

#### The temporal dependence
plot(prop.eu15.2001, prop.eu15.2011, col="red", pch=19, xlab="2001", ylab="2011")
abline(0,1, col="blue")
cor.test(prop.eu15.2001,prop.eu15.2011)




################################
#### Fit the model
################################
#### MCMC quantities
burnin <- 10000
n.sample <- 20000
thin <- 10
n.keep <- (n.sample - burnin)/thin

#### Format the data
Y.mat <- cbind(EU15_2001, EU15_2011)
Y <- as.numeric(t(Y.mat))
N.mat <- cbind(all_2001, all_2011)
N <- as.numeric(t(N.mat))

#### Run the model
model <- binomial.MCARleroux(formula=Y~1, trials=N, W=W.city, burnin=burnin, n.sample=n.sample, thin=thin)
model$summary.results

#### Compute the coordinates and ordering from the city centre
coords <- coordinates(sheffield.map.final)
dist.cc <- sqrt((coords[1] - city.centre.middle[1])^2 + (coords[ ,2] - city.centre.middle[2])^2)
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
  indicators.post[i, 1] <-RCI(fitted.mat[ ,1], all_2001, dist.order)
  indicators.post[i, 2] <-RCI(fitted.mat[ ,2], all_2011, dist.order)
  
  ## Compute D for both years
  p.2001 <- prob.mat[ ,1]
  p.2001.av <- sum(p.2001 * all_2001) / sum(all_2001)
  indicators.post[i, 3] <- sum(all_2001 * abs(p.2001 - p.2001.av)) / (2 * sum(all_2001) * p.2001.av * (1-p.2001.av))   
  
  p.2011 <- prob.mat[ ,2]
  p.2011.av <- sum(p.2011 * all_2011) / sum(all_2011)
  indicators.post[i, 4] <- sum(all_2011 * abs(p.2011 - p.2011.av)) / (2 * sum(all_2011) * p.2011.av * (1-p.2011.av))   
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
rownames(RCI.local.K) <- sheffield.map$lsoa

for(k in 1:K.length) {
  for(i in 1:n.city)
  {
    ## Compute the ordering from the current DZ
    centre <- coords[i, ]
    dist.cc <- sqrt((coords[ ,1] - centre[1])^2 + (coords[ ,2] - centre[2])^2)
    dist.order <- order(dist.cc) 
    
    ## Compute the RCI
    RCI.local.K[i,1,k] <- RCI.path(fitted[ ,1], total2001, dist.order,K=K.range[k])
    RCI.local.K[i,2,k] <- RCI.path(fitted[ ,2], total2011, dist.order,K=K.range[k])
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
