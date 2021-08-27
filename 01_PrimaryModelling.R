## ---- message=F, warning=F---------------------------------------------------------------------------------------------------------------------------------
# install.packages(c("tidyverse", "raster", "st", "rgdal", "gbm", "dismo", "blockCV", "automap")) if necessary
library(tidyverse)
library(raster)
library(sf)
library(rgdal)
library(gbm)
library(dismo)
library(blockCV)
library(automap)


## ---- message=F, warning=F, fig.align="center", fig.height=8, fig.width=10---------------------------------------------------------------------------------
files <- list.files(path="../_data", pattern="*.tif$", full.names = T)
mystack <- raster::stack(files)
mystack <- projectRaster(mystack, crs=crs("+proj=longlat +datum=WGS84 +no_defs"), method = "ngb")
plot(mystack)


## ---- message=F, warning=F, cache=T------------------------------------------------------------------------------------------------------------------------
primary.all <- readOGR("../_data/EPFD_primaryForest_OA.shp") 
#visualize data
glimpse(primary.all@data)
#filter out all unnecessary fields
primary.sf <- primary.all %>% 
  st_as_sf() %>% 
  mutate(is.primary=1) %>% 
  dplyr::select(is.primary)


## ---- warning=F--------------------------------------------------------------------------------------------------------------------------------------------
carpathians.sf <- read_sf("../_data/EuropeanMountainAreas/m_massifs_v1.shp") %>% 
  filter(name_mm == "Carpathians") %>% 
  st_transform(crs(primary.sf))



## ---- message=F, warning=F, fig.align="center", fig.height=5, fig.width=6----------------------------------------------------------------------------------
contained <- st_contains(carpathians.sf, primary.sf)
primary.carp <- primary.sf[unique(unlist(contained)),]

#visualize
ggplot()  + 
  geom_sf(data=carpathians.sf, col=NA, fill=4, alpha=0.5) + 
  geom_sf(data=primary.carp) + 
  theme_bw()


## ---- message=F, warning=F---------------------------------------------------------------------------------------------------------------------------------
primary.raster <- rasterize(primary.carp, mystack)
mystack <- stack(primary.raster, mystack)
names(mystack)[1] <- "Y.Primary.forest"


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
mydata0 <- getValues(mystack) %>% 
  as_tibble() 
mydata.coord <- raster::coordinates(mystack) %>% 
  as_tibble()
mydata <- mydata.coord %>% 
  bind_cols(mydata0) %>% 
  mutate(Y.Primary.forest=ifelse(!is.na(Y.Primary.forest), 1, 0)) %>% 
  ### mask only forested pixels. Exclude pixels with NA in GS
  filter(!is.na(X0.0b.mask) & !is.na(X3.2.Growing.Stock)) %>% 
  ### define biogeoregions as a factors
  mutate(X0.4.BiogeoRegions2016 = factor(X0.4.BiogeoRegions2016))

#visualize
glimpse(mydata)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
mysample.naive <- mydata %>% 
  #count num of PF pixels
  mutate(sample.size=nrow(mydata %>% filter(Y.Primary.forest==1))) %>%  
  #set sample size for pixels to be randomly drawn. 
  #We keep all primary forests, and sample 10 times as many background points
  mutate(sample.size=ifelse(Y.Primary.forest==1, sample.size, 10*sample.size)) %>% 
  group_by(Y.Primary.forest) %>%
  sample_n(sample.size) %>% 
  ungroup() %>% 
  # split in k=5 folds
  mutate(kfold=kfold(Y.Primary.forest, k=5)) %>% 
  mutate(Y.Primary.forest=as.logical(Y.Primary.forest)) %>% 
  as.data.frame()



## ----fig.align="center", fig.height=8, fig.width=8---------------------------------------------------------------------------------------------------------
train.data <- mysample.naive %>% 
  filter(kfold != 5)
test.data <- mysample.naive %>% 
  filter(kfold == 5)

#### run BRTs and explore output
brt.naive <- gbm.step(data=train.data, 
                      gbm.x = 5:19, gbm.y = "Y.Primary.forest",  
                      family = "bernoulli", tree.complexity = 5, 
                      learning.rate = 0.05, bag.fraction = 0.5)


summary.gbm(brt.naive, plotit = F)
gbm.plot(brt.naive, n.plots=6, plot.layout=c(3, 2), write.title = FALSE)


## ----fig.align="center", fig.height=8, fig.width=8---------------------------------------------------------------------------------------------------------
# evaluate our model based on the test dataset.
(e <- dismo::evaluate(p = test.data %>% 
                filter(Y.Primary.forest==1), 
              a = test.data %>% 
                filter(Y.Primary.forest==0), 
              brt.naive, 
              n.trees=brt.naive$gbm.call$best.trees))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------
## predict to the whole landscape
myp <- gbm::predict.gbm(brt.naive, 
                           getValues(mystack[[-1]]) %>% 
                             as.data.frame(), 
                           type="response")

myp.raster.naive <- raster(mystack[[1]])
myp.raster.naive <- setValues(myp.raster.naive, myp)
myp.raster.naive <- mask(myp.raster.naive, mask = mystack$X0.0b.mask)

plot(myp.raster.naive)


## ---- fig.height=4, fig.width=5----------------------------------------------------------------------------------------------------------------------------
mysample.naive$predicted <- gbm::predict.gbm(brt.naive, mysample.naive, type="response")
boxplot(predicted ~ Y.Primary.forest, data=mysample.naive)



## ---- cache=T, warning=F-----------------------------------------------------------------------------------------------------------------------------------
# make mysample.naive a spatial 'sf' object
mysample.naive.sf <- mysample.naive
coordinates(mysample.naive.sf) <- ~x+y
crs(mysample.naive.sf) <- crs(mystack)
mysample.naive.sf <- mysample.naive.sf %>% 
  st_as_sf() 


## Calculate spatial aucotorrelation of continuous predictors
sac <- spatialAutoRange(rasterLayer = mystack[[-c(1,5,14, 16)]],
                        sampleNumber = 5000,
                        doParallel = F,
                        showPlots = TRUE)
autocorrelation.i <- sac$range
print(paste("Autocorrelation is", round(autocorrelation.i/1000), "km"))


## ---- cache=T, warning=F-----------------------------------------------------------------------------------------------------------------------------------
#spatial blocking with randomly assigned grid cells, 
sb2 <- spatialBlock(speciesData = mysample.naive.sf, 
                    species = "Y.Primary.forest",
                    rasterLayer = mystack[[2]],
                    theRange=autocorrelation.i, 
                    k = 5,
                    iteration=49,
                    selection = "random", 
                    seed=100,
                    progress = T)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
train.data.aware <- mysample.naive %>% 
  filter(sb2$foldID != 5) 
test.data.aware <- mysample.naive %>% 
  filter(sb2$foldID == 5) 


#### run BRTs and explore output
brt.aware <- gbm.step(data=train.data.aware,
                        gbm.x = 5:19, gbm.y = "Y.Primary.forest",  
                        family = "bernoulli", tree.complexity = 5, 
                        learning.rate = 0.05, bag.fraction = 0.5)

(e <- evaluate(test.data.aware %>% 
                filter(Y.Primary.forest==1), 
              test.data.aware %>% 
                filter(Y.Primary.forest==0), 
              brt.aware, 
              n.trees=brt.aware$gbm.call$best.trees))

