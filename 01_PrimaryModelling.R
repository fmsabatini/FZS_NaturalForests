## ----message=F, warning=F---------------------------------------------------------------------
# install.packages(c("tidyverse", "raster", "st", "rgdal", "gbm", "dismo", "blockCV", "automap", "cowplot", "rgeos")) if necessary
library(tidyverse)
library(sf)
library(terra)

#library(rgeos)


## ----message=F, warning=F, fig.align="center", fig.height=8, fig.width=10, cache=T------------
files <- list.files(path="../_data", pattern="*.tif$", full.names = T)
mystack <- rast(files)
mystack <- project(mystack, y="epsg:4326", method = "near")
plot(mystack)


## ----message=F, warning=F, cache=T------------------------------------------------------------
primary.all <- read_sf("../_data/EPFD_primaryForest_OA.shp") 
#visualize data
primary.all
#filter out all unnecessary fields
primary.sf <- primary.all %>% 
  mutate(is.primary=1) %>% 
  dplyr::select(is.primary) %>% 
  st_make_valid() #Fixes some topological errors in the input shp
#Drop remaining invalid polygons (n=3)
primary.sf <- primary.sf %>% 
  slice(which(st_is_valid(primary.sf)))


## ----warning=F--------------------------------------------------------------------------------
carpathians.sf <- read_sf("../_data/EuropeanMountainAreas/m_massifs_v1.shp") %>% 
  filter(name_mm == "Carpathians") %>% 
  st_transform(crs(primary.sf)) %>% 
  st_make_valid() 



## ----message=F, warning=F, fig.align="center", fig.height=5, fig.width=6----------------------
contained <- st_contains(carpathians.sf, primary.sf)
primary.carp <- primary.sf[unique(unlist(contained)),]

#visualize
ggplot()  + 
  geom_sf(data=carpathians.sf, col=NA, fill=4, alpha=0.5) + 
  geom_sf(data=primary.carp) + 
  theme_bw()


## ----message=F, warning=F---------------------------------------------------------------------
primary.raster <- rasterize(primary.carp, mystack)
mystack <- c(primary.raster, mystack)
names(mystack)[1] <- "Y.Primary.forest"


## ---------------------------------------------------------------------------------------------
mydata0 <- values(mystack) %>% 
  as_tibble() 
mydata.coord <- crds(mystack, na.rm=F) %>% 
  as_tibble()
mydata <- mydata.coord %>% 
  bind_cols(mydata0) %>% 
  ### mask only forested pixels. Exclude pixels with NA in GS
  filter(!is.na(X0.0b.mask) & !is.na(X3.2.Growing.Stock)) %>% 
  ## Define all remaining forest pixels to NON-primary if NA
  mutate(Y.Primary.forest=ifelse(!is.na(Y.Primary.forest), 1, 0)) %>% 
  ### define biogeoregions as a factors
  mutate(X0.4.BiogeoRegions2016 = factor(X0.4.BiogeoRegions2016))

#visualize
glimpse(mydata)


## ---------------------------------------------------------------------------------------------
library(gbm)
library(dismo)

set.seed(2907)
#count num of PF pixels
n.pf <- nrow(mydata %>% filter(Y.Primary.forest==1))
mysize  <- 0.5 ## reduce sample size to increase speed of model 

mysample.naive <- mydata %>% 
  #define sample.size
  #set sample size for pixels to be randomly drawn. 
  #We keep all primary forests, and sample 10 times as many background points
  mutate(sample.size=ifelse(Y.Primary.forest==1, n.pf*mysize, 10*n.pf*mysize)) %>% 
  group_by(Y.Primary.forest) %>%
  sample_n(size = first(sample.size)) %>% 
  ungroup() %>% 
  # split in k=5 folds
  mutate(kfold=kfold(Y.Primary.forest, k=5)) %>% 
  as.data.frame()



## ----fig.align="center", fig.height=8, fig.width=8--------------------------------------------
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


## ----fig.align="center", fig.height=8, fig.width=8--------------------------------------------
# evaluate our model based on the test dataset.
(e <- dismo::evaluate(p = test.data %>% 
                filter(Y.Primary.forest==1), 
              a = test.data %>% 
                filter(Y.Primary.forest==0), 
              brt.naive, 
              n.trees=brt.naive$gbm.call$best.trees))



## ---------------------------------------------------------------------------------------------
## predict to the whole landscape
myp <- gbm::predict.gbm(brt.naive, 
                           values(mystack[[-1]], na.rm=T) %>% 
                             as.data.frame(), 
                           type="response")
myp.coords <- crds(mystack[[-1]], na.rm=T) 


# Create a new raster using the original as a template, and replace non-na values with predictions based on their coordinates
myp.raster.naive <- mystack[[1]]
cell_indices <- cellFromXY(myp.raster.naive, myp.coords)
myp.raster.naive[cell_indices] <- myp

## mask non forest pixels
myp.raster.naive <- mask(myp.raster.naive, mask = mystack$X0.0b.mask)

plot(myp.raster.naive)


## ----fig.height=4, fig.width=5----------------------------------------------------------------
mysample.naive$predicted <- gbm::predict.gbm(brt.naive, mysample.naive, type="response")
boxplot(predicted ~ Y.Primary.forest, data=mysample.naive)



## ----cache=T, warning=F-----------------------------------------------------------------------
library(blockCV)
library(automap)
library(cowplot)

# make mysample.naive a spatial 'sf' object
mysample.naive.sf <- mysample.naive %>% 
  st_as_sf(coords=c("x","y"),
           crs=crs(mystack))

## Calculate spatial aucotorrelation of continuous predictors
## Takes quite long-time to run!
## suggest to skip and replace with 
## autocorrelation.i <- 138000 #120 km
sac <- cv_spatial_autocor(r = mystack[[-c(5,14, 16)]],
                        column = "Y.Primary.forest",
                        num_sample = 5000,
                        plot = TRUE)
autocorrelation.i <- sac$range
print(paste("Autocorrelation is", round(autocorrelation.i/1000), "km"))


## ----cache=T, warning=F-----------------------------------------------------------------------
#spatial blocking with randomly assigned grid cells, 
sb2 <- cv_spatial(x = mysample.naive.sf, 
                    column = "Y.Primary.forest",
                    r = mystack[[2]],
                    size=autocorrelation.i, 
                    k = 5,
                    iteration=49,
                    selection = "random", 
                    seed=100,
                    progress = T)




## ---------------------------------------------------------------------------------------------
train.data.aware <- mysample.naive %>% 
  filter(sb2$folds_ids != 5) 
test.data.aware <- mysample.naive %>% 
  filter(sb2$folds_ids == 5) 


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


## ---------------------------------------------------------------------------------------------
sessionInfo()


## ----echo=F, eval=F---------------------------------------------------------------------------
## knitr::purl("01_PrimaryModelling.Rmd")

