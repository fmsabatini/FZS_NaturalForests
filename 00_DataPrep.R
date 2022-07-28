library(tidyverse)
library(raster)
library(sf)
library(rgdal)
library(gbm)
library(dismo)
library(blockCV)
library(automap)


## Import AOI - Carpathians
carpathians.sf <- read_sf("../_data/EuropeanMountainAreas/m_massifs_v1.shp") %>% 
  filter(name_mm == "Carpathians") %>% 
  st_set_crs("+init=epsg:3035")

## get bounding box for carpathians
bbox <- st_as_sfc(st_bbox(carpathians.sf))

### Import predictors and crop on bbox
files <- list.files(path="../_tmp", pattern="*.tif$", full.names = T)

### crop rasters and resave
mystack <- raster::stack(files)
names(mystack) <- files0 <- str_replace(names(mystack), pattern="\\.UKR|\\.UKR2|_UKR|UKR", replacement="") #clean up names
mystack <- raster::crop(mystack, as_Spatial(bbox)) ## need to convert to 'sp' object for compatibility reasons

raster::projectRaster(mystack, crs="+proj=longlat +datum=WGS84 +no_defs")

writeRaster(mystack, filename=file.path("../_data/", names(mystack)), bylayer=TRUE, format="GTiff", overwrite=T)








#### 1. Import data on primary forests ####
pgdb <- "../_data/EPFDv2.0_DatabaseOA/EPFD_v2.0.mdb"
fc_list0 <- ogrListLayers(pgdb)
topen <- which(fc_list0=="EU_PrimaryForests_Polygons_OA_v20" )

primary.all <- readOGR(pgdb, fc_list0[topen]) 

writeOGR(obj = primary.all, dsn = "../_data", layer = "EPFD_primaryForest_OA", driver="ESRI Shapefile" )
primary.all2 <- readOGR("../_data/EPFD_primaryForest_OA.shp")

primary.sf <- primary.all %>% 
  st_as_sf() %>% 
  st_transform(crs=crs(mystack)) %>% 
  mutate(is.primary=1) %>% 
  dplyr::select(is.primary)


### filter primary forests in Carpathians
contained <- st_contains(carpathians.sf, primary.sf)
primary.carp <- primary.sf[unique(unlist(contained)),]

### rasterize primary.carp and add to stack
primary.raster <- rasterize(primary.carp, mystack)


mystack <- stack(primary.raster, mystack)
names(mystack)[1] <- "Y.Primary.forest"



mydata0 <- getValues(mystack) %>% 
  as_tibble() 
mydata.coord <- raster::coordinates(mystack) %>% 
  as_tibble()
mydata <- mydata.coord %>% 
  bind_cols(mydata0) %>% 
  mutate(Y.Primary.forest=ifelse(!is.na(Y.Primary.forest), 1, 0)) %>% 
  filter(!is.na(X0.0b.mask.figs.UKR) & !is.na(X3.2.Growing.Stock.UKR)) %>% 
  ### assign nominal predictors to factors
  mutate(X0.4.BiogeoRegions2016 = factor(X0.4.BiogeoRegions2016))


### sample presence and absence data
### using a stratified resample based on a 10 x10 km grid
mydata <- mydata %>% 
  mutate(x.cut=floor(x/10000), 
         y.cut=floor(y/10000)) %>% 
  group_by(x.cut, y.cut) 

## explore distribution 
mydata %>% 
  summarize(n=sum(X0.0b.mask.figs.UKR, na.rm=T), 
            n.prim=sum(Y.Primary.forest)) %>% 
  arrange(desc(n.prim)) %>% 
  filter(n.prim>0)

## sample 2 presence and two absence points per grid (if available)
## create summary for each grid cell
train.data <- mydata %>% 
  left_join(mydata %>% 
              summarize(n=sum(X0.0b.mask.figs.UKR, na.rm=T),
                        n.prim=sum(Y.Primary.forest)), 
            by=c("x.cut", "y.cut")) %>% 
  mutate(sample.size=min(n.prim, 3)) %>% #ifelse(n.prim>3 & n > (n.prim*2), n.prim, 0)) %>% 
  mutate(sample.size=replace(sample.size, 
                             list = Y.Primary.forest ==0 & n>6,
                             values = 6)) %>% 
  filter(n.prim!=0) %>% 
  group_by(Y.Primary.forest, .add=T) %>% 
  sample_n(size = sample.size) %>% 
  ungroup() %>% 
  mutate(Y.Primary.forest=as.logical(Y.Primary.forest)) %>% 
  as.data.frame()

myk <- kfold(train.data, k=5)





## convert back to sf for checking
train.data.sf <- train.data
sp::coordinates(train.data.sf) <- ~x+y
crs(train.data.sf) <- crs(mystack)

train.data.sf <- train.data.sf %>% 
#  st_set_crs(crs(mystack)) %>% 
#  st_transform(crs(mystack)) %>% 
  st_as_sf()
  


##check with blockCV
sac <- spatialAutoRange(rasterLayer = mystack[[-c(1,5,14, 16)]],
                        sampleNumber = 5000,
                        #speciesData=pa_data %>% 
                        #  filter(complete.cases(mydata.world)),
                        doParallel = F,
                        showPlots = TRUE)
autocorrelation.i <- sac$range
print(paste("Autocorrelation is", round(autocorrelation.i/1000), "km"))

#spatial blocking with randomly assigned grid cells, 
# having the of autocorrelation.i
sb2 <- spatialBlock(speciesData = train.data.sf %>% 
                      st_transform(crs(primary.all)), # presence-background data
                    species = "Y.Primary.forest",
                    rasterLayer = projectRaster(mystack[[1]], crs = crs(primary.all)),
                    theRange=autocorrelation.i, ## median of the range of the variograms of the predictors
                    k = 5,
                    iteration=99,
                    selection = "random", 
                    seed=1,
                    progress = T)



plot(carpathians.sf[,1], col=NA)
plot(primary.carp, add=T)
plot(train.data.sf[,1], add=T, col=2, cex=0.6)
## all background points are around the prsence points


#### run BRTs and explore output
brt.primary <- gbm.step(data=train.data %>% 
                          filter(sb2$foldID != 5), gbm.x = 5:19, gbm.y = "Y.Primary.forest",  ##1600 trees
                                family = "bernoulli", tree.complexity = 5, # factors="X2.4a.GensBiome",
                                learning.rate = 0.01, bag.fraction = 0.5)


summary(brt.primary)
gbm.plot(brt.primary, n.plots=6, plot.layout=c(3, 2), write.title = FALSE)

e <- evaluate(train.data %>% 
                filter(sb2$foldID==5) %>% 
                filter(Y.Primary.forest==1), 
              train.data %>% 
                filter(sb2$foldID==5) %>% 
                filter(Y.Primary.forest==0), 
              brt.primary, 
              n.trees=brt.primary$gbm.call$best.trees)


## predict to the whole landscape
myp <- gbm::predict.gbm(brt.primary, 
                           getValues(mystack[[-1]]) %>% 
                             as.data.frame(), 
                           type="response")

myp.df <- raster::coordinates(mystack) %>% 
  as_tibble() %>% 
  bind_cols(brt.pred=myp) %>% 
  bind_cols(getValues(mystack[[c("X0.0b.mask.figs.UKR", "X3.2.Growing.Stock.UKR")]]) %>% 
              as_tibble()) %>% 
  mutate(brt.pred=ifelse( (!is.na(X0.0b.mask.figs.UKR) & !is.na(X3.2.Growing.Stock.UKR)), brt.pred, NA))


myp.raster <- raster(mystack[[1]])
myp.raster <- setValues(myp.raster, myp.df$brt.pred)







## plotting
plot(bbox)
plot(carpathians.sf[,1], col=NA, add=T)
plot(primary.raster, add=T)




