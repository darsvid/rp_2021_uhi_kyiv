#https://gis.stackexchange.com/questions/277575/how-to-do-spatial-correlation-between-two-sets-of-rasters-in-r

#example 1 - simple correlation between two rasters using moving window approach
#question - how to define mowing window size
library(raster) 
library(spatialEco)

b <- brick(system.file("external/rlogo.grd", package="raster"))
x <- b[[1]]
y <- b[[3]]
r.cor <- rasterCorrelation(x, y, s = 11, type = "pearson")
plot(r.cor)

#example 2 - Dutilleul's modified t-test. This method corrects the correlation for
#the presence of spatial autocorrelation, functionally addressing pseudo-replication
#in the data and providing a correct p-value

library(gstat)                                         
library(sp)
library(spatialEco)                                            

data(meuse)                                            
data(meuse.grid)                                       
coordinates(meuse) <- ~x + y                           
coordinates(meuse.grid) <- ~x + y                      

# GRID-1 log(copper):                                              
v1 <- variogram(log(copper) ~ 1, meuse)                  
x1 <- fit.variogram(v1, vgm(1, "Sph", 800, 1))           
G1 <- krige(zinc ~ 1, meuse, meuse.grid, x1, nmax = 30)
gridded(G1) <- TRUE                                      
G1@data = as.data.frame(G1@data[,-2])

# GRID-2 log(elev):                                              
v2 <- variogram(log(elev) ~ 1, meuse)                  
x2 <- fit.variogram(v2, vgm(.1, "Sph", 1000, .6))        
G2 <- krige(elev ~ 1, meuse, meuse.grid, x2, nmax = 30)
gridded(G2) <- TRUE    
G2@data <- as.data.frame(G2@data[,-2])
G2@data[,1] <- G2@data[,1]

corr <- raster.modified.ttest(G1, G2)    
plot(raster(corr,1), main="spatially adjusted raster correlation") #use 3 for p.value

# random points
corr.rand <- raster.modified.ttest(G1, G2, sub.sample = TRUE, type = "random")
plot(corr.rand)

# hexagonal 
corr.hex <- raster.modified.ttest(G1, G2, sub.sample = TRUE, d = 500, size = 1000)	
head(corr.hex@data)
bubble(corr.hex, "corr")


#now lets try the same approach with lst and ndvi trends correlation
#https://cran.r-project.org/web/packages/spatialEco/spatialEco.pdf
library(raster)
library(sp)
library(rgdal)
library(spatialEco)

#load rasters
ndvi_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/ndvi_slope_p01_rois/ndvi_slope_p01_01.tif")
lst_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/lst_slope_p01_rois/lst_slope_p01_01.tif")

#convert them to spatial pixels df
#ndvi_spdf <- as(ndvi_ts_r, 'SpatialGridDataFrame')
ndvi_spxdf <- as(ndvi_ts_r, 'SpatialPixelsDataFrame')
lst_spxdf <- as(lst_ts_r, 'SpatialPixelsDataFrame')
#x <- SpatialPixelsDataFrame(points = ndvi_ts_r[c("x", "y")], data = ndvi_ts_r)
#lst_spdf <- as(lst_ts_r, 'SpatialGridDataFrame')

#https://rdrr.io/cran/SpatialPack/src/R/cor.spatial.R
corr.rand <- raster.modified.ttest(ndvi_spxdf, lst_spxdf, sub.sample = T, type = "random")
head(corr.rand)
plot(corr.rand[1])
writeOGR(obj=corr.rand, dsn="F:/2021_lst_vs_veg/raster_test",
         layer="corr_rand", driver="ESRI Shapefile")

corr.hex <- raster.modified.ttest(ndvi_spxdf, lst_spxdf, sub.sample = T,
                                  type = "hexagon")
writeOGR(obj=corr.hex, dsn="F:/2021_lst_vs_veg/raster_test",
         layer="corr_hex", driver="ESRI Shapefile")

#think about applying parallel

corr <- raster.modified.ttest(ndvi_spxdf, lst_spxdf)
plot(raster(corr, 3))

#Error in if (p < 2) stop("'coords' must be a matrix with two columns") : 
    #argument is of length zero

corr.hex <- raster.modified.ttest(ndvi_spxdf, lst_spxdf, d = 1000, sub.sample = T,
                                  type = "hexagon", size = 1000)
head(corr.hex@data)
bubble(corr.hex, "corr")

cond.msk.fun <- function(x, y) {ifelse((!is.na(x) & !is.na(y)), 1, NA)}
both_r <- overlay(lst_ts_r, ndvi_ts_r, fun = cond.msk.fun, forcefun = T)
writeRaster(both_r, "F:/2021_lst_vs_veg/raster_test/lst_ndvi_true.tif", format='GTiff',overwrite=TRUE)

