library(raster)
library(sp)
library(rgdal)
library(spatialEco)

#load rasters
ndvi_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/ndvi_slope_p01_rois/ndvi_slope_p01_06.tif")
lst_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/lst_slope_p01_rois/lst_slope_p01_06.tif")

#convert them to spatial pixels df
ndvi_spxdf <- as(ndvi_ts_r, 'SpatialPixelsDataFrame')
lst_spxdf <- as(lst_ts_r, 'SpatialPixelsDataFrame')

corr <- raster.modified.ttest(ndvi_spxdf, lst_spxdf)
plot(raster(corr, 3)) #1 - corr, 3 - pvalue

corr_r <- raster(corr, 1)
writeRaster(corr_r, "F:/2021_lst_vs_veg/test_rois/correlations_p01/corr_06.tif",
            format='GTiff',overwrite=TRUE)

pval_r <- raster(corr, 3)
writeRaster(pval_r, "F:/2021_lst_vs_veg/test_rois/correlations_p01/pval_06.tif",
            format='GTiff',overwrite=TRUE)

cond.msk.fun <- function(x, y) {ifelse((y > 0.01), NA, x)}
corr_pval <- overlay(corr_r, pval_r, fun = cond.msk.fun, forcefun = T)
writeRaster(corr_pval, "F:/2021_lst_vs_veg/test_rois/correlations_p01/corr_pval001_06.tif",
            format='GTiff',overwrite=TRUE)
