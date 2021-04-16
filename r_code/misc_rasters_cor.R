library(raster)
#library(sp)
#library(rgdal)
#library(spatialEco)
library(SpatialPack)


#lst_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/lst_slope_p01_rois/lst_slope_p01_01.tif")
#ndvi_ts_r <- raster("F:/2021_lst_vs_veg/test_rois/ndvi_slope_p01_rois/ndvi_slope_p01_01.tif")

lst_ts_r <- raster("F:/2021_lst_vs_veg/raster_test/lst_slope_p01.tif")
ndvi_ts_r <- raster("F:/2021_lst_vs_veg/raster_test/ndvi_slope_p01.tif")
true_ts_r <- raster("F:/2021_lst_vs_veg/raster_test/lst_ndvi_true.tif")

#https://stackoverflow.com/questions/48720606/correlation-between-2-rasters-accounting-for-spatial-autocorrelation

#simple correlation
cor_simp <- cor(getValues(lst_ts_r), getValues(ndvi_ts_r), use = "complete.obs")
#-0.614971969490554 

#modified ttest
coord_df <- as.data.frame(true_ts_r, xy = T, na.rm = T)[, 1:2]
ndvi    <- extract(ndvi_ts_r, coord_df)
lst     <- extract(lst_ts_r, coord_df)
cor_mod <- modified.ttest(x = ndvi, y = lst, coords = coord_df,
                          nclass = NULL) 