# The function takes Landsat Collection 2 Level 2 data as an input, clips it to
# the region of interest, masks out invalid values using QA rasters, and converts
# raw vegetation index values to NDVI.
# For detailed information on the data, please refer to:
# Landsat 8 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-8-collection-2-level-2-science-product-guide
# Landsat 4-7 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-4-7-collection-2-level-2-science-product-guide
# Landsat Surface Reflectance-derived Spectral Indices
# https://www.usgs.gov/landsat-missions/landsat-surface-reflectance-derived-spectral-indices
# Landsat Normalized Difference Vegetation Index
# https://www.usgs.gov/landsat-missions/landsat-normalized-difference-vegetation-index

landsatC2L2NDVI.batch.fun <- function(input_dir  = "F:/rp_2021_uhi_kyiv/data_raw/Landsat_C2L2_veg/",
                                      output_dir = "F:/rp_2021_uhi_kyiv/data_output/Landsat_C2L2_NDVI/",
                                      roi_utm35_file = "F:/rp_2021_uhi_kyiv/data_raw/roi_shapes/buffer2025_1km_32635.shp",
                                      roi_utm36_file = "F:/rp_2021_uhi_kyiv/data_raw/roi_shapes/buffer2025_1km_32636.shp"){
    
    library(raster)
    library(rgdal)
    library(binaryLogic)
    
#---- DEFINE INITIAL INPUTS ----
    
    # define input-output directories
    setwd(input_dir)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = T)
    }
    
    # read clipping files
    roi_utm35 <- readOGR(roi_utm35_file)
    roi_utm36 <- readOGR(roi_utm36_file)
    
    # read Landsat files
    
    # read Landsat 8 UTM35 files
    ndvi_835_file <- list.files(pattern = glob2rx("LC08\\_*\\_182025\\_*\\_SR_NDVI.tif$"))
    qa_835_file   <- list.files(pattern = glob2rx("LC08\\_*\\_182025\\_*\\_QA_PIXEL.tif$"))
    date_835      <- substr(ndvi_835_file, 18, 25)
    
    # read Landsat 8 UTM36 files
    ndvi_836_file <- list.files(pattern = glob2rx("LC08\\_*\\_181025\\_*\\_SR_NDVI.tif$"))
    qa_836_file   <- list.files(pattern = glob2rx("LC08\\_*\\_181025\\_*\\_QA_PIXEL.tif$"))
    date_836      <- substr(ndvi_836_file, 18, 25)
    
    # read Landsat 4-7 UTM35 files
    ndvi_45735_file <- list.files(pattern = glob2rx("LT04\\_*\\_182025\\_*\\_SR_NDVI.tif$|LT05\\_*\\_182025\\_*\\_SR_NDVI.tif$|LE07\\_*\\_182025\\_*\\_SR_NDVI.tif$"))
    qa_45735_file   <- list.files(pattern = glob2rx("LT04\\_*\\_182025\\_*\\_QA_PIXEL.tif$|LT05\\_*\\_182025\\_*\\_QA_PIXEL.tif$|LE07\\_*\\_182025\\_*\\_QA_PIXEL.tif$"))
    date_45735      <- substr(ndvi_45735_file, 18, 25)
    
    # read Landsat 4-7 UTM36 files
    ndvi_45736_file <- list.files(pattern = glob2rx("LT04\\_*\\_181025\\_*\\_SR_NDVI.tif$|LT05\\_*\\_181025\\_*\\_SR_NDVI.tif$|LE07\\_*\\_181025\\_*\\_SR_NDVI.tif$"))
    qa_45736_file   <- list.files(pattern = glob2rx("LT04\\_*\\_181025\\_*\\_QA_PIXEL.tif$|LT05\\_*\\_181025\\_*\\_QA_PIXEL.tif$|LE07\\_*\\_181025\\_*\\_QA_PIXEL.tif$"))
    date_45736      <- substr(ndvi_45736_file, 18, 25)
    
#---- DEFINE CALCULATION FUNCTIONS ----
    
    qa8.fun   <- function(x) {ifelse((x>=21824 & x<=65534), x, NA)}
    qa457.fun <- function(x) {ifelse((x>=5440 & x<=16382), x, NA)}
    ndvi.fun  <- function(x, y) {ifelse(((x<=10000 & x>=-10000) & (y == 1)), x*y*0.0001, NA)}
    
#---- CALCULATE LANDSAT 8 NDVI ----
    
    # process Landsat 8 UTM36 NDVI ----
    
    print("Starting Landsat 8 UTM36 files")
    
    for (i in 1:length(ndvi_836_file)){
        # read rasters
        ndvi_836_raster <- raster(ndvi_836_file[i])
        qa_836_raster   <- raster(qa_836_file[i])
        
        # clip rasters to roi
        ndvi_836_pix <- mask(crop(ndvi_836_raster, extent(roi_utm36)), roi_utm36)
        qa_836_pix   <- mask(crop(qa_836_raster, extent(roi_utm36)), roi_utm36)
        
        # process qa raster to a binary mask
        qa_valid <- calc(qa_836_pix, fun = qa8.fun)
        # decode binary values
        qa_code   <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[1:6]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        ndvi_val <- overlay(ndvi_836_pix, qa_reclass, fun = ndvi.fun)
        
        writeRaster(ndvi_val, filename = paste0(output_dir,"ndvi_", date_836[i], ".tif"), format = "GTiff")
        print(paste("Processed", i, "out of", length(ndvi_836_file)))
        
    }
    
    print("Completed Landsat 8 UTM36 files")

    # process Landsat 8 UTM35 NDVI ----
    
    template_utm36 <- ndvi_val #save template for reprojection
    
    print("Starting Landsat 8 UTM35 files")
    
    for (i in 1:length(ndvi_835_file)){
        # read rasters
        ndvi_835_raster <- raster(ndvi_835_file[i])
        qa_835_raster   <- raster(qa_835_file[i])
        
        # clip rasters to roi
        ndvi_835_pix <- mask(crop(ndvi_835_raster, extent(roi_utm35)), roi_utm35)
        qa_835_pix   <- mask(crop(qa_835_raster, extent(roi_utm35)), roi_utm35)
        
        #process qa raster to a binary mask
        qa_valid <- calc(qa_835_pix, fun = qa8.fun)
        # decode binary values
        qa_code   <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[1:6]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        ndvi_val <- overlay(ndvi_835_pix, qa_reclass, fun = ndvi.fun)
        
        ndvi_val_36 <- projectRaster(from = ndvi_val, to = template_utm36)
        writeRaster(ndvi_val_36, filename = paste0(output_dir,"ndvi_", date_835[i], ".tif"), format = "GTiff")
        print(paste("Processed", i, "out of", length(ndvi_835_file)))
        
    }
    
    print("Completed Landsat 8 UTM35 files")
    
#---- CALCULATE LANDSAT 4-7 NDVI ----
    
    # process Landsat 4-7 UTM36 NDVI ----
    
    print("Starting Landsat 4-7 UTM36 files")
    
    for (i in 1:length(ndvi_45736_file)){
        # read rasters
        ndvi_45736_raster <- raster(ndvi_45736_file[i])
        qa_45736_raster   <- raster(qa_45736_file[i])
        
        # clip rasters to roi
        ndvi_45736_pix <- mask(crop(ndvi_45736_raster, extent(roi_utm36)), roi_utm36)
        qa_45736_pix   <- mask(crop(qa_45736_raster, extent(roi_utm36)), roi_utm36)
        
        #process qa raster to a binary mask
        qa_valid <- calc(qa_45736_pix, fun = qa457.fun)
        # decode binary values
        qa_code   <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[c(1, 2, 4:6)]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        ndvi_val <- overlay(ndvi_45736_pix, qa_reclass, fun = ndvi.fun)
        
        writeRaster(ndvi_val, filename = paste0(output_dir,"ndvi_", date_45736[i], ".tif"), format = "GTiff")
        print(paste("Processed", i, "out of", length(ndvi_45736_file)))
        
    }
    
    print("Completed Landsat 4-7 UTM36 files")
    
    template_utm36 <- ndvi_val #save template for reprojection
    
    # process Landsat 4-7 UTM35 NDVI ----
    
    print("Starting Landsat 4-7 UTM35 files")
    
    for (i in 1:length(ndvi_45735_file)){
        # read rasters
        ndvi_45735_raster <- raster(ndvi_45735_file[i])
        qa_45735_raster   <- raster(qa_45735_file[i])
        
        # clip rasters to roi
        ndvi_45735_pix <- mask(crop(ndvi_45735_raster, extent(roi_utm35)), roi_utm35)
        qa_45735_pix   <- mask(crop(qa_45735_raster, extent(roi_utm35)), roi_utm35)
        
        # process qa raster to a binary mask
        qa_valid <- calc(qa_45735_pix, fun = qa457.fun)
        # decode binary values
        qa_code   <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[c(1, 2, 4:6)]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        ndvi_val <- overlay(ndvi_45735_pix, qa_reclass, fun = ndvi.fun)
        ndvi_val_36 <- projectRaster(from = ndvi_val, to = template_utm36)
        
        writeRaster(ndvi_val_36, filename = paste0(output_dir,"ndvi_", date_45735[i], ".tif"), format = "GTiff")
        print(paste("Processed", i, "out of", length(ndvi_45735_file)))
        
    }
    
    print("Completed Landsat 4-7 UTM36 files")
    
    print("Completed overall processing!")
    
}