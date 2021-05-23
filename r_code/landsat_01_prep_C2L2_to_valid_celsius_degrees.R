# The function takes Landsat Collection 2 Level 2 data as an input, clips it to
# the region of interest, masks out invalid values using QA rasters, and converts
# raw surface temperature values to the land surface temperatures in Celsius degrees.
# For detailed information on the data please refer to:
# Landsat 8 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-8-collection-2-level-2-science-product-guide
# Landsat 4-7 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-4-7-collection-2-level-2-science-product-guide

landsatC2L2LST.batch.fun <- function(input_dir      = "F:/rp_2021_uhi_kyiv/data_raw/Landsat_C2L2/",
                                     output_dir     = "F:/rp_2021_uhi_kyiv/data_output/Landsat_C2L2_LST/",
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
    
    #read clipping files
    roi_utm35 <- readOGR(roi_utm35_file)
    roi_utm36 <- readOGR(roi_utm36_file)
    
    #read Landsat files
    
    # read Landsat 8 UTM36 files
    qa_836_file  <- list.files(pattern = glob2rx("LC08\\_*\\_181025\\_*\\_QA_PIXEL.TIF$"))
    t10_836_file <- list.files(pattern = glob2rx("LC08\\_*\\_181025\\_*\\_ST_B10.TIF$"))
    date_836     <- substr(qa_836_file, 18, 25)
    
    # read Landsat 8 UTM35 files
    qa_835_file  <- list.files(pattern = glob2rx("LC08\\_*\\_182025\\_*\\_QA_PIXEL.TIF$"))
    t10_835_file <- list.files(pattern = glob2rx("LC08\\_*\\_182025\\_*\\_ST_B10.TIF$"))
    date_835     <- substr(qa_835_file, 18, 25)
    
    # read Landsat 4-7 UTM36 files
    qa_45736_file <- list.files(pattern = glob2rx("LT04\\_*\\_181025\\_*\\_QA_PIXEL.TIF$|LT05\\_*\\_181025\\_*\\_QA_PIXEL.TIF$|LE07\\_*\\_181025\\_*\\_QA_PIXEL.TIF$"))
    t6_45736_file <- list.files(pattern = glob2rx("LT04\\_*\\_181025\\_*\\_ST_B6.TIF$|LT05\\_*\\_181025\\_*\\_ST_B6.TIF$|LE07\\_*\\_181025\\_*\\_ST_B6.TIF$"))
    date_45736    <- substr(qa_45736_file, 18, 25)
    
    # read Landsat 4-7 UTM35 files
    qa_45735_file <- list.files(pattern = glob2rx("LT04\\_*\\_182025\\_*\\_QA_PIXEL.TIF$|LT05\\_*\\_182025\\_*\\_QA_PIXEL.TIF$|LE07\\_*\\_182025\\_*\\_QA_PIXEL.TIF$"))
    t6_45735_file <- list.files(pattern = glob2rx("LT04\\_*\\_182025\\_*\\_ST_B6.TIF$|LT05\\_*\\_182025\\_*\\_ST_B6.TIF$|LE07\\_*\\_182025\\_*\\_ST_B6.TIF$"))
    date_45735    <- substr(qa_45735_file, 18, 25)
    
    
    
    
   #---- DEFINE CALCULATION FUNCTIONS ----
    
    qa8.fun   <- function(x) {ifelse((x>=21824 && x<=65534), x, NA)}
    qa457.fun <- function(x) {ifelse((x>=5440 && x<=16382), x, NA)}
    lst.fun   <- function(x, y) {((x*y*0.00341802) + 149) - 273.15}
    
    #---- CALCULTATE LANDSTAT 8 LST ----
    
    #process Landsat 8 UTM36
    
    print("Starting to calculate Landsat 8 UTM36 files")
    
    for (i in 1:length(qa_836_file)){
        # read rasters
        t10_836_raster <- raster(t10_836_file[i])
        qa_836_raster  <- raster(qa_836_file[i])
        
        #clip rasters to roi
        t10_836_pix  <- mask(crop(t10_836_raster, extent(roi_utm36)), roi_utm36)
        qa_836_pix   <- mask(crop(qa_836_raster, extent(roi_utm36)), roi_utm36)
        
        
        #process qa raster to a binary mask
        qa_valid <- calc(qa_836_pix, fun = qa8.fun)
        #decode binary values
        qa_code   <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[1:6]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        lst_c <- overlay(t10_836_pix, qa_reclass, fun = lst.fun)
        
        writeRaster(lst_c, filename = paste0(output_dir,"lst_", date_836[i], ".tif"), format = "GTiff")
        print(paste("Processed", i, "out of", length(qa_836_file)))
        
    }
    
    print("Completed to calculate Landsat 8 UTM36 files")
    
    #process Landsat 8 UTM35
    
    template_utm36 <- lst_c #save template for reprojection
    
    print("Starting to calculate Landsat 8 UTM35 files")
    
    for (i in 1:length(qa_835_file)){
        # read rasters
        t10_835_raster <- raster(t10_835_file[i])
        qa_835_raster  <- raster(qa_835_file[i])
        
        #clip rasters to roi
        t10_835_pix  <- mask(crop(t10_835_raster, extent(roi_utm35)), roi_utm35)
        qa_835_pix   <- mask(crop(qa_835_raster, extent(roi_utm35)), roi_utm35)
        
        
        #process qa raster to a binary mask
        qa_valid <- calc(qa_835_pix, fun = qa8.fun)
        #decode binary values
        qa_code <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[1:6]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m      <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        
        lst_c    <- overlay(t10_835_pix, qa_reclass, fun = lst.fun)
        lst_c_36 <- projectRaster(from = lst_c, to = template_utm36)
        writeRaster(lst_c_36, filename = paste0(output_dir,"lst_", date_835[i], ".tif"), format = "GTiff")
        
        print(paste("Processed", i, "out of", length(qa_835_file)))
        
    }
    
    print("Completed to calculate Landsat 8 UTM35 files")
    
    
    #---- CALCULTATE LANDSTAT 4-7 LST ----
    
    #process Landsat 4-7 UTM36
    
    print("Starting to calculate Landsat 4-7 UTM36 files")
    
    for (i in 1:length(qa_45736_file)){
        # read files
        t6_45736_raster <- raster(t6_45736_file[i])
        qa_45736_raster <- raster(qa_45736_file[i])
        
        #clip rasters to roi
        t6_45736_pix  <- mask(crop(t6_45736_raster, extent(roi_utm36)), roi_utm36)
        qa_45736_pix  <- mask(crop(qa_45736_raster, extent(roi_utm36)), roi_utm36)
        #process qa raster to a binary mask
        qa_valid <- calc(qa_45736_pix, fun = qa457.fun)
        #decode binary values
        qa_code <- unique(qa_valid)
        qa_decode <- vector()
        for (n in 1:length(qa_code)){
            x <- as.binary(qa_code[n], littleEndian = T)
            y <- ifelse(any(x[c(1, 2, 4:6)]), NA, 1)
            qa_decode <- c(qa_decode, y)
        }
        rcl_m <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
        qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
        
        lst_c <- overlay(t6_45736_pix, qa_reclass, fun = lst.fun)
        writeRaster(lst_c, filename = paste0(output_dir,"lst_", date_45736[i], ".tif"), format = "GTiff")
        
        print(paste("Processed", i, "out of", length(qa_45736_file)))
        
    }
    
    print("Completed to calculate Landsat 4-7 UTM36 files")
    
    #process Landsat 4-7 UTM35
    
    template_utm36 <- lst_c #save template for reprojection
    
    print("Starting to calculate Landsat 4-7 UTM35 files")
    
    for (i in 1:length(qa_45735_file)){
      # read files
      t6_45735_raster <- raster(t6_45735_file[i])
      qa_45735_raster <- raster(qa_45735_file[i])
      
      #clip rasters to roi
      t6_45735_pix  <- mask(crop(t6_45735_raster, extent(roi_utm35)), roi_utm35)
      qa_45735_pix  <- mask(crop(qa_45735_raster, extent(roi_utm35)), roi_utm35)
      
      #process qa raster to a binary mask
      qa_valid <- calc(qa_45735_pix, fun = qa457.fun)
      #decode binary values
      qa_code <- unique(qa_valid)
      qa_decode <- vector()
      for (n in 1:length(qa_code)){
        x <- as.binary(qa_code[n], littleEndian = T)
        y <- ifelse(any(x[c(1, 2, 4:6)]), NA, 1)
        qa_decode <- c(qa_decode, y)
      }
      rcl_m <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
      qa_reclass <- reclassify(x = qa_valid, rcl = rcl_m, include.lowest = T)
      
      lst_c <- overlay(t6_45735_pix, qa_reclass, fun = lst.fun)
      lst_c_36 <- projectRaster(from = lst_c, to = template_utm36)
      writeRaster(lst_c_36, filename = paste0(output_dir,"lst_", date_45735[i], ".tif"), format = "GTiff")
      
      print(paste("Processed", i, "out of", length(qa_45735_file)))
      
    }
    print("Completed to calculate Landsat 4-7 UTM35 files")
    
    print("Completed overall processing!")
  
}