# The script takes Landsat Collection 2 Level 2 data as an input, clips it to
# the region of interest, masks out invalid values using QA rasters, and converts
# values to the Celsius degrees.
# For detailed information on the data please refer to:
# Landsat 8 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-8-collection-2-level-2-science-product-guide
# Landsat 4-7 Collection 2 Level 2 Science Product Guide
# https://www.usgs.gov/media/files/landsat-4-7-collection-2-level-2-science-product-guide

# some advice on project organization
# https://richpauloo.github.io/2018-10-17-How-to-keep-your-R-projects-organized/

library(raster)
library(rgdal)
library(binaryLogic)

wd <- list()
wd$landsat <- "F:/rp_2021_uhi_kyiv/data_raw/Landsat_C2L2/"
wd$rois    <- "F:/rp_2021_uhi_kyiv/data_raw/roi_shapes/"

#read clipping files
roi_utm35 <- readOGR(paste0(wd$rois, "buffer2025_1km_32635.shp"))
roi_utm36 <- readOGR(paste0(wd$rois, "buffer2025_1km_32636.shp"))


# read Landsat 8 UTM35 files
pix_qa_835_file <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_182025\\_*\\_QA_PIXEL.TIF$"))
#rad_qa_835_file <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_182025\\_*\\_QA_RADSAT.TIF$"))
t10_835_file    <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_182025\\_*\\_ST_B10.TIF$"))
#tun_835_file    <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_182025\\_*\\_ST_QA.TIF$"))
date_835 <- substr(pix_qa_835_file, 18, 25)

# read Landsat 8 UTM36 files
pix_qa_836_file <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_181025\\_*\\_QA_PIXEL.TIF$"))
#rad_qa_836_file <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_181025\\_*\\_QA_RADSAT.TIF$"))
t10_836_file    <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_181025\\_*\\_ST_B10.TIF$"))
#tun_836_file    <- list.files(wd$landsat, pattern = glob2rx("LC08\\_*\\_181025\\_*\\_ST_QA.TIF$"))
date_836 <- substr(pix_qa_836_file, 18, 25)

# read Landsat 4-7 UTM35 files
pix_qa_45735_file <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_182025\\_*\\_QA_PIXEL.TIF$|LT05\\_*\\_182025\\_*\\_QA_PIXEL.TIF$|LE07\\_*\\_182025\\_*\\_QA_PIXEL.TIF$"))
#rad_qa_45735_file <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_182025\\_*\\_QA_RADSAT.TIF$|LT05\\_*\\_182025\\_*\\_QA_RADSAT.TIF$|LE07\\_*\\_182025\\_*\\_QA_RADSAT.TIF$"))
t6_45735_file     <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_182025\\_*\\_ST_B6.TIF$|LT05\\_*\\_182025\\_*\\_ST_B6.TIF$|LE07\\_*\\_182025\\_*\\_ST_B6.TIF$"))
#tun_45735_file    <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_182025\\_*\\_ST_QA.TIF$|LT05\\_*\\_182025\\_*\\_ST_QA.TIF$|LE07\\_*\\_182025\\_*\\_ST_QA.TIF$"))
date_45735 <- substr(pix_qa_45735_file, 18, 25)

# read Landsat 4-7 UTM36 files
pix_qa_45736_file <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_181025\\_*\\_QA_PIXEL.TIF$|LT05\\_*\\_181025\\_*\\_QA_PIXEL.TIF$|LE07\\_*\\_181025\\_*\\_QA_PIXEL.TIF$"))
#rad_qa_45736_file <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_181025\\_*\\_QA_RADSAT.TIF$|LT05\\_*\\_181025\\_*\\_QA_RADSAT.TIF$|LE07\\_*\\_181025\\_*\\_QA_RADSAT.TIF$"))
t6_45736_file     <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_181025\\_*\\_ST_B6.TIF$|LT05\\_*\\_181025\\_*\\_ST_B6.TIF$|LE07\\_*\\_181025\\_*\\_ST_B6.TIF$"))
#tun_45736_file    <- list.files(wd$landsat, pattern = glob2rx("LT04\\_*\\_181025\\_*\\_ST_QA.TIF$|LT05\\_*\\_181025\\_*\\_ST_QA.TIF$|LE07\\_*\\_181025\\_*\\_ST_QA.TIF$"))
date_45736 <- substr(pix_qa_45736_file, 18, 25)

#clip landsat 8 utm 35 st raster
t10_835_raster <- raster(paste0(wd$landsat, t10_835_file[1]))
pix_t10 <- mask(crop(t10_835_raster, extent(roi_utm35)), roi_utm35)

#clip landsat 8 utm 35 qa raster
pix_qa_raster <- raster(paste0(wd$landsat, pix_qa_835_file[1]))
pix_qa <- mask(crop(pix_qa_raster, extent(roi_utm35)), roi_utm35)

#mask out valid data range 21824 – 65534
pix.qa.fun <- function(x) {ifelse((x>=21824 && x<=65534), x, NA)}
pix_qa_valid <- calc(pix_qa, fun = pix.qa.fun)

#decode binary values
qa_code <- unique(pix_qa_valid)
qa_decode <- vector()

for (i in 1:length(qa_code)){
    x <- as.binary(qa_code[i], littleEndian = T)
    y <- ifelse(any(x[1:6]), NA, 1)
    qa_decode <- c(qa_decode, y)
}

#reclassify raster
rcl_m <- as.matrix(cbind(qa_code, qa_decode), ncol = 2, byrow = T)
qa_reclass <- reclassify(x = pix_qa_valid, rcl = rcl_m, include.lowest = T)

#calculate celsius temperature
lst_c <- ((pix_t10*qa_reclass*0.00341802) + 149) - 273.15



# L4-7

#qa valid values 5440–16382
# check the values - they are slightly other, there is an unused pixel 3
