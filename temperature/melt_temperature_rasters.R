rm(list=ls())

library(raster)
library(data.table)
library(zoo)
library(parallel)
library(Rfast)
library(sf)

source("FILEPATH/get_ids.R")


INTERACTIVE = commandArgs()[2]=='--interactive'


if (INTERACTIVE == FALSE) { # If launched from cmd line use args
  args <- commandArgs(trailingOnly = T) 
  
  loc_id <- args[1]
  year <- args[2]
  release <- args[3]
  n_draws <- args[4]
  outDir <- args[5]
  job.name <- args[6]
  debug <- FALSE
  
} else { # otherwise, this is an interactive run for testing/development and hardcode arguments
  loc_id <- 349
  year <-  2020
  release <- 9
  n_draws <- 1000  
  outDir <- "FILEPATH_REPLACED_FOR_SECURITY"
  job.name <- paste0("melt_", loc_id, "_", year)
  debug <- FALSE # make true for debugging mode that will retain data to help debugging and development
}



max_draw <- n_draws - 1
release_name <- get_ids('release')[release_id == release, release_name]


message(paste0("Location   = ", loc_id))
message(paste0("Year       = ", year))
message(paste0("Output dir = ", outDir))
message(paste0("Job name   = ", job.name))
message(paste0("Draws      = ", max_draw))

# Establish directories
era_dir <- 'FILEPATH_REPLACED_FOR_SECURITY'
pop_dir <- 'FILEPATH_REPLACED_FOR_SECURITY'
shp_dir <- file.path("FILEPATH_REPLACED_FOR_SECURITY", gsub(" ", "_", release_name), "master/shapefiles/")



### READ IN SHAPEFILE ###
shp_version <- paste0(gsub(" ", "", release_name), "_analysis_final")
shape <- read_sf(dsn = shp_dir, layer = shp_version)
message("Shapefile loaded")


# Subset shapefile to target location, unless the location is global (location_id 1)
if (loc_id==1) {
  sub_shape <- shape
} else {
  sub_shape <- shape[shape$loc_id == loc_id, ] 
}
message("Shapefile subsetted")



### POPULATION AND TEMPERATURE DATA ARE NOT AVAILABLE FOR ALL YEARS -- FIND AND SUB WITH CLOSEST YEAR IF NEEDED ###
temp_files <- list.files(file.path(era_dir, 'temp'), 'era5_temp_daily_[0-9]{4}.gri')
temp_file_years <- as.integer(regmatches(temp_files, regexpr('[0-9]{4}', temp_files)))
latest_temp_year <- max(temp_file_years) 
tempYear <- min(year, latest_temp_year)


pop_files <- list.files(pop_dir, 'world_pop_temp[0-9]{4}.csv')
pop_file_years <- as.integer(regmatches(pop_files, regexpr('[0-9]{4}', pop_files)))
latest_pop_year <- max(pop_file_years) 
popYear <- min(year, latest_pop_year)


### CHECK DATE OF TEMP FILE TO ENSURE IT WAS CREATED AFTER THE YEAR IT REPRESENTS ###
temp_file_mtime <- file.info(file.path(era_dir, 'temp', paste0("era5_temp_daily_", tempYear, ".gri")))$mtime
earliest_possible <- as.Date(paste0(tempYear+1, '-01-01'))

if (temp_file_mtime < earliest_possible) {
  warning('Temperature file created before the time it represents!')
} else {
  message('Temperature file modified after end of year it represents')
}



### READ IN POP, TEMP, TEMP SD, & TEMP ZONE DATA ###
start <- Sys.time()

pop           <- fread(file.path(pop_dir, paste0('world_pop_temp', popYear, '.csv')))
dailyTempCat  <- brick(file.path(era_dir, 'temp', paste0('era5_temp_daily_', tempYear, '.gri')), overwrite = TRUE, format = "raster")
sd            <- brick(file.path(era_dir, 'temp_spread', paste0('spread_daily_', tempYear, '.gri')), overwrite = TRUE, format = "raster")
meanTempCat   <- raster(file.path(era_dir, "annual_mean.gri"))
message("Rasters loaded")


# CONVERT POP FROM XY TO RASTER #
pop <- rotate(rasterFromXYZ(pop[,c("longitude","latitude","population")]))
message("Population converted to raster")


# ENSURE ALL RASTERS HAVE SAME PROJECTION #
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(pop)   <- projection 
crs(dailyTempCat)  <- projection
crs(sd)    <- projection
crs(meanTempCat) <- projection
message("Rasters reprojected")


# ENSURE ALL RASTERS HAVE SAME PIXEL SIZE AND ALIGNMENT #
pop <- raster::resample(pop, dailyTempCat, method = "bilinear")
sd <- raster::resample(sd, dailyTempCat, method = "bilinear")
meanTempCat <- raster::resample(meanTempCat, dailyTempCat, method = "bilinear")
message("Rasters resampled")

(dataReadTime <- difftime(Sys.time(), start))


### DEFINE FUNCTION TO EXTRACT THE DATA FROM PIXELS IN LOCATION AND FORMAT FOR MERGE ###
extract2shape <- function(x, shape, melt = FALSE, mask = TRUE, crop = TRUE) {
  tmp <- eval(parse(text = x)) # x is a character containing the name of the object -- getting the object itself here
  if (crop==TRUE) {
    tmp <- crop(tmp, shape)
  } 
  
  if (mask==TRUE) {
    tmp <- raster::mask(tmp, shape)
  }
  
  tmp <- as.data.table(rasterToPoints(tmp, spatial = TRUE))
  
  if (melt==TRUE) {
    tmp <- melt(data = tmp, id.vars = c("x","y"))
    setnames(tmp, "variable", "layer")
    setnames(tmp, "value", x)
  } else {
    names(tmp)[1] <- x
  }
  
  tmp[, `:=` (x = round(x, digits = 2), y = round(y, digits = 2))]
  
  return(tmp)
}

message("Extract function defined")  

### EXTRACT DATA FROM RASTERS AND CREATE DATA TABLE ###  
start <- Sys.time()

# FOR VERY SMALL LOCATIONS WE MAY LOSE ALL PIXELS IF WE MASK TO THE EXACT LOCATION POLYGON, FOR THESE LOCS WE'LL BUFFER THE POLYGON #
# TEST FOR THAT HERE
if (loc_id==1) {
  crop <- F
  mask <- F

} else {
  crop <- T
  mask <- T
  
  testPixels <- raster::mask(crop(pop, sub_shape), sub_shape)
  message("Test pixels extracted")

  testVals <- as.numeric(values(testPixels))
  message("Test values extracted")
  

  while (sum(testVals, na.rm = T)==0) {
    message("Zero population found in the polygon. Buffering.")
    sub_shape <- buffer(sub_shape, 1)
    testPixels <- raster::mask(crop(pop, sub_shape), sub_shape)
    testVals <- as.numeric(values(testPixels))
  }

  rm(testVals, testPixels)
}

message("Pixel test completed")



# EXTRACT DATA #  
for (x in c("pop", "dailyTempCat", "sd", "meanTempCat")) {
  message(x)
  if (x %in% c("dailyTempCat", "sd")) {
    melt = T
  } else {
    melt = F
  }
  assign(x, extract2shape(x, sub_shape, melt = melt, mask = mask))
}
message("Rasters extracted to data tables")


### CLEAN AND FIND UNIQUE VALUES OF TEMPERATURE ZONE ###
meanTempCat[, meanTempCat  := as.integer(round(meanTempCat-273.15, digits = 0))]
zones <- unique(meanTempCat$meanTempCat)


### MERGE COMPONENTS ###
master <- merge(pop, meanTempCat,  by=c("x","y"), all = TRUE)
if (debug==F) rm(pop, meanTempCat)
master <- merge(master, dailyTempCat, by=c("x","y"), all = TRUE)
if (debug==F) rm(dailyTempCat)
master <- merge(master, sd, by=c("x", "y", "layer"), all = TRUE)
if (debug==F) rm(sd)

message("Components merged")

### WARN IF THERE ARE MISSING VALUES OF ANY NEEDED VARIABLES ###    
if (sum(is.na(master$pop))>0) message("Population values are missing for some pixels")
if (sum(is.na(master$dailyTempCat))>0) message("Daily temperature values are missing for some pixels")
if (sum(is.na(master$meanTempCat))>0) message("Temperature zone values are missing for some pixels")
if (sum(is.na(master$sd))>0) message("Temperature SD values are missing for some pixels")


### ENSURE THAT POPULATION VALUES ARE NOT ALL NA OR ZERO (E.G. MISASLIGNMENT WITH SMALL ISLAND)
# If all values of population are missing replace with one (this should never happen, but is from an older code version and retained for insurance/neurosis)
if (sum(!is.na(master$pop))==0) {
  master[, pop := 1]
  
  # If all values of population sum to < 1, replace population with 1 to avoid floating point precision issues
} else if (sum(master$pop, na.rm=TRUE) < 1) {
  master[, pop := 1]
  
  # if there are missing population pixels fill them through interpolation  
} else if (sum(is.na(master$pop)) > 0) {
  master[, pop := na.approx(pop, na.rm=FALSE), by="layer"]
} 
message("Missing and zero populations resolved")

# Drop any pixels with zero or missing populations #
master <- master[pop > 0 & is.na(pop)==F, ]

if (debug==T) masterBkup <- copy(master)




# ### CONVERT TEMP VARIABLES AND CREATE TEMPERATURE DRAWS ###
master[, dailyTempCat := 10 * (dailyTempCat - 273.1)]
master[, sd := 10 * sd]
master <- master[, .SD, .SDcols = c("meanTempCat", "dailyTempCat", "sd", "pop")]
setkeyv(master, c("meanTempCat", "dailyTempCat"))
  
wideStart <- Sys.time()

  
### CREATE DRAWS AND COLLAPSE POPULATION BY DAILY TEMP AND TEMP ZONE ###
  collapse.draws <- function(draw) { 
    longTemp <- copy(master)[, dailyTempCat := as.integer(round(sd * Rnorm(.N, 0, 1) + dailyTempCat))]  
    
    longTemp <- longTemp[, lapply(.SD, sum), by = .(meanTempCat, dailyTempCat), .SDcols = "pop"] 
    longTemp[, draw := as.integer(draw)]
    
    message("Reshape to long & collapse completed for draw ", draw)
    
    return(longTemp)
    }

  long <- do.call(rbind, mclapply(0:max_draw, collapse.draws))
  
  wideEnd <- Sys.time()
  difftime(wideEnd, wideStart)


message("Reshape to long & collapse completed for all zones")

# Save the output #
write.csv(long, file = paste0(outDir, "/", job.name, ".csv"), row.names = F)

end <- Sys.time()
(locProcessTime <- difftime(end, start))
message("File exported -- done")

shpReadTime
dataReadTime
locProcessTime




