############################################################################ ###
#                                                                              #
# Elliot Shayle | Philipps Universit?t Marburg | 26/09/2023                    #
#                                                                              #
# Conduct statistical analysis on Kannenberg's field sites using PlanetScope   #
#                                                                              #
############################################################################ ###

### Set-up the working environment ####

## Step 0 is packages et al, OFC

# Packages

install.packages("sf")
install.packages("sp")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("raster")
install.packages("tidyr")
install.packages("lubridate")
install.packages("TSstudio")
install.packages("plotly")
install.packages("readxl")
install.packages("betareg")
install.packages("blockCV")
install.packages("car")
install.packages("spatialRF")
install.packages("fields")
install.packages("terra")
install.packages("exactextractr")
install.packages("ranger")
install.packages("scales")
install.packages("ggspatial")
install.packages("ggmap")
install.packages("patchwork")

library(sf)
library(sp)
library(ggplot2)
library(dplyr)
library(raster)
library(tidyr)
library(lubridate)
library(TSstudio)
library(plotly)
library(readxl)
library(stringr)
library(betareg)
library(blockCV)
library(car)
library(spatialRF)
library(fields)
library(terra)
library(exactextractr)
library(ranger)
library(scales)
library(ggspatial)
library(ggmap)
library(patchwork)

# Working directory

setwd("C:/Users/ellio/OneDrive/Documents/Work work/Philipps Universität Marburg/Defoliation Timing and Tree Mortality/Code") # On my MS Surface
setwd("C:/Users/Elliot Samuel Shayle/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Code") # On my home PC
setwd("~/Defoliation Timing and Tree Mortality/Code") # On RStudio server IDE

# Optionally load .RData file as some of this code is computationally intensive

load("~/Work work/Philipps Universität Marburg/Defoliation Timing and Tree Mortality/Code/.RData") # Should work on either PC

### Import and aggregate my data because it's split into 32 different folders ####
## Step 1 is to load the manifest file from each folder

# Load the manifest files from all my directories

manifest_data_dir <- list() # create a list to store all the manifest data

Path.to.KBerg.Data <- "C:/Users/ellio/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps" # This line is just here to make my code look tidy (MS Surface)
Path.to.KBerg.Data <- "C:/Users/Elliot Samuel Shayle/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps" # (UK Desktop)

# This for loop extracts the manifest JSON from each folder
for (i in 1:32) # This will break if I change the number of PlanetScope data folders
  {# The line below won't work correctly if I change the directory order or the name of the manifest files
  manifest_data_dir <- c(manifest_data_dir, paste0(Path.to.KBerg.Data,"/Tree_defoliation_study-Kannenberg_order_", i, "/manifest.json"))}

combined_manifest_data <- list() # Create a list to combine manifest data within

for (i in 1:length(manifest_data_dir)) {
  # Read and parse the JSON data from each manifest file
  manifest_data <- jsonlite::fromJSON(manifest_data_dir[[i]], flatten = TRUE)
  
  # Add the data to the combined list
  combined_manifest_data <- c(combined_manifest_data, list(manifest_data))
} 
### This won't work if my file paths don't have the umlauts correctly encoded

# Next, put the directory location of each file in each dataframe

for (i in 1:length(combined_manifest_data)) {
tmp.file.dirs <- as.data.frame(rep( # Make a lone DF column to attach to manifest DF
    paste0(Path.to.KBerg.Data,"/Tree_defoliation_study-Kannenberg_order_", i, "/"), # If I ever change the file locations, it will break this code
                     length(combined_manifest_data[[i]][["files"]]$path))) # Make the column the length of the manifest DF
names(tmp.file.dirs)[1] <- "Root to dir" # Rename the DF's column to something human readable
combined_manifest_data[[i]][["files"]] <- bind_cols(tmp.file.dirs, combined_manifest_data[[i]][["files"]]) # Join my column onto the manifest DF
rm(tmp.file.dirs) # Clean up when you're done :)
}

## Now I want to combine all the manifests into one easily searchable DF object

# Remove the extraneous "names" variable from each list

for (i in 1:length(combined_manifest_data)) {
  combined_manifest_data[[i]] <- combined_manifest_data[[i]][["files"]]  
}

# Append each dataframe together into one super dataframe

combined_manifest_DF <- do.call(rbind, combined_manifest_data)

### Start processing my metadata ####
## Get a list of dates which I have data for (this will be my X axis)

## Option 1: Use string operators to extract dates from the file names

# Get dates

subset_manifest_TIFFs <- combined_manifest_DF[combined_manifest_DF$media_type == "image/tiff", ]

subset_manifest_TIFFs$date <- substr(subset_manifest_TIFFs$path, 1, 10)

subset_manifest_TIFFs$date <- as.Date(subset_manifest_TIFFs$date, format = "%Y-%m-%d")

## Create new columns for each site to record presence or absence within the GeoTIFF

# List of site names
site_names <- c("AR1", "AR2", "AR3", "AR4", "CM1", "CM2", "CM3", "CM4", "MD1", "MD2", "MD3", "MD4")

# Create new columns with site names followed by "present" and fill with NA values
for (site in site_names) {
  col_name <- paste(site, "present", sep = " ")
  subset_manifest_TIFFs[, col_name] <- NA
}

## Option 2: Import KBerg_eAOI_Results_DF.csv and match strip IDs to dates

KBerg_eAoI_Search_Results <- read.csv("KBerg_eAOI_Results_DF.csv")

### Intercept Kannenberg sites with PlanetScope Strips to determine site coverage ####

## Import and prep Kannenberg's field sites GeoJSON

KBerg_eAoI <- st_read("C:/Users/ellio/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Analyses/Kannenberg Margins over Fieldsites NAD 83 CRS as GeoJSON.geojson") # (MS Surface)
"KBerg_co-ords" <-st_read("C:/Users/ellio/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Analyses/Kannenberg Sites Unbuffered NAD 83 CRS.geojson") # (MS Surface)

KBerg_AoI <- st_read("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Analyses/Kannenberg Buffered Sites NAD 83 CRS as Shape File.shp")

KBerg_eAoI <- st_read("C:/Users/Elliot Samuel Shayle/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Analyses/Kannenberg Margins over Fieldsites NAD 83 CRS as GeoJSON.geojson") # (UK Desktop)
"KBerg_co-ords" <-st_read("C:/Users/Elliot Samuel Shayle/OneDrive/Documents/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Analyses/Kannenberg Sites Unbuffered NAD 83 CRS.geojson") # (UK Desktop)

## Split off each site
# These are the extended AoI sites (with margins around the fieldsite)

site_eAR1 <- KBerg_eAoI %>% filter(Site == "AR1")
site_eAR2 <- KBerg_eAoI %>% filter(Site == "AR2")
site_eAR3 <- KBerg_eAoI %>% filter(Site == "AR3")
site_eAR4 <- KBerg_eAoI %>% filter(Site == "AR4")
site_eCM1 <- KBerg_eAoI %>% filter(Site == "CM1")
site_eCM2 <- KBerg_eAoI %>% filter(Site == "CM2")
site_eCM3 <- KBerg_eAoI %>% filter(Site == "CM3")
site_eCM4 <- KBerg_eAoI %>% filter(Site == "CM4")
site_eMD1 <- KBerg_eAoI %>% filter(Site == "MD1")
site_eMD2 <- KBerg_eAoI %>% filter(Site == "MD2")
site_eMD3 <- KBerg_eAoI %>% filter(Site == "MD3")
site_eMD4 <- KBerg_eAoI %>% filter(Site == "MD4")

# These are just Kannenberg's fieldsites, no extra margins, circular

site_AR1 <- KBerg_AoI %>% filter(Site == "AR1")
site_AR2 <- KBerg_AoI %>% filter(Site == "AR2")
site_AR3 <- KBerg_AoI %>% filter(Site == "AR3")
site_AR4 <- KBerg_AoI %>% filter(Site == "AR4")
site_CM1 <- KBerg_AoI %>% filter(Site == "CM1")
site_CM2 <- KBerg_AoI %>% filter(Site == "CM2")
site_CM3 <- KBerg_AoI %>% filter(Site == "CM3")
site_CM4 <- KBerg_AoI %>% filter(Site == "CM4")
site_MD1 <- KBerg_AoI %>% filter(Site == "MD1")
site_MD2 <- KBerg_AoI %>% filter(Site == "MD2")
site_MD3 <- KBerg_AoI %>% filter(Site == "MD3")
site_MD4 <- KBerg_AoI %>% filter(Site == "MD4")

## intercept with PlanetScope data to determine availability

# AR1: Do the overlap / intersect / intercept check!

subset_manifest_TIFFs$`AR1 present` <- FALSE  # Initialize the column with FALSE values

#AR1_present_vector <- subset_manifest_TIFFs$`AR1 present`

# Number of GeoTIFF files
num_files <- nrow(subset_manifest_TIFFs)

# Create my for loop to determine presence or absence

for (i in 1:num_files) {
  # Start by loading in the GeoTIFF as a Raster object
  tmp.GeoTIFF <- raster(file.path(subset_manifest_TIFFs$`Root to dir`[i], subset_manifest_TIFFs$path[i]))
  tmp.GeoTIFF <- reclassify(tmp.GeoTIFF, cbind(0, NA), include.lowest = T)
  
  # Create a nested for loop to check each site for each row!
  for (j in 1:dim(`KBerg_co-ords`)[1]) {
    
    # Check whether the coordinate exists within the Raster
    tmp.presence <- raster::extract(tmp.GeoTIFF, `KBerg_co-ords`[j,]) # Checks if the site in the integer is within the raster
  
  # If the site is within the raster, then return true into the table
  if (is.na(tmp.presence)) {
    subset_manifest_TIFFs[i, c(7+j)] <- F} # If it's true that tmp.presence contains NA, then put a false (absent) value in the cell
  else {subset_manifest_TIFFs[i, c(7+j)] <- T}
    print(c(i, j))
  }
  
  # Tidy up afterwards :)
  rm(tmp.GeoTIFF) 
  rm(tmp.presence)
  print("Tidied up :)")
}

# Save output, because this code above takes hours to run!

write.csv(subset_manifest_TIFFs, file = "KBerg_manifest_TIFFs_DF.csv")

### Tidy up "subset_manifest_TIFFs" to be clearer

## First, subset it down to just the GeoTIFFs (no UDMs)

subset_manifest_TIFFs <- subset_manifest_TIFFs %>% filter(!grepl("_udm2", path)) # This subsets the dataframe

# Now I want to tidy it up a bit

names(subset_manifest_TIFFs)[names(subset_manifest_TIFFs) == "path"] <- "GeoTIFF_FileName" # The way this line works is it gets a vector of column names, and then subsets by a vector of Booleans.
  
# Then give each row it's unique strip ID

pattern <- "(?<=^.{17}).{6,7}(?=_composite\\.tif$)" # Create a pattern to recognise the strip ID within the filename

single.StripID.column <- data.frame(Strip_ID = c(regmatches(subset_manifest_TIFFs$GeoTIFF_FileName, regexpr(pattern, subset_manifest_TIFFs$GeoTIFF_FileName, perl = TRUE)))) # Extract the strip ID to a new column using regular expressions functions!

subset_manifest_TIFFs <- bind_cols(single.StripID.column, subset_manifest_TIFFs) # Bind my column of strip IDs to my pre-existing table

## Secondly, add back in the UDMs!

# I'll rename the media type column because it is redundant

names(subset_manifest_TIFFs)[names(subset_manifest_TIFFs) == "media_type"] <- "UDM_FileName"

# Now I need to create the UDM filenames. Good thing they all follow a consistent naming scheme!

subset_manifest_TIFFs$UDM_FileName <- sub(".{4}$", "_udm2.tif", tmp.test.df$GeoTIFF_FileName) # Replace the last 4 characters of each GeoTIFF filename with "_udm2.tif"

## I may want to amend the `Root to dir` column to make it match the system I'm running the code on

# I can make it system agnostic by replacing the WD bits with a "~"

subset_manifest_TIFFs$`Root to dir` <- sub(".*\\/Work work\\/", "~/Work work/", subset_manifest_TIFFs$`Root to dir`)

### Create Cloud coverage metadata for my sites which I can use to enrich my visualisations ####

## Step 1: Create an empty DF to record cloud coverage

Cloud_Coverage <- data.frame(
  Strip_ID = subset_manifest_TIFFs$Strip_ID,
  AR1_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR2_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR3_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR4_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)), 
  CM1_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)), 
  CM2_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  CM3_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  CM4_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD1_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD2_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD3_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD4_Cloud = rep(NA, length(subset_manifest_TIFFs$Strip_ID))
)

# Whilst I'm here, I'll reproject KBerg's AoI for analysis with the PlanetScope Usable Data Mask

KBerg_AoI_Reproj <- st_transform(KBerg_AoI, crs(raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[1],subset_manifest_TIFFs$UDM_FileName[1]))))
# The above line takes the CRS of a UDM file and then transforms Kannenberg's co-ordinates to match that CRS

## Step 2: Intersect each KBerg site over each eAoI data mask's cloud coverage layer
# But only for sites which are present in that strip
## Step 3: See what % of pixels are masked and save it to a DF
# Step 2 and 3 are combined into one big nested for loop

for (q in 1:length(KBerg_AoI_Reproj$Site)){
  for (i in 1:nrow(subset_manifest_TIFFs)) { # Going down through each column
    if (subset_manifest_TIFFs[i,8+q] == F) # If the ??? of columns in subset_manifest_TIFFs changes, then this code will break
    {Cloud_Coverage[i,1+q] = NA} # If site not present in strip, then "NA"
    else {c(
      tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$UDM_FileName[i])), # Load the relevant UDM (including 0 values)
      tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[q,]), # Crop it down to the Kannenberg site in this iteration of the for loop
      tmp.site.sum.cloud <- cellStats(tmp.UDM, "sum", na.rm = TRUE), # Get the sum of cloud covered pixels (each covered pixel = 1)
      tmp.cloud.coverage <- (tmp.site.sum.cloud[6] / length(tmp.UDM@data@values[,6])), # Divide the total number of pixels by the number of covered pixels
      Cloud_Coverage[i,1+q] <- tmp.cloud.coverage, # Put the value for the amount of cloud coverage into its cell in the cloud coverage DF (as a decimal number)
      print(c(q, i, tmp.cloud.coverage)))}}
}

# Append the date information to the Cloud_Coverage DF to enable time series
Cloud_Coverage <- bind_cols(subset_manifest_TIFFs$date, Cloud_Coverage) # Simply bind the date to Cloud_Coverage (make sure order is correct or it will misalign)
names(Cloud_Coverage)[names(Cloud_Coverage) == "...1"] <- "Date" # Rename the date column to something clearer

# Check the results
head(Cloud_Coverage)

# Testing bits VVV ##############################################################
# Chat GPT 1st go:

#for (i in 1:nrow(subset_manifest_TIFFs)) { # Going down through each column
#  if (subset_manifest_TIFFs$`AR1 present`[i] == F) {Cloud_Coverage$AR1_Cloud[i] = NA} # If site not present in strip, then "NA"
#  else {
#    tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$UDM_FileName[i])) # Load the relevant UDM (including 0 values)
#    tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[1,]) # Crop it down to just the first Kannenberg site
#    tmp.cloud.mask <- tmp.UDM$cloud > 0 # Turn all pixels with any amount of cloud into a "TRUE" value
#    tmp.cloud.pixels <- sum(tmp.cloud.mask, na.rm = TRUE) # Take the sum of how many cloud pixels there are in the site
#    tmp.pixel.count <- cellStats(tmp.UDM, "sum", na.rm = TRUE) # Count all available pixels (should be 100, but useful for sites split across strips)
#    tmp.cloud.coverage <- (length(tmp.UDM@data@values[,6]) / tmp.pixel.count[6]) * 100 # Divide cloud pixels by total pixels to get a percentage of cloud cover
#    Cloud_Coverage$AR1_Cloud[i] <- tmp.cloud.coverage # Finally, save that value to the cloud coverage table
 #   print(i)
#   }
#}

# All my attempt: (1 site only)

#for (i in 1:nrow(subset_manifest_TIFFs)) { # Going down through each column
#  if (subset_manifest_TIFFs$`AR1 present`[i] == F) {Cloud_Coverage$AR1_Cloud[i] = NA} # If site not present in strip, then "NA"
#  else {c(
#    tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$UDM_FileName[i])), # Load the relevant UDM (including 0 values)
#    tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[1,]), # Crop it down to just the first Kannenberg site
#    tmp.site.sum.cloud <- cellStats(tmp.UDM, "sum", na.rm = TRUE), # Get the sum of cloud covered pixels (each covered pixel = 1)
#    tmp.cloud.coverage <- (tmp.site.sum.cloud[6] / length(tmp.UDM@data@values[,6])), # Divide the total number of pixels by the number of covered pixels
#    Cloud_Coverage$AR1_Cloud[i] <- tmp.cloud.coverage, # Put the value for the amount of cloud coverage into its cell in the cloud coverage DF
#    print(c(i, tmp.cloud.coverage)))}}

# My attempt: (all sites)

#for (q in 1:length(KBerg_AoI_Reproj$Site)){
#  for (i in 1:nrow(subset_manifest_TIFFs)) { # Going down through each column
#    if (subset_manifest_TIFFs[i,8+q] == F) # If the ??? of columns in subset_manifest_TIFFs changes, then this code will break
#      {Cloud_Coverage[i,1+q] = NA} # If site not present in strip, then "NA"
#    else {c(
#      tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$UDM_FileName[i])), # Load the relevant UDM (including 0 values)
#      tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[q,]), # Crop it down to the Kannenberg site in this iteration of the for loop
#      tmp.site.sum.cloud <- cellStats(tmp.UDM, "sum", na.rm = TRUE), # Get the sum of cloud covered pixels (each covered pixel = 1)
#      tmp.cloud.coverage <- (tmp.site.sum.cloud[6] / length(tmp.UDM@data@values[,6])), # Divide the total number of pixels by the number of covered pixels
#      Cloud_Coverage[i,1+q] <- tmp.cloud.coverage, # Put the value for the amount of cloud coverage into its cell in the cloud coverage DF
#      print(c(q, i, tmp.cloud.coverage)))}}
#}

# More testing bits VVV

# Iterate through each row
#for (i in 1:nrow(subset_manifest_TIFFs)) {
#  if (subset_manifest_TIFFs$`AR1 present`[i]) {
#    tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$UDM_FileName[i])) # Load the relevant UDM (including 0 values)
#    tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[1,]) # Crop it down to just the Kannenberg site
#    tmp.cloud.mask <- tmp.UDM$cloud > 0
#    tmp.cloud.pixels <- sum(tmp.cloud.mask, na.rm = TRUE)
    
    # Count all available pixels
#    tmp.pixel.count <- sum(!is.na(tmp.UDM[[1]]))
    
    # Calculate cloud coverage if there are available pixels
#    if (as.logical(tmp.pixel.count > 0)) {  # Convert to logical here
#      tmp.cloud.coverage <- (tmp.cloud.pixels / tmp.pixel.count) * 100
#      Cloud_Coverage$AR1_Cloud[i] <- tmp.cloud.coverage
#    } else {
#      Cloud_Coverage$AR1_Cloud[i] <- NA
#    }
#    print(i)
#  }
#}

#tmp.GeoTIFF <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[1623],subset_manifest_TIFFs$GeoTIFF_FileName[1623]))
#tmp.GeoTIFF <- reclassify(tmp.GeoTIFF, cbind(0, NA), include.lowest = T)
#tmp.GeoTIFF <- raster::crop(tmp.GeoTIFF, KBerg_AoI_Reproj[1,])
#tmp.UDM <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[1623],subset_manifest_TIFFs$UDM_FileName[1623]))
#tmp.UDM <- raster::crop(tmp.UDM, KBerg_AoI_Reproj[1,])
#tmp.cloud.mask <- tmp.UDM$cloud > 0
#plotRGB(tmp.GeoTIFF, r = 3, g = 2, b = 1)

# Testing bits ^^^ ###############################

### Extract NDVI & separate RGB values and put these in a dataframe ####

## Step 1: Create a dataframe to store mean NDVI values for each site

Mean_Site_NDVI <- data.frame(
  Strip_ID = subset_manifest_TIFFs$Strip_ID,
  Date = subset_manifest_TIFFs$date,
  AR1_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR2_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR3_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  AR4_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)), 
  CM1_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)), 
  CM2_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  CM3_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  CM4_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD1_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD2_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD3_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID)),
  MD4_NDVI = rep(NA, length(subset_manifest_TIFFs$Strip_ID))
)

## Step 2, 3, & 4: Using KBerg's study sites as masks, mask my GeoTIFFs and calculate mean NDVI
# This will be with the mask function from the raster package
# The for loop will only run for cells with cloud coverage of 0 from the cloud coverage DF

for (n in 1:length(KBerg_AoI_Reproj$Site)){ # Iterates for each Kannenberg site
  for (i in 1:nrow(Mean_Site_NDVI)) { # Going down through each column
    if (is.na(Cloud_Coverage[i, 2 + n]) || Cloud_Coverage[i, 2 + n] != 0) # Checks for both site presence and cloud coverage. If NOT present or WITH cloud, then gives NA. Otherwise, runs the function. If there is cloud coverage If the ??? of columns in Mean_Site_NDVI changes, then this code will break
    {Mean_Site_NDVI[i,2+n] = NA} # If site is NA or <0, then put "NA" in Mean_Site_NDVI
    else {c(
      tmp.masked.site.GeoTIFF <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i],subset_manifest_TIFFs$GeoTIFF_FileName[i])), # Load the whole GeoTIFF containing Kannenberg's site
      tmp.masked.site.GeoTIFF <- raster::crop(tmp.masked.site.GeoTIFF, KBerg_AoI_Reproj[n,]), # Crop it down to just the relevant site (saves computation)
      tmp.masked.site.GeoTIFF[[5]] <- (tmp.masked.site.GeoTIFF[[4]] - tmp.masked.site.GeoTIFF[[3]]) / (tmp.masked.site.GeoTIFF[[4]] + tmp.masked.site.GeoTIFF[[3]]), # This calculates the site's NDVI and writes it to the (empty) 5th layer of my GeoTIFF
      tmp.masked.site.GeoTIFF <- raster(tmp.masked.site.GeoTIFF, 5), # Load just layer 5 (the NDVI values)
      tmp.masked.site.GeoTIFF <- raster::mask(tmp.masked.site.GeoTIFF, KBerg_AoI_Reproj[n,]),
      tmp.mean.NDVI <- mean(tmp.masked.site.GeoTIFF@data@values, na.rm = TRUE),
      Mean_Site_NDVI[i,2+n] <- tmp.mean.NDVI, # Put the site's mean NDVI value into its cell in the sites' NDVI DF (as a decimal number)
      print(c(n, i, tmp.mean.NDVI)))}}
}

## Step 5: Extract separate RGB values for each site
# This can basically repeat the previous code but without calculating NDVI

## Extracting mean Red, Green, and Blue values for each site

tst.Mean_Site_NDVI <- Mean_Site_NDVI[1:14]
tmp.loop.counter <- 0

for (n in 1:length(KBerg_AoI_Reproj$Site)){ # Iterates for each Kannenberg site
  for (i in 1:nrow(Mean_Site_NDVI)) { # Iterates through each row
    if (is.na(Cloud_Coverage[i, 2 + n]) || Cloud_Coverage[i, 2 + n] != 0) { 
      # If site is NA or cloud coverage > 0, assign NA for RGB values
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Red")] <- NA  
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Green")] <- NA  
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Blue")] <- NA 
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_NIR")] <- NA
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_EVI")] <- NA
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_SAVI")] <- NA
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_GNDVI")] <- NA
      tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_MSAVI2")] <- NA
      tmp.loop.counter <- tmp.loop.counter + 1
      print(tmp.loop.counter)
    } else {
      c(
        tmp.masked.site.GeoTIFF <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[i], subset_manifest_TIFFs$GeoTIFF_FileName[i])), # Load the GeoTIFF
        tmp.masked.site.GeoTIFF <- raster::crop(tmp.masked.site.GeoTIFF, KBerg_AoI_Reproj[n,]), # Crop to site
        
        # Extract and mask Red (Layer 3), Green (Layer 2), Blue (Layer 1), and Near Infra-Red (Layer 4) 
        tmp.masked.Red <- raster::mask(raster(tmp.masked.site.GeoTIFF, 3), KBerg_AoI_Reproj[n,]),
        tmp.masked.Green <- raster::mask(raster(tmp.masked.site.GeoTIFF, 2), KBerg_AoI_Reproj[n,]),
        tmp.masked.Blue <- raster::mask(raster(tmp.masked.site.GeoTIFF, 1), KBerg_AoI_Reproj[n,]),
        tmp.masked.NIR <- raster::mask(raster(tmp.masked.site.GeoTIFF, 4), KBerg_AoI_Reproj[n,]),
        
        # Compute mean values for each band
        tmp.mean.Red <- mean(tmp.masked.Red@data@values, na.rm = TRUE),
        tmp.mean.Green <- mean(tmp.masked.Green@data@values, na.rm = TRUE),
        tmp.mean.Blue <- mean(tmp.masked.Blue@data@values, na.rm = TRUE),
        tmp.mean.NIR <- mean(tmp.masked.NIR@data@values, na.rm = TRUE),
        
        ## I'm also going to calculate different indices to test their performance
        # Starting with EVI
        tmp.mean.EVI <- mean(2.5 * (tmp.mean.NIR - tmp.mean.Red) / (tmp.mean.NIR + 6 * tmp.mean.Red - 7.5 * tmp.mean.Blue + 1)),
        
        # Soil Adjusted Vegetation Index (SAVI) - This uses an adjustable value to compensate for soil. I will try with 0.5
        tmp.mean.SAVI <- mean(((tmp.mean.NIR - tmp.mean.Red) * (1 + 0.5)) / (tmp.mean.NIR + tmp.mean.Red + 0.5)),
        
        # Green Normalised Difference Vegetation Index (GNDVI)
        tmp.mean.GNDVI <- mean((tmp.mean.NIR - tmp.mean.Green) / (tmp.mean.NIR + tmp.mean.Green)),
        
        # Multiple Soil Adjusted Vegetation Index 2 (MSAVI2)
        tmp.mean.MSAVI2 <- mean((2 * tmp.mean.NIR + 1 - sqrt((2 * tmp.mean.NIR + 1)^2 - 8 * (tmp.mean.NIR - tmp.mean.Red))) / 2),
        
        # Save values into Mean_Site_NDVI using dynamic column names
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Red")] <- tmp.mean.Red,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Green")] <- tmp.mean.Green,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_Blue")] <- tmp.mean.Blue,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_NIR")] <- tmp.mean.NIR,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_EVI")] <- tmp.mean.EVI,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_SAVI")] <- tmp.mean.SAVI,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_GNDVI")] <- tmp.mean.GNDVI,
        tst.Mean_Site_NDVI[i, paste0(sub("_NDVI", "", names(tst.Mean_Site_NDVI)[2 + n]), "_MSAVI2")] <- tmp.mean.MSAVI2,
        
        tmp.loop.counter <- tmp.loop.counter + 1,
        
        print(tst.Mean_Site_NDVI[i,]),
        print(tmp.loop.counter)
      )
    }
  }
}

View(tst.Mean_Site_NDVI) # Check that it all calculated correctly
Mean_Site_NDVI <- tst.Mean_Site_NDVI # If it did calculate correctly (it did), then overwrite the previous Mean_Site_NDVI object

# ChatGPT testing section VVV ###############################################################

# Define a function to calculate mean NDVI for a given site and row
calculate_mean_NDVI <- function(site_index, row_index) {
  if (is.na(Cloud_Coverage[row_index, 2 + site_index]) || Cloud_Coverage[row_index, 2 + site_index] != 0) {
    return(NA)
  }
  
  tmp.masked.site.GeoTIFF <- raster::brick(paste0(subset_manifest_TIFFs$`Root to dir`[row_index], subset_manifest_TIFFs$GeoTIFF_FileName[row_index]))
  tmp.masked.site.GeoTIFF <- raster::crop(tmp.masked.site.GeoTIFF, KBerg_AoI_Reproj[site_index,])
  tmp.masked.site.GeoTIFF[[5]] <- (tmp.masked.site.GeoTIFF[[4]] - tmp.masked.site.GeoTIFF[[3]]) / (tmp.masked.site.GeoTIFF[[4]] + tmp.masked.site.GeoTIFF[[3]])
  tmp.masked.site.GeoTIFF <- raster(tmp.masked.site.GeoTIFF, 5)
  tmp.masked.site.GeoTIFF <- raster::mask(tmp.masked.site.GeoTIFF, KBerg_AoI_Reproj[site_index,])
  
  tmp.mean.NDVI <- mean(tmp.masked.site.GeoTIFF@data@values, na.rm = TRUE)
  return(tmp.mean.NDVI)
}

# Apply the function across all combinations of sites and rows
result_matrix <- t(sapply(seq_along(KBerg_AoI_Reproj$Site), function(site_index) {
  sapply(seq_len(nrow(tst.Mean_Site_NDVI)), function(row_index) {
    calculate_mean_NDVI(site_index, row_index)
  })
}))

# Convert the matrix to a data frame
tst.Mean_Site_NDVI[, (2 + seq_along(KBerg_AoI_Reproj$Site))] <- result_matrix

# End of testing section  ^^^ #####################################################################

### Import mortality data + generate cool statistics! ####

## Step 1: Import the .CSV which Kannenberg gave me
# I'm importing 2 tables from my Excel document: The raw data, and the stand composition.

KBerg_Mortality_Data <- read_excel("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/Kannenberg 2021 Utah Visualisations.xlsx", # File path & name
                                   sheet = 1, range = "A1:F1189", progress = T)

# As the dates are written in American (YYYY-DD-MM), the code below reformats it
# This takes the strings from the columns, and reorders it in the correct format (YYYY-MM-DD)
KBerg_Mortality_Data$`Survey date` <- sub("([0-9]{4})-([0-9]{2})-([0-9]{2})", "\\1-\\3-\\2", KBerg_Mortality_Data$`Survey date`)

# This tells R that this is a date format column, not merely a character column
KBerg_Mortality_Data$`Survey date` <- as.Date(KBerg_Mortality_Data$`Survey date`)

# Calculate the mean of canopy dieback for each individual tree
KBerg_Mortality_Data$`Mortality Mean` <- rowMeans(KBerg_Mortality_Data[, c("Mortality 1", "Mortality 2")], na.rm = TRUE)

#KBerg_Stand_Composition <- read_excel("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/Kannenberg 2021 Utah Visualisations.xlsx", # File path & name
#                                      sheet = 2, progress = T) # I trimmed the stand composition data in Excel so it would read into R easily

## Step 2: I now need to generate my stand composition statistics
# Because each site was completely sampled twice, I'm separating my data into the 2 sampling instances

KBerg_Mortality_Data_May <- KBerg_Mortality_Data %>%
  filter(format(`Survey date`, "%m") == "05")  # Selects the May sampling round

KBerg_Mortality_Data_Oct <- KBerg_Mortality_Data %>%
  filter(format(`Survey date`, "%m") == "10")  # Selects the October sampling round

# This custom function can be used to calculate important stand composition statistics for each sampling round
calculate_species_statistics <- function(df) { # The input object should be a dataframe
  df %>% 
    group_by(Site) %>% # The dataframe will then be grouped by site
    summarise( # The data will then be summarised as follows:
      Junipers_in_Site = sum(Species == "Juniper"), # Count all Juniperus osteosperma, put it in the relevant column
      Juniper_Site_Percent = (Junipers_in_Site / n())*100, # Calculate the percentage of the stand which is Juniper
      Juniper_Mean_Dieback_Percent = mean(`Mortality Mean`[Species == "Juniper"], na.rm = TRUE), # Calculate the mean canopy dieback of Juniper using the mean value of each fieldworker's guess
      Pinons_in_Site = sum(Species == "Pinon"), # Repeat for Pinon
      Pinon_Site_Percent = (Pinons_in_Site / n())*100, 
      Pinon_Mean_Dieback_Percent = mean(`Mortality Mean`[Species == "Pinon"], na.rm = TRUE), # So this takes the mean mortality value from just the site in question for just rows where the species is pinon
      Total_Tree_Count = n(), # "n" just means the count
      Total_Mean_Dieback_Percent = mean(`Mortality Mean`, na.rm = TRUE)
    )
}

# Apply the function to May and October's data
KBerg_Stand_Composition_May <- calculate_species_statistics(KBerg_Mortality_Data_May)
KBerg_Stand_Composition_Oct <- calculate_species_statistics(KBerg_Mortality_Data_Oct)

# Export these dataframes as .CSVs so that they can go in my manuscript!
write.csv(KBerg_Stand_Composition_May, "Kannenberg 2021 - Stand Composition in May.csv", row.names = FALSE)
write.csv(KBerg_Stand_Composition_Oct, "Kannenberg 2021 - Stand Composition in October.csv", row.names = FALSE)

## Deprecated section VVV ##########################################################
## Step 3: Generate relevant mortality statistics

# Calculate mean tree mortality for each site and by species
# Starting with mean Juniper and Pinon mortality for each site
#mortality_means <- KBerg_Mortality_Data %>%
#  group_by(Site, Species) %>%
#  summarize(Mean_Mortality = mean(`Mortality Mean`, na.rm = TRUE))

# Put the calculated means into KBerg_Stand_Composition
#KBerg_Stand_Composition <- KBerg_Stand_Composition %>%
#  left_join(mortality_means %>% filter(Species == "Juniper") %>%
#              select(Site, Mean_Mortality) %>%
#              rename(Juniper_Mean_Mortality = Mean_Mortality), by = "Site") %>%
#  left_join(mortality_means %>% filter(Species == "Pinon") %>%
#              select(Site, Mean_Mortality) %>%
#              rename(Pinon_Mean_Mortality = Mean_Mortality), by = "Site")

#KBerg_Stand_Composition <- KBerg_Stand_Composition %>% # Ordering the columns nicely
#  select(Site, 
#         "Junipers in Site", "Juniper as Site %", Juniper_Mean_Mortality,
#         "Pinons in Site", "Pinon as Site %", Pinon_Mean_Mortality,
#         "Total Trees in Site")

# Calculate mean tree mortality for each site, but not grouped by species
#mean_mortality_by_site <- KBerg_Mortality_Data %>%
#  group_by(Site) %>%
#  summarize(Mean_Mortality = mean(`Mortality Mean`, na.rm = TRUE))

# Merge the mean mortality data into the KBerg_Stand_Composition dataframe
#KBerg_Stand_Composition <- left_join(KBerg_Stand_Composition, mean_mortality_by_site, by = "Site")
## End of deprecated section ^^^ #########################################################

## Step 3: Run T-tests to explore differences in mortality

# 2 sample T-test canopy dieback comparison of Pinon vs Juniper individual dieback using May and October samplings together
t_test_entire_dataset <- t.test(KBerg_Mortality_Data$`Mortality Mean`[KBerg_Mortality_Data$Species == "Pinon"], # Just the Pinons
                                KBerg_Mortality_Data$`Mortality Mean`[KBerg_Mortality_Data$Species == "Juniper"], # Just the Junipers
                                paired = F) # Performs a 2 sample T-test because each species is one sample
print("2 Sample T-test for the Entire Dataset: [DON'T USE; each tree sampled twice, once in May, once in Oct]")
print(t_test_entire_dataset)
# Significant

# 2 sample T-test canopy dieback comparison of Pinon vs Juniper individual dieback using May sampling only
T.test_All.Trees_May <- t.test(KBerg_Mortality_Data_May$`Mortality Mean`[KBerg_Mortality_Data_May$Species == "Pinon"], # Just the Pinons
                                KBerg_Mortality_Data_May$`Mortality Mean`[KBerg_Mortality_Data_May$Species == "Juniper"], # Just the Junipers
                                paired = F) 
print("2 Sample T-test Comparing Canopy Dieback Between Pinon and Juniper in May:")
print(T.test_All.Trees_May)
# Significant

# 2 sample T-test canopy dieback comparison of Pinon vs Juniper individual dieback using October sampling only
T.test_All.Trees_Oct <- t.test(KBerg_Mortality_Data_Oct$`Mortality Mean`[KBerg_Mortality_Data_Oct$Species == "Pinon"], # Just the Pinons
                               KBerg_Mortality_Data_Oct$`Mortality Mean`[KBerg_Mortality_Data_Oct$Species == "Juniper"], # Just the Junipers
                               paired = F) 
print("2 Sample T-test Comparing Canopy Dieback Between Pinon and Juniper in October:")
print(T.test_All.Trees_Oct)
# Significant

# Paired sample T-tests to see if there is a significant difference between May and October samplings
# Are there significant differences in canopy dieback for all trees between May and October?
# Note that I'm using the average dieback per site  as the exact number of individuals varies between May and Oct
T.test_All.Trees_May.Oct <- t.test(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent, # All species together, May data only
                                   KBerg_Stand_Composition_Oct$Total_Mean_Dieback_Percent, # All species together, October data only
                                   paired = T) # Paired sample T-test because we're comparing the same samples before and after treatment (the addition of time)
print("Paired Sample T-test Comparing Canopy Dieback for All Species Between May and October:")
print(T.test_All.Trees_May.Oct)
# Not significant (but nearly)

# Are there significant differences in canopy dieback for just Juniper between May and October?
# Note that I'm using the average dieback per site  as the exact number of individuals varies between May and Oct
T.test_Juniper_May.Oct <- t.test(KBerg_Stand_Composition_May$Juniper_Mean_Dieback_Percent, # Juniper only, May data only
                                   KBerg_Stand_Composition_Oct$Juniper_Mean_Dieback_Percent, # Juniper only, October data only
                                   paired = T)
print("Paired Sample T-test Comparing May and October Canopy Dieback for Juniperus osteosperma:")
print(T.test_Juniper_May.Oct)
# Not significant (but nearly)

# Are there significant differences in canopy dieback for just Juniper between May and October?
# Note that as sites AR4, CM4, & MD3 have the same count of Junipers, I can do a paired T test between individuals of these sites only
T.test_Juniper.Individuals_May.Oct <- t.test(KBerg_Mortality_Data_May$`Mortality Mean`[KBerg_Mortality_Data_May$Species == "Juniper" & KBerg_Mortality_Data_May$Site %in% c("AR4", "CM4", "MD3")], # Individual Juniper only, May data only, Sites AR4, CM4, & MD3
                                             KBerg_Mortality_Data_Oct$`Mortality Mean`[KBerg_Mortality_Data_Oct$Species == "Juniper" & KBerg_Mortality_Data_Oct$Site %in% c("AR4", "CM4", "MD3")], # Individual Juniper only, October data only, Sites AR4, CM4, & MD3
                                 paired = T)
print("Paired Sample T-test Comparing May and October Canopy Dieback for Individual Juniperus osteosperma:")
print(T.test_Juniper.Individuals_May.Oct)
# Not significant

# Are there significant differences in canopy dieback for just Pinon between May and October?
# Note that I'm using the average dieback per site  as the exact number of individuals varies between May and Oct
T.test_Pinon_May.Oct <- t.test(KBerg_Stand_Composition_May$Pinon_Mean_Dieback_Percent,
                                 KBerg_Stand_Composition_Oct$Pinon_Mean_Dieback_Percent,
                                 paired = T)
print("Paired Sample T-test Comparing May and October Canopy Dieback for Pinus edulis:")
print(T.test_Pinon_May.Oct)
# Not significant

# Step 4: ANOVA testing to explore dieback differences across sites
# Was Juniper's dieback extent significantly different across sites in May?
ANOVA_Juniper_May <- aov(`Mortality Mean` ~ Site, 
                         data = KBerg_Mortality_Data_May[KBerg_Mortality_Data_May$Species == "Juniper",])
print(c("1-Way ANOVA Test Comparing Juniperus osteosperma's Canopy Dieback Across the 12 Sites in May",
        summary(ANOVA_Juniper_May)))
# Significant
# Tukey's Honestly-Significant Differences test to explore the differences further
TukeyHSD(ANOVA_Juniper_May)

# Was Juniper's dieback extent significantly different across sites in October?
ANOVA_Juniper_Oct <- aov(`Mortality Mean` ~ Site, 
                         data = KBerg_Mortality_Data_Oct[KBerg_Mortality_Data_Oct$Species == "Juniper",])
print(c("1-Way ANOVA Test Comparing Juniperus osteosperma's Canopy Dieback Across the 12 Sites in October",
        summary(ANOVA_Juniper_Oct)))
# Significant
# Tukey's Honestly-Significant Differences test to explore the differences further
TukeyHSD(ANOVA_Juniper_Oct)

# Was Pinon's dieback extent significantly different across sites in May?
ANOVA_Pinon_May <- aov(`Mortality Mean` ~ Site, 
                         data = KBerg_Mortality_Data_May[KBerg_Mortality_Data_May$Species == "Pinon",])
print(c("1-Way ANOVA Test Comparing Pinus edulis's Canopy Dieback Across the 12 Sites in May",
        summary(ANOVA_Pinon_May)))
# Significant 
# Tukey's Honestly-Significant Differences test to explore the differences further
TukeyHSD(ANOVA_Pinon_May)

# Was Pinon's dieback extent significantly different across sites in October?
ANOVA_Pinon_Oct <- aov(`Mortality Mean` ~ Site, 
                         data = KBerg_Mortality_Data_Oct[KBerg_Mortality_Data_Oct$Species == "Pinon",])
print(c("1-Way ANOVA Test Comparing Pinus edulis's Canopy Dieback Across the 12 Sites in October",
        summary(ANOVA_Pinon_Oct)))
# Significant 
# Tukey's Honestly-Significant Differences test to explore the differences further
TukeyHSD(ANOVA_Pinon_Oct)

### Because I know that Pinon and Juniper dieback does significantly differ from one another, I won't do any ANOVAs with them combined

### Q1: Linear regression NDVI and Ground Truthed Mortality ####

## Using all trees in site for May only

# Initialize a vector to store NDVI values
tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))

# Loop through each row of KBerg_Stand_Composition_May
for (i in 1:nrow(KBerg_Stand_Composition_May)) {
  # Get the site name and survey date
  tmp.site <- KBerg_Stand_Composition_May$Site[i]
  tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] # First occurrence
  
  # Find the nearest available NDVI date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]]),] # Create a dataframe of dates with data for the relevant site
  
  tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),] # Subsets to just the row with the nearest matching date
  
  # Get the NDVI value for the nearest date
  tmp.site_NDVI[i] <- Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
}

# Create a dataframe with the variables
LM_Data_May_NDVI.SiteMeanDieback <- data.frame(Total_Mean_Dieback = KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent, Mean_NDVI = tmp.site_NDVI)

## I also want to know what dates and strip ID were actually selected in the end 
## Although I don't want to screw up any subsequent code by messing with my for loops and variables
## So I'm duplicating and amending the code above to work towards that purpose

# Initialize a vector to store NDVI metadata
tmp.site_NDVI_meta <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
tmp.site_stripID_meta <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
tmp.site_Image.date_meta <- c(rep(NA,length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent)))

# Loop through each row of KBerg_Stand_Composition_May
for (i in 1:nrow(KBerg_Stand_Composition_May)) {
  # Get the site name and survey date
  tmp.site <- KBerg_Stand_Composition_May$Site[i]
  tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] # First occurrence
  
  # Find the nearest available NDVI date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]]),] # Create a dataframe of dates with data for the relevant site
  
  tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),] # Subsets to just the row with the nearest matching date
  
  # Get the NDVI value for the nearest date
  tmp.site_NDVI_meta[i] <- Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  
  # Get the Strip ID for the selected image
  tmp.site_stripID_meta[i] <- Mean_Site_NDVI$Strip_ID[Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  
  # Get the date of the Strip ID
  tmp.site_Image.date_meta[i] <- as.character(Mean_Site_NDVI$Date[Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID])
}

# Create a dataframe with the variables
May_Sampling_NDVI_Metadata <- data.frame(row.names = KBerg_Stand_Composition_May$Site, Mean_NDVI = tmp.site_NDVI_meta, Strip_ID = tmp.site_stripID_meta, Imaging_Date = tmp.site_Image.date_meta)

## The site MD4 NDVI value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip ???2319101 from 2019/04/27

LM_Data_May_NDVI.SiteMeanDieback[12,2] <- Mean_Site_NDVI$MD4_NDVI[which(Mean_Site_NDVI$Strip_ID == 2319101)]

# Fit a linear regression model
LM_Model_May_NDVI.SiteMeanDieback <- lm(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_May_NDVI.SiteMeanDieback)

# Summary of the regression model
summary(LM_Model_May_NDVI.SiteMeanDieback) # Summarises model

plot(LM_Model_May_NDVI.SiteMeanDieback) # Plots descriptive figures about model fit

plot(LM_Data_May_NDVI.SiteMeanDieback) # Plots simple X,Y scatter of points

plot(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_May_NDVI.SiteMeanDieback) # Same as before, but with axes swapped

# Plotting the linear regression model with enhanced features
plot(LM_Data_May_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     xlab = "Sitewide Mean of NDVI", ylab = "Sitewide Mean of Observed Relative Canopy Dieback (%)", 
     main = "Linear Regression: Exploring the Relationship Between NDVI and Relative Canopy Dieback in May 2019",
     pch = 16, col = "blue") 

# Adding regression line
abline(LM_Model_May_NDVI.SiteMeanDieback, col = "red")

# Adding site labels
rownames(LM_Data_May_NDVI.SiteMeanDieback) <- paste0(KBerg_Stand_Composition_May$Site, "_May")
text(LM_Data_May_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     labels = rownames(LM_Data_May_NDVI.SiteMeanDieback), pos = 3, col = "black", cex = 0.8)

# Adding formula of the regression line
text(0.32, 20, 
     paste("Regression line: y =", round(coef(LM_Model_May_NDVI.SiteMeanDieback)[1], 2), 
           "+", 
           round(coef(LM_Model_May_NDVI.SiteMeanDieback)[2], 2), "x"), 
     col = "red")

## Using all trees in site for Oct only

# Initialize a vector to store NDVI values
tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_Oct$Total_Mean_Dieback_Percent))

# Loop through each row of KBerg_Stand_Composition_Oct
for (i in 1:nrow(KBerg_Stand_Composition_Oct)) {
  # Get the site name and survey date
  tmp.site <- KBerg_Stand_Composition_Oct$Site[i]
  tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] # First occurrence
  
  # Find the nearest available NDVI date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]]),] # Create a dataframe of dates with data for the relevant site
  
  tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),] # Subsets to just the row with the nearest matching date
  
  # Get the NDVI value for the nearest date
  tmp.site_NDVI[i] <- Mean_Site_NDVI[[paste0(tmp.site, "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
}

# Create a dataframe with the variables
LM_Data_Oct_NDVI.SiteMeanDieback <- data.frame(Total_Mean_Dieback = KBerg_Stand_Composition_Oct$Total_Mean_Dieback_Percent, Mean_NDVI = tmp.site_NDVI)

# Fit a linear regression model
LM_Model_Oct_NDVI.SiteMeanDieback <- lm(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_Oct_NDVI.SiteMeanDieback)

# Summary of the regression model
summary(LM_Model_Oct_NDVI.SiteMeanDieback)

plot(LM_Model_Oct_NDVI.SiteMeanDieback)

plot(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_Oct_NDVI.SiteMeanDieback)

# Plotting the linear regression model with enhanced features
plot(LM_Data_Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     xlab = "Sitewide Mean of NDVI", ylab = "Sitewide Mean of Observed Relative Canopy Dieback (%)", 
     main = "Linear Regression: Exploring the Relationship Between NDVI and Relative Canopy Dieback in October 2019",
     pch = 16, col = "blue") 

# Adding regression line
abline(LM_Model_Oct_NDVI.SiteMeanDieback, col = "red")

# Adding site labels
rownames(LM_Data_Oct_NDVI.SiteMeanDieback) <- paste0(KBerg_Stand_Composition_Oct$Site, "_Oct")
text(LM_Data_Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     labels = rownames(LM_Data_Oct_NDVI.SiteMeanDieback), pos = 3, col = "black", cex = 0.8)

# Adding formula of the regression line
text(0.32, 20, 
     paste("Regression line: y =", round(coef(LM_Model_Oct_NDVI.SiteMeanDieback)[1], 2), 
           "+", 
           round(coef(LM_Model_Oct_NDVI.SiteMeanDieback)[2], 2), "x"), 
     col = "red")

## Using all trees in site for May and October together
# Append the May groundtruth & NDVI dataframe to the October one

LM_Data_May.Oct_NDVI.SiteMeanDieback <- rbind(LM_Data_May_NDVI.SiteMeanDieback, LM_Data_Oct_NDVI.SiteMeanDieback)

# Fit a linear regression model
LM_Model_May.Oct_NDVI.SiteMeanDieback <- lm(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_May.Oct_NDVI.SiteMeanDieback)

summary(LM_Model_May.Oct_NDVI.SiteMeanDieback)

plot(LM_Model_May.Oct_NDVI.SiteMeanDieback)

plot(Total_Mean_Dieback ~ Mean_NDVI, data = LM_Data_May.Oct_NDVI.SiteMeanDieback)

## Plotting the linear regression model with enhanced features
plot(LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     xlab = "Sitewide Mean of NDVI", ylab = "Sitewide Mean of Observed Relative Canopy Dieback (%)", 
     title(main = "Linear Regression of NDVI and the Mean of Sitewide Relative Canopy Dieback"),
     pch = 16, col = "blue") 

# Adding regression line
abline(LM_Model_May.Oct_NDVI.SiteMeanDieback, col = "red")

# Adding site labels
text(LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     labels = rownames(LM_Data_May.Oct_NDVI.SiteMeanDieback), pos = 3, col = "black", cex = 0.8)

# Adding formula of the regression line
text(0.32, 20, 
     paste("Regression line: y =", round(coef(LM_Model_May.Oct_NDVI.SiteMeanDieback)[1], 2), 
           "+", 
           round(coef(LM_Model_May.Oct_NDVI.SiteMeanDieback)[2], 2), "x"), 
     col = "red")

## Plot a professional chart for use in my manuscript (Figure 3)
# Set the plot margins
par(
  mar = c(5, 5, 4, 2),
  cex.lab = 1.4,         # Axis label size
  cex.axis = 1.2,        # Tick label size
  cex.main = 1.6,        # Title size
  lab = c(5, 5, 7)       # Axis tick count: x, y, label
)

# Base scatterplot
plot(
  LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI,
  LM_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback,
  xlab = "Sitewide Mean of NDVI",
  ylab = "Sitewide Mean of Observed Relative Canopy Dieback (%)",
  main = "Linear Relationship Between NDVI and Canopy Dieback",
  pch = 4,
  col = ifelse(ifelse(grepl("May", rownames(LM_Data_May.Oct_NDVI.SiteMeanDieback)), "May", "October") == "May", "blue", "red"),
  cex = 1.75,
  lwd = 2
)

# Add regression line clipped to x-range
regLine(LM_Model_May.Oct_NDVI.SiteMeanDieback,
        col = "black", lwd = 2)

# Add legend (cross symbol)
legend("topright",
       legend = c("May sampling", "October sampling"),
       col = c("blue", "red"),
       pch = 4,
       pt.cex = 2,
       pt.lwd = 2,
       bty = "o")

# Add regression equation text below the legend
text(
  x = 0.32,
  y = 20.5,
  labels = paste0("Regression line: y = ", round(coef(LM_Model_May.Oct_NDVI.SiteMeanDieback)[1], 2),
                  ifelse(round(coef(LM_Model_May.Oct_NDVI.SiteMeanDieback)[2], 2) >= 0, " + ", " − "), 
                  abs(round(coef(LM_Model_May.Oct_NDVI.SiteMeanDieback)[2], 2)), " × x"),
  cex = 1.2,
  adj = 0
)

### Q2.1: Enhanced modelling approaches and testing ####

## Beta-regression of all trees in site for May and October together
# As my data are bounded (percentage dieback) beta regression is a relevant analytical tool to circumvent this constraint
# Further supported by DOI: 10.1002/ecs2.3940

# Creating the dataframe for analysis
BetaReg_Data_May.Oct_NDVI.SiteMeanDieback <- LM_Data_May.Oct_NDVI.SiteMeanDieback

# Taking the canopy dieback data and converting it to a proportion between 0 and 1
BetaReg_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback <- LM_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback/100

# Running the beta regression
BetaReg_Model_May.Oct_NDVI.SiteMeanDieback <- betareg(Total_Mean_Dieback ~ Mean_NDVI, data = BetaReg_Data_May.Oct_NDVI.SiteMeanDieback)

# Now see the results!
summary(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback)

## Plots of the model (adapted from ChatGPT - double check)
plot(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback) # basic overview of the model

# Predicted Vs. Actual dieback plot
plot(BetaReg_Data_May.Oct_NDVI.SiteMeanDieback$Total_Mean_Dieback, 
     predict(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback, type = "response"), 
     xlab = "Observed Canopy Dieback", 
     ylab = "Predicted Canopy Dieback", 
     main = "Observed vs. Predicted Canopy Dieback",
     pch = 16, col = "blue", xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, col = "red")

# Residuals against predicted values
plot(predict(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback, type = "response"), 
     residuals(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback, type = "response"), 
     xlab = "Predicted Canopy Dieback", 
     ylab = "Residuals", 
     main = "Residuals vs. Predicted Values",
     pch = 16, col = "blue")
abline(h = 0, col = "red")

# NDVI's effect on canopy dieback
plot(BetaReg_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, 
     fitted(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback), 
     xlab = "Mean NDVI", 
     ylab = "Predicted Canopy Dieback", 
     main = "Effect of NDVI on Canopy Dieback",
     pch = 16, col = "blue")
lines(lowess(BetaReg_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, 
             fitted(BetaReg_Model_May.Oct_NDVI.SiteMeanDieback)), 
      col = "red")

## Multiple regression of just Juniper dieback within sites (May & Oct sampling together)
# First I need to add variables to the LM_Data_May.Oct_NDVI.SiteMeanDieback dataframe

KBerg_Stand_Composition_May.Oct <- rbind(KBerg_Stand_Composition_May, KBerg_Stand_Composition_Oct) # Combine my stand composition DFs so that I can bind columns onto other DF

KBerg_Stand_Composition_May.Oct$Site <- paste0( # Add collection month to site ID
  KBerg_Stand_Composition_May.Oct$Site, c(rep("_May", nrow(KBerg_Stand_Composition_May)), rep("_Oct", nrow(KBerg_Stand_Composition_Oct))))

LM_Data_May.Oct_NDVI.SiteMeanDieback <- bind_cols(LM_Data_May.Oct_NDVI.SiteMeanDieback, # Append sitewide average of Juniper dieback
                                                  select(KBerg_Stand_Composition_May.Oct, Juniper_Mean_Dieback_Percent))

LM_Data_May.Oct_NDVI.SiteMeanDieback <- bind_cols(LM_Data_May.Oct_NDVI.SiteMeanDieback, # Append proportion of Junipers in site
                                                  select(KBerg_Stand_Composition_May.Oct, Juniper_Site_Percent))

# Now I'm going to conduct the multiple regression

GLM_Model_May.Oct_NDVI.JuniperMeanDieback <- lm(Mean_NDVI ~ Juniper_Mean_Dieback_Percent + Juniper_Site_Percent, 
                                                data = LM_Data_May.Oct_NDVI.SiteMeanDieback)

summary(GLM_Model_May.Oct_NDVI.JuniperMeanDieback) 

# Plotting the linear regression model with enhanced features
plot(LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May.Oct_NDVI.SiteMeanDieback$Juniper_Mean_Dieback_Percent, 
     xlab = "Mean NDVI", ylab = "Juniper Specific Mean Canopy Dieback (%)", 
     main = "Multiple Linear Regression: Mean NDVI vs Juniper Canopy Dieback and Site Proportion with May & Oct Data",
     pch = 16, col = "blue") 

# Adding regression line
abline(GLM_Model_May.Oct_NDVI.JuniperMeanDieback, col = "red") # Doesn't render properly

# Adding site labels
text(LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI, LM_Data_May.Oct_NDVI.SiteMeanDieback$Juniper_Mean_Dieback_Percent, 
     labels = rownames(LM_Data_May.Oct_NDVI.SiteMeanDieback), pos = 3, col = "black", cex = 0.8)

# Adding formula of the regression line
text(0.32, 20, 
     paste("Regression line: y =", round(coef(GLM_Model_May.Oct_NDVI.JuniperMeanDieback)[1], 2), 
           "+", 
           round(coef(GLM_Model_May.Oct_NDVI.JuniperMeanDieback)[2], 2), "x"), 
     col = "red")

### Q2.2: Can I use remote sensing data to spot the difference in dieback between species? ####
# The premise of this analysis is that I want to test the resolving power of my remote sensing dataset
# I know that there is a biologically significant difference between dieback extent for each species
# But does it also appear when I test using Remote Sensing data?

## Step 0: Get all the important data together

# Add a column for site-wide Mean NDVI data
KBerg_Stand_Composition_May.Oct <- bind_cols(KBerg_Stand_Composition_May.Oct, LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI) # Check data is in the same order!
names(KBerg_Stand_Composition_May.Oct)[10] <- names(LM_Data_May.Oct_NDVI.SiteMeanDieback[names(LM_Data_May.Oct_NDVI.SiteMeanDieback) == "Mean_NDVI"]) # Name the newly appended column

## Step 1: Check the suitability of the data for this analysis

# Linearity
plot(Total_Mean_Dieback_Percent ~ Mean_NDVI, data = KBerg_Stand_Composition_May.Oct, col = rainbow(nrow(KBerg_Stand_Composition_May.Oct)/2))
abline(lm(Total_Mean_Dieback_Percent ~ Mean_NDVI, data = KBerg_Stand_Composition_May.Oct))

# Homogeneity of Regression Slopes
SpeciesComp.NDVI_Interaction_Model <- lm(Total_Mean_Dieback_Percent ~ Mean_NDVI * Juniper_Site_Percent, data = KBerg_Stand_Composition_May.Oct)
anova(SpeciesComp.NDVI_Interaction_Model)  # Check if the interaction term is significant

## Step 2: Run the analysis as an ANCOVA test

# Fit the ANCOVA model
SpeciesComp.NDVI_ANCOVA_Model <- lm(Total_Mean_Dieback_Percent ~ Mean_NDVI + Juniper_Site_Percent, 
                   data = KBerg_Stand_Composition_May.Oct)

## Step 3: Run the analysis as a Beta Regression

# Rescale the response variable to be between 0 and 1
KBerg_Stand_Composition_May.Oct$Total_Mean_Dieback_Decimal <- 
  (KBerg_Stand_Composition_May.Oct$Total_Mean_Dieback_Percent / 100) # Only the response variable needs to be rescaled

# Fit the beta regression model
SpeciesComp.NDVI_BetaReg_Model <- betareg(
  Total_Mean_Dieback_Decimal ~ Mean_NDVI * Juniper_Site_Percent, 
  data = KBerg_Stand_Composition_May.Oct)

## Step 3.5: Run the analysis with Dirk's suggested amendment to the ratio variable
# Dirk thought that if I multiply Juniper proportion by NDVI, and then use that variable in the model
# Then it might be more effective as a model, whilst still incorporating proportion data

# Create a variable for the combined ratio and Juniper site proportion
KBerg_Stand_Composition_May.Oct$Combined_Composition.NDVI_decimal <- # New variable name
  KBerg_Stand_Composition_May.Oct$Mean_NDVI * (1 - (KBerg_Stand_Composition_May.Oct$Juniper_Site_Percent/100)) # composition converted to decimal, decimal multiplied by NDVI

# Fit the beta regression model
SpeciesComp.NDVI_BetaReg_Model.v2 <- betareg(
  Total_Mean_Dieback_Decimal ~ KBerg_Stand_Composition_May.Oct$Combined_Composition.NDVI_decimal, 
  data = KBerg_Stand_Composition_May.Oct)

# Fit a "null hypothesis" model without species composition to compare against
MeanDieback.NDVI_BetaReg_Model <- betareg(
  Total_Mean_Dieback_Decimal ~ KBerg_Stand_Composition_May.Oct$Mean_NDVI, # Using the standard, unweighted NDVI values
  data = KBerg_Stand_Composition_May.Oct)

## Step 4: Read the results and determine the models' suitability
## ANCOVA first

# Summarize the ANCOVA results
summary(SpeciesComp.NDVI_ANCOVA_Model)

# Check the significance of the Species Composition
anova(SpeciesComp.NDVI_ANCOVA_Model)

# Check the normality of the residuals with this histogram
hist(residuals(SpeciesComp.NDVI_ANCOVA_Model), main = "Residuals", xlab = "Residuals", col = "skyblue")

# Check the normality of the residuals with this Shapiro-Wilk test
shapiro.test(residuals(SpeciesComp.NDVI_ANCOVA_Model))

# Homoscedasticity with a Non-constant variance test
ncvTest(SpeciesComp.NDVI_ANCOVA_Model)

## Beta regression results (Straightforward model)

# Print the model summary
summary(SpeciesComp.NDVI_BetaReg_Model)

# Calculate confidence intervals for the coefficients
confint(SpeciesComp.NDVI_BetaReg_Model)

# Pseudo-R^2 value
print(c("Pseudo R-squared:", SpeciesComp.NDVI_BetaReg_Model$pseudo.r.squared))

## Beta regression results (Dirk's suggested model)

# Print the model summary
summary(SpeciesComp.NDVI_BetaReg_Model.v2)

# Calculate confidence intervals for the coefficients
confint(SpeciesComp.NDVI_BetaReg_Model.v2)

# Pseudo-R^2 value
print(c("Pseudo R-squared:", SpeciesComp.NDVI_BetaReg_Model.v2$pseudo.r.squared))

## Beta regression results (No NDVI weighting (default) model)

# Print the model summary
summary(MeanDieback.NDVI_BetaReg_Model)

# Calculate confidence intervals for the coefficients
confint(MeanDieback.NDVI_BetaReg_Model)

# Pseudo-R^2 value
print(c("Pseudo R-squared:", MeanDieback.NDVI_BetaReg_Model$pseudo.r.squared))

### Q2.3: Modelling species-specific canopy dieback ####
## This follows the same principles as my straightforward beta models
## The primary difference is that I will create separate Juniper and Pi?on dieback models

# Step 1: Rescale the response variables to be a decimal
# Rescale Juniper to be a decimal
KBerg_Stand_Composition_May.Oct$Juniper_Mean_Dieback_Decimal <- 
  (KBerg_Stand_Composition_May.Oct$Juniper_Mean_Dieback_Percent / 100)

# Rescale Pi?on to be a decimal
KBerg_Stand_Composition_May.Oct$Pinon_Mean_Dieback_Decimal <- 
  (KBerg_Stand_Composition_May.Oct$Pinon_Mean_Dieback_Percent / 100)

# I must add a tiny amount to each of the Pi?on dieback values because the dataset contains 0 values
KBerg_Stand_Composition_May.Oct$Pinon_Mean_Dieback_Decimal <- KBerg_Stand_Composition_May.Oct$Pinon_Mean_Dieback_Decimal + 0.00000000001

# Step 2: Run a beta regression with species composition predicting only Juniper dieback
BetaRegModel_Juniper.SpeciesComp.NDVI_May.Oct <- betareg(
  Juniper_Mean_Dieback_Decimal ~ Mean_NDVI * Juniper_Site_Percent, 
  data = KBerg_Stand_Composition_May.Oct)

# Step 3: Run a beta regression with species composition predicting only Pi?on dieback
BetaRegModel_Pinon.SpeciesComp.NDVI_May.Oct <- betareg(
  KBerg_Stand_Composition_May.Oct$Pinon_Mean_Dieback_Decimal ~ Mean_NDVI * Juniper_Site_Percent, 
  data = KBerg_Stand_Composition_May.Oct)

# Step 4: compare results of these beta regressions

summary(SpeciesComp.NDVI_BetaReg_Model) # Both tree species combined
summary(BetaRegModel_Juniper.SpeciesComp.NDVI_May.Oct) # Juniper dieback results only
summary(BetaRegModel_Pinon.SpeciesComp.NDVI_May.Oct) # Pi?on dieback results only

### Use my data to plot more charts (data availability and descriptive) ####

### Starting with a chart of data availability (for descriptive purposes)

## A simple X+Y for each site with availability on the X axis, and site on Y
# It might be a good idea to make one for each year

# Reshape the data from wide to long format
subset_manifest_long <- subset_manifest_TIFFs %>%
  pivot_longer(cols = ends_with("present"), names_to = "Site", values_to = "Present")

# Remove the word 'present' from strings

subset_manifest_long$Site <- substr(subset_manifest_long$Site, 1, nchar(subset_manifest_long$Site)-8)

# Filter data for a specific year, e.g., 2022
year_to_plot <- 2019

# Create the boxplot for the selected year
ggplot(subset_manifest_long %>% filter(year(date) == year_to_plot), aes(x = date, y = Site, fill = Present)) +
  geom_tile() +
  scale_fill_manual(values = c("TRUE" = "steel blue", "FALSE" = "transparent")) +
  labs(
    x = "Date",
    y = "Site",
    fill = "Presence"
  ) +
  theme_minimal()+
  scale_x_date(
    date_breaks = "3 months",      # Major tickmarks every 3 months
    date_minor_breaks = "1 month" # Minor tickmarks every month
  ) # This creates 1 combined chart

for (i in 1:length(unique(subset_manifest_long$Site))) {
  # Create the boxplot for the selected year
  ggplot(subset_manifest_long %>% filter(year(date) == year_to_plot) %>% filter(Site == unique(subset_manifest_long$Site)[i]), # Use the for loop to select a different unique site each time
         aes(x = date, y = unique(subset_manifest_long$Site)[i], fill = Present)) +
    geom_tile() +
    scale_fill_manual(values = c("TRUE" = "steel blue", "FALSE" = "transparent")) +
    labs(
      x = "Date",
      y = "Site",
      fill = "Presence"
    ) +
    theme_minimal()+
    scale_x_date(
      date_breaks = "3 months",      # Major tickmarks every 3 months
      date_minor_breaks = "1 month" # Minor tickmarks every month
    ) # This creates a separate chart for each year
}

## Add a VLine for the the sampling date

## Add a callout for mean site mortality

### Make charts showing cloud coverage through time

# This function plots each site's time series in the same window
ts_plot(Cloud_Coverage,
        line.mode = "markers",
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = F, # Slider is disabled for the "multiple" plot type
        type = "multiple",
        Xtitle = "Date of Imaging",
        Ytitle = "Cloud Coverage (0 = none, 1 = total coverage)", 
        title = "A Time Series Plot Illustrating Cloud Coverage Over Time for Each of Kannenberg's 12 Sites",
        Xgrid = T)

# This function plots each site's time series as a new series within the same plot
ts_plot(Cloud_Coverage,
        line.mode = "markers",
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = TRUE,
        type = "single",
        Xtitle = "Date of Imaging",
        Ytitle = "Cloud Coverage (0 = none, 1 = total coverage)", 
        title = "A Time Series Plot Illustrating Cloud Coverage Over Time for Each of Kannenberg's 12 Sites",
        Xgrid = TRUE)

# This function plots a site individually
ts_plot(Cloud_Coverage[,c(1,3)], # Column 1 is the date, columns 3+ are specific sites
        line.mode = "lines+markers",
        dash = "dot",
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = TRUE,
        type = "single",
        Xtitle = "Date of Imaging",
        Ytitle = "Cloud Coverage (0 = none, 1 = total coverage)", 
        title = "A Time Series Plot Illustrating Cloud Coverage Over Time for Kannenberg Site: ",
        Xgrid = TRUE)

# This for loop plots each site's time series individually and automatically
for (i in 1:length(KBerg_AoI$Site)) {
  tmp.ts.plot <- ts_plot(Cloud_Coverage[,c(1, 2+i)], # Column 1 is the date, columns 3+ are specific sites
                         line.mode = "markers",
                         dash = "dot",
                         color = rainbow(length(KBerg_AoI$Site))[i],
                         slider = TRUE,
                         type = "single",
                         Xtitle = "Date of Imaging",
                         Ytitle = "Cloud Coverage (0 = none, 1 = total coverage)", 
                         title = paste0("A Time Series Plot Illustrating Cloud Coverage Over Time for Kannenberg Site: ", KBerg_AoI$Site[i]),
                         Xgrid = TRUE)
  print(tmp.ts.plot)
  
  # This nested for loop plots a vertical line for each unique instance of KBerg_Mortality_Data$`Survey date`
  for (j in 1:length(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])))
  {tmp.ts.plot <- tmp.ts.plot %>% 
    add_trace(
      type = "scatter",
      mode = "lines",
      x = unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j],
      y = c(0, 1),
      line = list(color = "black", width = 2),
      showlegend = FALSE
    )
  print(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j])}
  print(tmp.ts.plot)
}

### Charts for NDVI through time

# This function plots each site's time series in the same window
ts_plot(Mean_Site_NDVI,
        line.mode = "markers",
        dash = NULL,
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = F, # Slider is disabled for the "multiple" plot type
        type = "multiple",
        Xtitle = "Date of Imaging",
        Ytitle = "NDVI", 
        title = "A Time Series Plot Illustrating the Mean Sitewide NDVI Over Time for Each of Kannenberg's 12 Sites",
        Xgrid = T) %>%
  layout(xaxis = list(zeroline=T), 
         yaxis = list(range = list(-0.3, 0.8)))

# This function plots each site's time series as a new series within the same plot
ts_plot(Mean_Site_NDVI,
        line.mode = "markers",
        dash = NULL,
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = TRUE,
        type = "single",
        Xtitle = "Date of Imaging",
        Ytitle = "NDVI", 
        title = "A Time Series Plot Illustrating the Mean Sitewide NDVI Over Time for Each of Kannenberg's 12 Sites",
        Xgrid = TRUE) %>%
  layout(xaxis = list(zeroline=T), 
         yaxis = list(range = list(-0.3, 0.8))) #%>%
#add_trace(
#    type = "scatter",
#    mode = "lines",
#    x = unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j],
#    y = c(-0.3, 0.8),
#    line = list(color = "black", width = 2),
#    showlegend = FALSE
#  )

# This function plots a site individually
ts_plot(Mean_Site_NDVI[,c(2,3)], # Column 1 is the date, columns 3+ are specific sites
        line.mode = "markers",
        color = rainbow(length(KBerg_AoI$Site), alpha = 0.35),
        slider = TRUE,
        type = "single",
        Xtitle = "Date of Imaging",
        Ytitle = "NDVI", 
        title = "A Time Series Plot Illustrating the Mean Sitewide NDVI Over Time for Kannenberg Site: ",
        Xgrid = TRUE)

# This for loop plots each site's time series individually and automatically
for (i in 1:length(KBerg_AoI$Site)) {
  tmp.ts.NDVI.plot <- ts_plot(Mean_Site_NDVI[,c(2, 2+i)], # Column 2 is the date, columns 3+ are specific sites
                              line.mode = "markers",
                              dash = NULL,
                              color = rainbow(length(KBerg_AoI$Site))[i],
                              slider = TRUE,
                              type = "single",
                              Xtitle = "Date of Imaging",
                              Ytitle = "NDVI", 
                              title = paste0("A Time Series Plot Illustrating the Mean Sitewide NDVI Over Time for Kannenberg Site: ", KBerg_AoI$Site[i]),
                              Xgrid = TRUE) %>%
    layout(xaxis = list(zeroline=T), 
           yaxis = list(range = list(-0.3, 0.8)))
  
  # This nested for loop plots a vertical line for each unique instance of KBerg_Mortality_Data$`Survey date`
  for (j in 1:length(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])))
  {tmp.ts.NDVI.plot <- tmp.ts.NDVI.plot %>% 
    add_trace(
      type = "scatter",
      mode = "lines",
      x = unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j],
      y = c(-0.3, 0.8),
      line = list(color = "black", width = 2),
      showlegend = FALSE
    )
  print(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j])}
  print(tmp.ts.NDVI.plot)
}

# Make a multipanel plot showcasing individual site NDVI
ts.NDVI.Multiplot <- as.list(rep(NA, length(KBerg_AoI$Site)))

for (i in 1:length(KBerg_AoI$Site)) {
  ts.NDVI.Multiplot[[i]] <- ts_plot(Mean_Site_NDVI[,c(2, 2+i)], # Column 2 is the date, columns 3+ are specific sites
                                    line.mode = "markers",
                                    dash = NULL,
                                    color = rainbow(length(KBerg_AoI$Site))[i],
                                    slider = F,
                                    type = "single",
                                    Xtitle = "Date of Imaging",
                                    Ytitle = "NDVI", 
                                    Xgrid = TRUE) %>%
    layout(xaxis = list(zeroline=T), 
           yaxis = list(range = list(-0.3, 0.8)))
  
  # This nested for loop plots a vertical line for each unique instance of KBerg_Mortality_Data$`Survey date`
  for (j in 1:length(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])))
  {ts.NDVI.Multiplot[[i]] <- ts.NDVI.Multiplot[[i]] %>% 
    add_trace(
      type = "scatter",
      mode = "lines",
      x = unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j],
      y = c(-0.3, 0.8),
      line = list(color = "black", width = 2),
      showlegend = FALSE
    )
  print(unique(KBerg_Mortality_Data$`Survey date`[KBerg_Mortality_Data$Site == KBerg_AoI$Site[i]])[j])
  print(ts.NDVI.Multiplot[[i]])}
}

subplot(ts.NDVI.Multiplot,nrows = 12, shareX = TRUE) %>% 
  layout(title = list(text = "NDVI Time Series for Each Site"),
         plot_bgcolor='#e5ecf6', 
         xaxis = list( 
           zerolinecolor = '#ffff', 
           zerolinewidth = 2, 
           gridcolor = 'ffff'), 
         yaxis = list( 
           zerolinecolor = '#ffff', 
           range = list(-0.3, 0.8),
           zerolinewidth = 2, 
           gridcolor = 'ffff')) 

### Q3: Time series of multiple regressions (Can I predict the future or assess current forest health?) ####

### Approach using classic linear regressions
## This was the first approach I tried using linear regressions

# Make an object to store results
TSRegression_Results_SiteMeanDieback.MayGT <- list() # The naming syntax is "TSRegression" = Time Series Regression, "SiteMeanDieback" = all trees' canopy dieback, ".MayGT" = Using May ground truthing

## Segment my mean site NDVI dataframe into 10 day windows
# Create an empty list to store subsets
NDVI_Subsets <- list()

# Define the start and end dates for the 10-day windows
tmp.start_date <- min(Mean_Site_NDVI$Date)
tmp.end_date <- max(Mean_Site_NDVI$Date)

# Loop through each 10-day window
while (tmp.start_date <= tmp.end_date) {
  # Get the end date of the 10-day window
  tmp.window_end_date <- tmp.start_date + 9
  
  # Subset the dataframe for the current 10-day window
  tmp.subset_df <- Mean_Site_NDVI[Mean_Site_NDVI$Date >= tmp.start_date & Mean_Site_NDVI$Date <= tmp.window_end_date, ]
  
  # Append the subset to the list
  NDVI_Subsets[[paste(tmp.start_date, "-", tmp.window_end_date)]] <- tmp.subset_df
  
  # Update the start date for the next window
  tmp.start_date <- tmp.window_end_date + 1
}

# Now, NDVI_Subsets is a list containing subsets of Mean_Site_NDVI dataframe for every 10-day window
# You can access each subset using its name, e.g., NDVI_Subsets[["2022-01-01 - 2022-01-10"]] to access the subset for that window

## Take just one instance of each site from the NDVI subsets
## Run the model time series
# Loop through each NDVI subset
for (subset_name in names(NDVI_Subsets)) {
  # Extract the current NDVI subset
  tmp.current_subset <- NDVI_Subsets[[subset_name]]
  
  # Initialize a vector to store NDVI values for the current subset
  tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
  
  # Loop through each row of KBerg_Stand_Composition_May to get NDVI values from the current subset
  for (i in 1:nrow(KBerg_Stand_Composition_May)) {
    # Get the site name
    tmp.site <- KBerg_Stand_Composition_May$Site[i]
    
    # Extract the NDVI values for the current site from the current subset
    tmp.site_ndvi_values <- tmp.current_subset[[paste0(tmp.site, "_NDVI")]]
    
    # Calculate the mean NDVI for the site in the current subset
    if (length(tmp.site_ndvi_values) > 0) {
      tmp.site_NDVI[i] <- mean(tmp.site_ndvi_values, na.rm = TRUE)
    } else {
      tmp.site_NDVI[i] <- NA  # Set to NA if there are no NDVI values
    }
  }
  
  # Create a dataframe with the variables for the current subset
  tmp.LM_Data_TimeSeries_NDVI.SiteMeanDieback <- data.frame(
    Total_Mean_Dieback = KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent, 
    Mean_NDVI = tmp.site_NDVI
  )
  
  # Check if there are any non-NA values in Mean_NDVI
  if (sum(!is.na(tmp.LM_Data_TimeSeries_NDVI.SiteMeanDieback$Mean_NDVI)) > 0) {
    # Fit a linear regression model for the current subset
    tmp.LM_Model_TimeSeries_NDVI.SiteMeanDieback <- lm(
      Total_Mean_Dieback ~ Mean_NDVI, 
      data = tmp.LM_Data_TimeSeries_NDVI.SiteMeanDieback
    )
    
    # Store the regression model in the list with the subset name as the key
    TSRegression_Results_SiteMeanDieback.MayGT[[subset_name]] <- tmp.LM_Model_TimeSeries_NDVI.SiteMeanDieback
  } else {
    # If all values are NA, store a message instead of a model
    TSRegression_Results_SiteMeanDieback.MayGT[[subset_name]] <- "No non-NA NDVI values available"
  }
}

# TSRegression_Results_SiteMeanDieback.MayGT contains the linear regression models for each NDVI subset
# I can access each model using its subset name, e.g., TSRegression_Results_SiteMeanDieback.MayGT[["2022-01-01 - 2022-01-10"]]

## Put all of these data together into a dataframe for visualisation

# Initialize vectors to store the R-squared and p-values
tmp.r_squared_values <- numeric(length(TSRegression_Results_SiteMeanDieback.MayGT))
tmp.p_values <- numeric(length(TSRegression_Results_SiteMeanDieback.MayGT))
subset_names <- names(TSRegression_Results_SiteMeanDieback.MayGT)

# Loop through each result in the TSRegression_Results_SiteMeanDieback.MayGT list
for (i in seq_along(TSRegression_Results_SiteMeanDieback.MayGT)) {
  tmp.result <- TSRegression_Results_SiteMeanDieback.MayGT[[i]]
  
  if (is.character(tmp.result)) {
    # If the result is a message (e.g., "No non-NA NDVI values available"), set R-squared and p-value to NA
    tmp.r_squared_values[i] <- NA
    tmp.p_values[i] <- NA
  } else {
    # Extract R-squared value
    tmp.r_squared_values[i] <- summary(tmp.result)$r.squared
    
    # Extract p-value for the Mean_NDVI coefficient
    tmp.p_values[i] <- summary(tmp.result)$coefficients[2, 4]
  }
}

# Create a dataframe with the subset names, R-squared values, and p-values
TSRegression_SummaryDF_SiteMeanDieback.MayGT <- data.frame(
  Subset = subset_names,
  R_Squared = tmp.r_squared_values,
  P_Value = tmp.p_values
)

## Plot the R squareds and P values as a chart for observation

# Convert Subset to a date format if possible, otherwise treat as a factor
TSRegression_SummaryDF_SiteMeanDieback.MayGT$Subset <- as.Date(TSRegression_SummaryDF_SiteMeanDieback.MayGT$Subset, format = "%Y-%m-%d")
if (any(is.na(TSRegression_SummaryDF_SiteMeanDieback.MayGT$Subset))) {
  TSRegression_SummaryDF_SiteMeanDieback.MayGT$Subset <- as.factor(TSRegression_SummaryDF_SiteMeanDieback.MayGT$Subset)
}

# Create a graph of p-values over time
p_value_plot <- ggplot(data = TSRegression_SummaryDF_SiteMeanDieback.MayGT, aes(x = Subset, y = P_Value)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
  labs(title = "P-values Over Time", x = "Time Subset", y = "P-value") +
  theme_minimal()

# Create a graph of P-values and R^2 values over time
# Combined P value and R^2 plot:
P.R2_values_plot <- ggplot(data = TSRegression_SummaryDF_SiteMeanDieback.MayGT) +
  # P-value line and points (red, left y-axis)
  geom_line(aes(x = Subset, y = P_Value), color = "red", alpha = 0.3) +
  geom_point(aes(x = Subset, y = P_Value), color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  # R? line and points (blue, right y-axis)
  geom_line(aes(x = Subset, y = R_Squared), color = "blue", alpha = 0.3) +
  geom_point(aes(x = Subset, y = R_Squared), color = "blue", alpha = 0.7) +
  # Vertical reference line
  geom_vline(xintercept = as.Date("2019-05-07", format = "%Y-%m-%d"), linetype = "solid", color = "black") +
  # Axis labels and title
  labs(
    title = "P-values and R? Values of Linear Regressions Over Time",
    x = "Time Subset",
    y = "P-value"
  ) +
  # Secondary y-axis for Pseudo R?
  scale_y_continuous(
    name = "P-values (red points)",
    sec.axis = sec_axis(
      trans = ~ .,
      name = "R? values (blue points)"
    )
  ) +
  # Customized x-axis with yearly and quarterly tick marks
  scale_x_date(
    date_breaks = "1 year",  # Quarterly tick marks
    date_minor_breaks = "1 month"  # Minor tick marks for each month (optional)
    # Year and abbreviated month
  ) +
  theme_bw() +
  theme( # Legend code not working...
    legend.position = c(0.95, 0.95), # Top-right corner
    legend.justification = "top right",
    legend.background = element_rect(fill = "white", color = "black", size = 0.5))

# Display the plot
print(p_value_plot) # I visually counted 252 points (I might've over/undercounted by a couple)
print(P.R2_values_plot)

### Q3: Time series using Beta regressions (Can I predict the future or assess canopy greening in the present?) ####

### Approach using beta regressions to account for the fact that my data is bounded
## This approach is instead of using linear regressions

# Make an object to store results
TSBetaReg_Results_SiteMeanDieback.MayGT <- list() # The naming syntax is "TSBetaReg" = Time Series Beta Regression, "SiteMeanDieback" = all trees' canopy dieback, ".MayGT" = Using May ground truthing

## Segment my mean site NDVI dataframe into 10 day windows
# Create an empty list to store subsets
NDVI_Subsets <- list()

# Define the start and end dates for the 10-day windows
tmp.start_date <- min(Mean_Site_NDVI$Date)
tmp.end_date <- max(Mean_Site_NDVI$Date)

# Loop through each 10-day window
while (tmp.start_date <= tmp.end_date) {
  # Get the end date of the 10-day window
  tmp.window_end_date <- tmp.start_date + 9
  
  # Subset the dataframe for the current 10-day window
  tmp.subset_df <- Mean_Site_NDVI[Mean_Site_NDVI$Date >= tmp.start_date & Mean_Site_NDVI$Date <= tmp.window_end_date, ]
  
  # Append the subset to the list
  NDVI_Subsets[[paste(tmp.start_date, "-", tmp.window_end_date)]] <- tmp.subset_df
  
  # Update the start date for the next window
  tmp.start_date <- tmp.window_end_date + 1
}

# Now, NDVI_Subsets is a list containing subsets of Mean_Site_NDVI dataframe for every 10-day window
# You can access each subset using its name, e.g., NDVI_Subsets[["2022-01-01 - 2022-01-10"]] to access the subset for that window

## Take just one instance of each site from the NDVI subsets
## Run the model time series
# Loop through each NDVI subset
for (subset_name in names(NDVI_Subsets)) {
  # Extract the current NDVI subset
  tmp.current_subset <- NDVI_Subsets[[subset_name]]
  
  # Initialize a vector to store NDVI values for the current subset
  tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
  
  # Loop through each row of KBerg_Stand_Composition_May to get NDVI values from the current subset
  for (i in 1:nrow(KBerg_Stand_Composition_May)) {
    # Get the site name
    tmp.site <- KBerg_Stand_Composition_May$Site[i]
    
    # Extract the NDVI values for the current site from the current subset
    tmp.site_ndvi_values <- tmp.current_subset[[paste0(tmp.site, "_NDVI")]]
    
    # Calculate the mean NDVI for the site in the current subset
    if (length(tmp.site_ndvi_values) > 0) {
      tmp.site_NDVI[i] <- mean(tmp.site_ndvi_values, na.rm = TRUE)
    } else {
      tmp.site_NDVI[i] <- NA  # Set to NA if there are no NDVI values
    }
  }
  
  # Create a dataframe with the variables for the current subset
  tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback <- data.frame(
    Total_Mean_Dieback = KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent/100, # Converted to decimal for Beta Regression 
    Mean_NDVI = tmp.site_NDVI
  )
  
  # Check if there are any non-NA values in Mean_NDVI
  if (sum(!is.na(tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback$Mean_NDVI)) > 0) {
    # Fit a beta regression model for the current subset
    tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback <- betareg(
      Total_Mean_Dieback ~ Mean_NDVI, 
      data = tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback
    )
    
    # Store the beta regression model in the list with the subset name as the key
    TSBetaReg_Results_SiteMeanDieback.MayGT[[subset_name]] <- tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback
  } else {
    # If all values are NA, store a message instead of a model
    TSBetaReg_Results_SiteMeanDieback.MayGT[[subset_name]] <- "No non-NA NDVI values available"
  }
}

# TSBetaReg_Results_SiteMeanDieback.MayGT contains the Beta Regression models for each NDVI subset
# I can access each model using its subset name, e.g., TSBetaReg_Results_SiteMeanDieback.MayGT[["2022-01-01 - 2022-01-10"]]

## Put all of these data together into a dataframe for visualisation

# Initialize vectors to store the Pseudo-R-squared and p-values
tmp.pseudo_r_squared_values <- numeric(length(TSBetaReg_Results_SiteMeanDieback.MayGT))
tmp.p_values <- numeric(length(TSBetaReg_Results_SiteMeanDieback.MayGT))
subset_names <- names(TSBetaReg_Results_SiteMeanDieback.MayGT)

# Loop through each result in the TSBetaReg_Results_SiteMeanDieback.MayGT list
for (i in seq_along(TSBetaReg_Results_SiteMeanDieback.MayGT)) {
  tmp.result <- TSBetaReg_Results_SiteMeanDieback.MayGT[[i]]
  
  if (is.character(tmp.result)) {
    # If the result is a message (e.g., "No non-NA NDVI values available"), set R-squared and p-value to NA
    tmp.pseudo_r_squared_values[i] <- NA
    tmp.p_values[i] <- NA
  } else {
    # Extract R-squared value
    tmp.pseudo_r_squared_values[i] <- summary(tmp.result)$pseudo.r.squared
    
    # Extract p-value for the Mean_NDVI coefficient
    tmp.p_values[i] <- summary(tmp.result)$coefficients$precision[1, 4]
  }
}

# Create a dataframe with the subset names, R-squared values, and p-values
TSBetaReg_SummaryDF_SiteMeanDieback.MayGT <- data.frame(
  Subset = subset_names,
  Pseudo_R_Squared = tmp.pseudo_r_squared_values,
  P_Value = tmp.p_values
)

## Plot the pseudo-R-squareds and P values as a chart for observation

# Convert Subset to a date format if possible, otherwise treat as a factor
TSBetaReg_SummaryDF_SiteMeanDieback.MayGT$Subset <- as.Date(TSBetaReg_SummaryDF_SiteMeanDieback.MayGT$Subset, format = "%Y-%m-%d")
if (any(is.na(TSBetaReg_SummaryDF_SiteMeanDieback.MayGT$Subset))) {
  TSBetaReg_SummaryDF_SiteMeanDieback.MayGT$Subset <- as.factor(TSBetaReg_SummaryDF_SiteMeanDieback.MayGT$Subset)
}

# Plot p-values over time
P_value_plot <- ggplot(data = TSBetaReg_SummaryDF_SiteMeanDieback.MayGT, aes(x = Subset, y = P_Value)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
  labs(title = "P-values of Beta Regressions Over Time", x = "Time Subset", y = "P-value") +
  theme_minimal()

# Display the plot
print(P_value_plot) # I visually counted 252 points (I might've over/undercounted by a couple)

### Q3: Time series using Beta regressions with Juniper proportion (Can I predict the future or assess current canopy greening?) ####

### Approach using beta regressions to account for the fact that my data is bounded
### And Juniper as a proportion of the site because it previously increased model performance
## This will be the definitive version of this statistical test which I'll publish
## This approach will include the effects of stand composition via Dirk's suggested NDVI weighting
## But I will then also run the model with Juniper-prop as a separate interaction term

# Make an object to store results
TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT <- list() # The naming syntax is "TSBetaReg" = Time Series Beta Regression, "SiteMeanDieback" = all trees' canopy dieback, ".MayGT" = Using May ground truthing

## Segment my mean site NDVI dataframe into 10 day windows
# Create an empty list to store subsets
NDVI_Subsets <- list()

# Define the start and end dates for the 10-day windows
tmp.start_date <- min(Mean_Site_NDVI$Date)
tmp.end_date <- max(Mean_Site_NDVI$Date)

# Loop through each 10-day window
while (tmp.start_date <= tmp.end_date) {
  # Get the end date of the 10-day window
  tmp.window_end_date <- tmp.start_date + 9
  
  # Subset the dataframe for the current 10-day window
  tmp.subset_df <- Mean_Site_NDVI[Mean_Site_NDVI$Date >= tmp.start_date & Mean_Site_NDVI$Date <= tmp.window_end_date, ]
  
  # Append the subset to the list
  NDVI_Subsets[[paste(tmp.start_date, "-", tmp.window_end_date)]] <- tmp.subset_df
  
  # Update the start date for the next window
  tmp.start_date <- tmp.window_end_date + 1
}

# Now, NDVI_Subsets is a list containing subsets of Mean_Site_NDVI dataframe for every 10-day window
# You can access each subset using its name, e.g., NDVI_Subsets[["2022-01-01 - 2022-01-10"]] to access the subset for that window

## Take just one instance of each site from the NDVI subsets
## Run the model time series

## THIS IS THE VERSION OF THE CODE WHICH USES NDVI WEIGHTING
# Loop through each NDVI subset
for (subset_name in names(NDVI_Subsets)) {
  # Extract the current NDVI subset
  tmp.current_subset <- NDVI_Subsets[[subset_name]]
  
  # Initialize a vector to store NDVI values for the current subset
  tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
  
  # I'm also initialising this new vector below to hold Juniper proportion adjusted NDVI values for the current subset
  tmp.site_adjusted_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
  
  # Loop through each row of KBerg_Stand_Composition_May to get NDVI values from the current subset
  for (i in 1:nrow(KBerg_Stand_Composition_May)) {
    # Get the site name
    tmp.site <- KBerg_Stand_Composition_May$Site[i]
    
    # Extract the NDVI values for the current site from the current subset
    tmp.site_ndvi_values <- tmp.current_subset[[paste0(tmp.site, "_NDVI")]]
    
    # Calculate the mean NDVI for the site in the current subset
    if (length(tmp.site_ndvi_values) > 0) {
      
      ## This is where this version of the code deviates from the previous version
      # Instead of just taking the mean, it also multiplies it by the proportion of Junipers in the site
      tmp.site_NDVI[i] <- mean(tmp.site_ndvi_values, na.rm = TRUE) # Here, tmp.site_NDVI[i] is the site's mean NDVI for the time period being calculated
      # The line of code below is new:
      tmp.site_adjusted_NDVI[i] <- tmp.site_NDVI[i] * (1 - (KBerg_Stand_Composition_May[i,]$Juniper_Site_Percent/100)) # composition converted to decimal, decimal multiplied by tmp.site_NDVI
    } else {
      tmp.site_NDVI[i] <- NA  # Set to NA if there are no NDVI values
      tmp.site_adjusted_NDVI[i] <- NA # New line for this analysis
    }
  }
  
  # Create a dataframe with the variables for the current subset
  tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback <- data.frame(
    Total_Mean_Dieback = KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent/100, # Converted to decimal for Beta Regression 
    Mean_NDVI = tmp.site_NDVI,
    Adjusted_Mean_NDVI = tmp.site_adjusted_NDVI
  )
  
  # Check if there are any non-NA values in Mean_NDVI
  if (sum(!is.na(tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback$Adjusted_Mean_NDVI)) > 0) {
    # Fit a beta regression model for the current subset
    tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback <- betareg(
      Total_Mean_Dieback ~ Adjusted_Mean_NDVI, 
      data = tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback
    )
    
    # Store the beta regression model in the list with the subset name as the key
    TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[[subset_name]] <- tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback
  } else {
    # If all values are NA, store a message instead of a model
    TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[[subset_name]] <- "No non-NA NDVI values available"
  }
}

## THIS IS THE VERSION OF THE CODE WHICH USES Juniper proportion as a separate interaction term

tmp.loop.counter <- 0 # Because I was debugging the code quite a lot, I wanted this to keep track of how it's going

# Loop through each NDVI subset
for (subset_name in names(NDVI_Subsets)) {
  # Extract the current NDVI subset
  tmp.current_subset <- NDVI_Subsets[[subset_name]]
  
  # Initialize a vector to store NDVI values for the current subset
  tmp.site_NDVI <- numeric(length(KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent))
  
  # Loop through each row of KBerg_Stand_Composition_May to get NDVI values from the current subset
  for (i in 1:nrow(KBerg_Stand_Composition_May)) {
    # Get the site name
    tmp.site <- KBerg_Stand_Composition_May$Site[i]
    
    # Extract the NDVI values for the current site from the current subset
    tmp.site_ndvi_values <- tmp.current_subset[[paste0(tmp.site, "_NDVI")]]
    
    # Calculate the mean NDVI for the site in the current subset
    if (length(tmp.site_ndvi_values) > 0) {
      tmp.site_NDVI[i] <- mean(tmp.site_ndvi_values, na.rm = TRUE)
    } else {
      tmp.site_NDVI[i] <- NA  # Set to NA if there are no NDVI values
    }
  }
  
  # Create a dataframe with the variables for the current subset
  tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback <- data.frame(
    Total_Mean_Dieback = KBerg_Stand_Composition_May$Total_Mean_Dieback_Percent/100, # Converted to decimal for Beta Regression 
    Mean_NDVI = tmp.site_NDVI,
    Juniper_Site_Proportion = KBerg_Stand_Composition_May$Juniper_Site_Percent
  )
  tmp.loop.counter <- tmp.loop.counter+1
  print(tmp.loop.counter)
  print(subset_name)
  
  # Check if there are any non-NA values in Mean_NDVI
  if (sum(!is.na(tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback$Mean_NDVI)) > 0) {
    # Fit a beta regression model for the current subset
    tryCatch({
      tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback <- betareg(
        Total_Mean_Dieback ~ Mean_NDVI * Juniper_Site_Proportion, 
        data = tmp.BetaReg_Data_TimeSeries_NDVI.SiteMeanDieback
      )
      TSBetaReg_Results_SiteMeanDieback.MayGT[[subset_name]] <- tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback
    }, error = function(e) {
      TSBetaReg_Results_SiteMeanDieback.MayGT[[subset_name]] <- paste("Model failed:", e$message)
    })
    
    # Store the beta regression model in the list with the subset name as the key
    TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[[subset_name]] <- tmp.BetaReg_Model_TimeSeries_NDVI.SiteMeanDieback
  } else {
    # If all values are NA, store a message instead of a model
    TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[[subset_name]] <- "No non-NA NDVI values available"
  }
}

tmp.loop.counter <- 0 

## TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT contains the Beta Regression models for each NDVI subset
# I can access each model using its subset name, e.g., TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[["2022-01-01 - 2022-01-10"]]

## Put all of these data together into a dataframe for visualisation

# Initialize vectors to store the Pseudo-R-squared and p-values
tmp.pseudo_r_squared_values <- numeric(length(TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT))
tmp.p_values <- numeric(length(TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT))
subset_names <- names(TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT)

# Loop through each result in the TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT list
for (i in seq_along(TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT)) {
  tmp.result <- TSBetaReg_Results_SiteMeanDieback.JuniperProp.MayGT[[i]]
  
  if (is.character(tmp.result)) {
    # If the result is a message (e.g., "No non-NA NDVI values available"), set R-squared and p-value to NA
    tmp.pseudo_r_squared_values[i] <- NA
    tmp.p_values[i] <- NA
  } else {
    # Extract R-squared value
    tmp.pseudo_r_squared_values[i] <- summary(tmp.result)$pseudo.r.squared
    
    # Extract p-value for the Mean_NDVI coefficient
    tmp.p_values[i] <- summary(tmp.result)$coefficients$precision[1, 4]
  }
}

# Create a dataframe with the subset names, R-squared values, and p-values
TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT <- data.frame(
  Subset = subset_names,
  Pseudo_R_Squared = tmp.pseudo_r_squared_values,
  P_Value = tmp.p_values
)

## Plot the pseudo-R-squareds and P values as a chart for observation

# Convert Subset to a date format if possible, otherwise treat as a factor
TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT$Subset <- as.Date(TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT$Subset, format = "%Y-%m-%d")
if (any(is.na(TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT$Subset))) {
  TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT$Subset <- as.factor(TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT$Subset)
}

# Plot p-values over time
P_value_plot <- ggplot(data = TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT, aes(x = Subset, y = P_Value)) +
  # coord_cartesian(xlim = (as.Date(c("2018-06-01", "2023-06-01"), format = "%Y-%m-%d")) + ylim = c(0,1)) +
  geom_line(color = "red", alpha=0.3) +
  geom_point(color = "red", alpha=0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = as.Date("2019-05-07", format = "%Y-%m-%d"), linetype = "solid", color = "black") +
  labs(title = "P-values of Beta Regressions Over Time", x = "Time Subset", y = "P-value") +
  theme_bw()

# Plot R^2 values over time
R2_value_plot <- ggplot(data = TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT, aes(x = Subset, y = Pseudo_R_Squared)) +
  geom_line(color = "blue", alpha=0.3) +
  geom_point(color = "blue", alpha=0.7) +
  geom_vline(xintercept = as.Date("2019-05-07", format = "%Y-%m-%d"), linetype = "solid", color = "black") +
  labs(title = "Pseudo R? values of Beta Regressions Over Time", x = "Time Subset", y = "Pseudo R?") +
  theme_bw()

# Combined P value and R^2 plot:
P.R2_values_plot <- ggplot(data = TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT) +
  # P-value line and points (red, left y-axis)
  geom_line(aes(x = Subset, y = P_Value), color = "red", alpha = 0.3) +
  geom_point(aes(x = Subset, y = P_Value), shape = "X", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  # R? line and points (blue, right y-axis)
  geom_line(aes(x = Subset, y = Pseudo_R_Squared), color = "blue", alpha = 0.3) +
  geom_point(aes(x = Subset, y = Pseudo_R_Squared), shape = "X", color = "blue", alpha = 0.7) +
  # Vertical reference line
  geom_vline(xintercept = as.Date("2019-05-07", format = "%Y-%m-%d"), linetype = "solid", color = "black") +
  # Axis labels and title
  labs(
    title = "P-values and Pseudo R? Values of Beta Regressions Over Time",
    x = "Time Subset",
    y = "P-value"
  ) +
  # Secondary y-axis for Pseudo R?
  scale_y_continuous(
    name = "P-values (red points)",
    sec.axis = sec_axis(
      trans = ~ .,
      name = "Pseudo R? values (blue points)"
    )
  ) +
  # Customized x-axis with yearly and quarterly tick marks
  scale_x_date(
    date_breaks = "1 year",  # Quarterly tick marks
    date_minor_breaks = "1 month"  # Minor tick marks for each month (optional)
      # Year and abbreviated month
  ) +
  theme_bw() #+
  #theme( # Legend code not working...
   # legend.position = c(0.95, 0.95), # Top-right corner
    #legend.justification = "top right",
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5))

# Display the plot
print(P_value_plot) 
print(R2_value_plot)
print(P.R2_values_plot)

## Final version for manuscript (Figure 4)
## Pseudo R² plot (now appears separately and on top)
p_combined_top <- ggplot(TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT, aes(x = Subset, y = Pseudo_R_Squared)) +
  geom_line(color = "blue", alpha = 0.4) +
  geom_point(aes(x = Subset, y = Pseudo_R_Squared), shape = "•", size = 2.5, color = "blue", alpha = 0.85) +
  geom_vline(xintercept = as.Date("2019-05-07"), linetype = "solid", color = "black") +
  labs(
    title = "Pseudo-R² Values of Beta Regression Models Over Time",
    x = NULL,
    y = expression(paste("Pseudo-", R^2))
  ) +
  scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), plot.margin = margin(10, 35, 0, 0)) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 19),
    axis.title.y = element_text(size = 18))

## P-value plot (now separated and on the bottom)
p_combined_bottom <- ggplot(TSBetaReg_SummaryDF_SiteMeanDieback.JuniperProp.MayGT, aes(x = Subset, y = P_Value)) +
  geom_line(color = "red", alpha = 0.4) +
  geom_point(aes(x = Subset, y = P_Value), shape = "•", size = 2.5, color = "red", alpha = 0.85) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  geom_vline(xintercept = as.Date("2019-05-07"), linetype = "solid", color = "black") +
  labs(
    title = "P-values of Beta Regression Models Over Time",
    x = "Time Subset",
    y = "P-value"
  ) +
  scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), plot.margin = margin(10, 35, 0, 0))+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 19),
    axis.title.y = element_text(size = 18))

# Combine and display: R² on top, P-value on bottom
(p_combined_top / p_combined_bottom) + plot_layout(heights = c(1, 1))

### Q3 Continued: How is relative canopy mortality progressing? ####
## I will use my SpeciesComp.NDVI_BetaReg_Model object as the model and NDVI_Subsets as my NDVI time series

# First I will initialise a list to store the results of my predictions
CanopyDieback_Timeseries_Predictions <- list()

# Get Juniper proportions (assumed fixed per site)
juniper.may.site.props <- KBerg_Stand_Composition_May$Juniper_Site_Percent
names(juniper.may.site.props) <- KBerg_Stand_Composition_May$Site

# Loop through each NDVI interval
for (subset_name in names(NDVI_Subsets)) {
  tmp.subset <- NDVI_Subsets[[subset_name]]
  
  # Initialize vector to store predicted dieback
  tmp.predictions <- numeric(length(juniper.may.site.props))
  names(tmp.predictions) <- names(juniper.may.site.props)
  
  for (site_name in names(juniper.may.site.props)) {
    tmp.NDVI.site <- paste0(site_name, "_NDVI")
    
    if (tmp.NDVI.site %in% names(tmp.subset)) {
      tmp.NDVI.values <- tmp.subset[[tmp.NDVI.site]]
      tmp.sitewide.timewide.NDVI.mean <- mean(tmp.NDVI.values, na.rm = TRUE)
      
      if (!is.na(tmp.sitewide.timewide.NDVI.mean)) {
        # Create a new data frame with predictor values
        tmp.prediction_input <- data.frame(
          Mean_NDVI = tmp.sitewide.timewide.NDVI.mean,
          Juniper_Site_Percent = juniper.may.site.props[[site_name]]
        )
        
        # Predict using the beta regression model
        tmp.predicted.value <- predict(SpeciesComp.NDVI_BetaReg_Model, newdata = tmp.prediction_input, type = "response")
        tmp.predictions[[site_name]] <- tmp.predicted.value
      } else {
        tmp.predictions[[site_name]] <- NA
      }
    } else {
      tmp.predictions[[site_name]] <- NA
    }
  }
  
  # Store predictions for this time slice
  CanopyDieback_Timeseries_Predictions[[subset_name]] <- tmp.predictions
}

# Now I will put these results into a single dataframe for easier reading and plotting
CanopyDieback_Timeseries_Predictions_DF <- do.call(rbind, CanopyDieback_Timeseries_Predictions)

# Convert rownames (time intervals) into a proper column
CanopyDieback_Timeseries_Predictions_DF <- data.frame(
  Interval_Start = as.Date(sub(" -.*", "", rownames(CanopyDieback_Timeseries_Predictions_DF)), format = "%Y-%m-%d"),
  CanopyDieback_Timeseries_Predictions_DF,
  row.names = NULL
)

## Now plot the relative canopy dieback predictions through time (Figure 5)

# Reshape the data: pivot longer for ggplot
CanopyDieback_Timeseries_Predictions_Long <- CanopyDieback_Timeseries_Predictions_DF %>%
  pivot_longer(
    cols = -Interval_Start,
    names_to = "Site",
    values_to = "Predicted_Dieback")

# Create a colourblind friendly colour palette 
colourblind_colours <- c(
  "AR1" = "#990000",  # dark red
  "AR2" = "#CC3333",  # medium red
  "AR3" = "#E66A6A",  # lighter red
  "AR4" = "#FF9999",  # pale red (pink)
  
  "CM1" = "#003366",  # dark blue
  "CM2" = "#336699",  # medium blue
  "CM3" = "#6699CC",  # lighter blue
  "CM4" = "#99CCFF",  # pale blue
  
  "MD1" = "#004d00",  # dark green
  "MD2" = "#339933",  # medium green
  "MD3" = "#66CC66",  # lighter green
  "MD4" = "#99FF99"   # pale green
)

ggplot(CanopyDieback_Timeseries_Predictions_Long, aes(x = Interval_Start, y = Predicted_Dieback, color = Site)) +
  geom_line(alpha = 0.4, size = 0.5) +  # Thin, faint lines
  geom_point(alpha = 0.8, size = 1.5) +  # Slightly larger points for clarity
  
  scale_color_manual(values = colourblind_colours) +
  
  labs(
    title = "Predicted Canopy Dieback Over Time",
    subtitle = "Site-level predictions from Beta Regression Model",
    x = "Date",
    y = "Sitewide Predicted Relative Canopy Dieback (%)",
    color = "Site",
    ylim = c(0, NA)
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(scale = 100)
  ) +
  
  scale_x_date(
    date_breaks = "3 months",
    date_minor_breaks = "1 month",
    date_labels = "%Y-%b"
  ) +
  
  theme_bw(base_size = 16) +  # Increase base font size for clarity at print scale
  
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.position = c(0.83, 0.98),  # Inside top-right corner
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"),
    legend.box = "vertical",
    legend.box.just = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  
  guides(color = guide_legend(ncol = 3))  # 3-column legend

# Save the plot at a size which fits the MDPI template's usable area
ggsave("Figure5_V2.png", width = 18.46, height = 23.6, units = "cm", dpi = 900)

### Q4: Can these sites predict other sites in the forest? (Lateral Analysis Approach 1: Leave One Out Cross Validation) ####

### The bulk of this code came from ChatGPT. Consequently, there's a couple of things to note:
## ChatGPT is only using May data. 
## ChatGPT created each LM whilst excluding one of the sites, and is then testing the LM's prediction against the excluded site
## Why is that, tho? Shouldn't I use all available data to create the best model I have, and then test its predictive power?
# ^^^ I subsequently resolved all of these queries and put my thoughts/learnings in the notes

# Initialize vectors to store results
CV_Predicted_Dieback <- numeric(nrow(KBerg_Stand_Composition_May.Oct)) # just a blank, empty vector for subsequent storage
CV_Actual_Dieback <- KBerg_Stand_Composition_May.Oct$Total_Mean_Dieback_Percent # This is the actual data we're comparing against [Why just May?]

# Loop through each site for leave-one-out cross-validation
for (i in 1:nrow(KBerg_Stand_Composition_May.Oct)) {
  # Exclude the current site
  CV_Training_Data <- KBerg_Stand_Composition_May.Oct[-i, ] # Creating a Dataframe of training data without the current site
  CV_Test_Data <- KBerg_Stand_Composition_May.Oct[i, ] # Creating a Dataframe with only the current site to compare against
  
  # Extract the corresponding NDVI values (These nested for loops are understandable but complex)
  CV_Training_NDVI <- numeric(length(CV_Training_Data$Total_Mean_Dieback_Percent)) # Vector creation
  CV_Test_NDVI <- numeric(1)
  
  # The for loop below populates each row of the CV_Training_NDVI Dataframe iteratively as it loops
  # It works by looking at the sample date for each site, and then finding the nearest available date with NDVI data
  # So it doesn't do much statistically, just makes the Dataframe used for model training
  
  for (j in 1:nrow(CV_Training_Data)) {
    tmp.site <- CV_Training_Data$Site[j] # Each "site" corresponds to a site from the KBerg dataset. This loops through each site.
    if(grepl("_May", tmp.site)) # This finds the dieback sampling date (ground truthing date)
    {tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site, end=-5)][1]} 
    else {tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site, end=-5)][1]}
    tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]]),] # In my Sites' mean NDVI Dataframe, this takes the dates which have non-NA NDVI values for site j
    tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),] # Finds the nearest date with NDVI to the ground truthing date
    CV_Training_NDVI[j] <- Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID] # Takes the NDVI value and puts it in the Dataframe
  }
  
  # This basically just works like the loop above, but because only one site is being left out, it needn't loop or iterate
  
  tmp.site.2 <- CV_Test_Data$Site
  if(grepl("_May", tmp.site.2)) 
  {tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site.2, end=-5)][1]} 
  else{tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site.2, end=-5)][1]}
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site.2, end=-5), "_NDVI")]]),]
  tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),]
  CV_Test_NDVI <- Mean_Site_NDVI[[paste0(str_sub(tmp.site.2, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID] # This finds the NDVI of the site the model will be tested on
  
  # Create training data frame
  # This is of the sites included in building the model
  tmp.CV.training.df <- data.frame(Total_Mean_Dieback = CV_Training_Data$Total_Mean_Dieback_Percent, Mean_NDVI = CV_Training_NDVI)
  
  # Train the linear regression model
  tmp.CV.linear.model <- lm(Total_Mean_Dieback ~ Mean_NDVI, data = tmp.CV.training.df)
  
  ## Predict the canopy dieback for the test site
  # This is the main part of the analysis. Each LM predicts a canopy dieback amount for the missing site. The output is a single numeric variable.
  # The predicted dieback value for each site left out is put into this vector for subsequent model evaluation with the actual data
  CV_Predicted_Dieback[i] <- predict(tmp.CV.linear.model, newdata = data.frame(Mean_NDVI = CV_Test_NDVI))
}

## Assess the results
## This includes site MD4 May (which is an outlier)

# Calculate performance metrics
CV_MAE <- mean(abs(CV_Predicted_Dieback - CV_Actual_Dieback)) # Mean Absolute Error (the average disparity between observed data and model predictions)
CV_MSE <- mean((CV_Predicted_Dieback - CV_Actual_Dieback)^2) # Mean Squared Error
CV_RMSE <- sqrt(CV_MSE) # Root Mean Squared Error
CV_R_Squared <- cor(CV_Predicted_Dieback, CV_Actual_Dieback)^2

# Print performance metrics
cat("Mean Absolute Error (MAE):", CV_MAE, "\n")
cat("Mean Squared Error (MSE):", CV_MSE, "\n")
print(CV_RMSE)
cat("R-squared:", CV_R_Squared, "\n")

# Plot predicted vs actual dieback with labels
plot(CV_Actual_Dieback, CV_Predicted_Dieback, xlab = "Actual Canopy Dieback", ylab = "Predicted Canopy Dieback", xlim = c(0,110), ylim = c(0,110),
     main = "Leave-One-Out Cross-Validation", pch = 16, col = "blue")
abline(0, 1, col = "red")  # Add a y=x reference line. Positive results would cling closely to the line

# Add labels (assuming KBerg_Stand_Composition_May$Site contains site names)
text(CV_Actual_Dieback, CV_Predicted_Dieback, labels = KBerg_Stand_Composition_May.Oct$Site, pos = 2, cex = 0.8, col = "black")

## This excludes site MD4 May (which is an outlier)

# Calculate performance metrics but excluding the site MD4 outlier
CV_MAE <- mean(abs(CV_Predicted_Dieback[-12] - CV_Actual_Dieback[-12])) # Mean Absolute Error (the average disparity between observed data and model predictions)
CV_MSE <- mean((CV_Predicted_Dieback[-12] - CV_Actual_Dieback[-12])^2) # Mean Squared Error
CV_RMSE <- sqrt(CV_MSE) # Root Mean Squared Error
CV_R_Squared <- cor(CV_Predicted_Dieback[-12], CV_Actual_Dieback[-12])^2

# Print performance metrics
cat("Mean Absolute Error (MAE):", CV_MAE, "\n")
cat("Mean Squared Error (MSE):", CV_MSE, "\n")
cat("R-squared:", CV_R_Squared, "\n")

# Plot predicted vs actual dieback with labels
plot(CV_Actual_Dieback[-12], CV_Predicted_Dieback[-12], xlab = "Actual Canopy Dieback", ylab = "Predicted Canopy Dieback", xlim = c(0,100), ylim = c(0,100),
     main = "Leave-One-Out Cross-Validation", pch = 16, col = "blue")
abline(0, 1, col = "red")  # Add a y=x reference line. Positive results would cling closely to the line

# Add labels (assuming KBerg_Stand_Composition_May$Site contains site names)
text(CV_Actual_Dieback[-12], CV_Predicted_Dieback[-12], labels = KBerg_Stand_Composition_May.Oct$Site[-12], pos = 2, cex = 0.8, col = "black")

### Q4: Can these sites predict other sites in the forest? (Lateral Analysis Approach 2: Spatial Cross Validation) ####

## Step 1: Prepare the data for analysis
# Load my dataset
SCV_Data <- KBerg_Stand_Composition_May.Oct

# These data don't have coordinates, so I'll add that
# I will use co-ordinate data from my KBerg AoI dataframe
# I'm doing this simply. If I change the site order, then the co-ordinates being appended will be mismatched
SCV_Data <- bind_cols(SCV_Data, c(KBerg_AoI$Latitude, KBerg_AoI$Latitude)) # add latitude c(repeated because I need 2*12 occurrences)
colnames(SCV_Data)[10] <- "Latitude"
SCV_Data <- bind_cols(SCV_Data, c(KBerg_AoI$Longitude, KBerg_AoI$Longitude)) # add longitude c(repeated because I need 2*12 occurrences)
colnames(SCV_Data)[11] <- "Longitude"

# I also need to add NDVI data
SCV_Data <- bind_cols(SCV_Data, LM_Data_May.Oct_NDVI.SiteMeanDieback$Mean_NDVI)
colnames(SCV_Data)[12] <- "Mean_NDVI"

# Now that my data has latitude and longitude, convert it to a spatial object
sp::coordinates(SCV_Data) <- ~Longitude + Latitude  # Pointing out the column names for coordinates
SCV_Data_sf <- st_as_sf(SCV_Data, coords = c("Longitude", "Latitude"), crs = crs(KBerg_AoI)) # CRS was matched to the CRS of KBerg_AoI

## Step 2: Create spatial blocks

# Extract coordinates for each block based on site names
block1_sites <- c("AR1", "AR2", "AR3", "AR4")
block2_sites <- c("CM1", "CM2", "CM3", "CM4")
block3_sites <- c("MD1", "MD2", "MD3", "MD4")

# Filter the coordinates for each block
block1_coords <- `KBerg_co-ords` %>% filter(Site %in% block1_sites)
block2_coords <- `KBerg_co-ords` %>% filter(Site %in% block2_sites)
block3_coords <- `KBerg_co-ords` %>% filter(Site %in% block3_sites)

# Function to create a convex hull polygon for a given set of points
create_convex_hull <- function(points) {
  st_convex_hull(st_union(points))
}

# Generate polygons for each block
block1_polygon <- create_convex_hull(block1_coords)
block2_polygon <- create_convex_hull(block2_coords)
block3_polygon <- create_convex_hull(block3_coords)

# Combine into a single sf object with an identifier for each block

# My crude attempt at combining into one object
spatial_blocks <- c(block1_polygon, block2_polygon, block3_polygon)

# Create an sf data frame with the block identifiers
spatial_blocks <- st_sf(
  Block = c("Block1", "Block2", "Block3"),  # Assign block names
  geometry = spatial_blocks,                # Use the sfc_POLYGON object
  crs = st_crs(`KBerg_co-ords`)             # Set the CRS to match your KBerg_co-ords
)
spatial_blocks$Fold <- 1:3  # Assign fold numbers 1, 2, and 3 to each block

# Check the result
View(spatial_blocks)

# Check the result
plot(st_geometry(spatial_blocks), col = c("red", "blue", "green"), main = "Manual Spatial Blocks")
plot(st_geometry(`KBerg_co-ords`), add = TRUE, pch = 19, col = "black")

# Ensure both spatial_blocks_sf and SCV_Data_sf have the same CRS
st_crs(SCV_Data_sf) <- st_crs(spatial_blocks_sf)

# Note that this function is deprecated and I should upgrade to the new function below
set.seed(42) # For reproducibility
spatial_folds_old <- spatialBlock(
  speciesData = SCV_Data_sf,
  theRange = 5000,  # Adjust spatial range as needed (e.g., 5000 meters)
  k = 3,            # Choose the number of folds (3, as that's the number of regions which KBerg sampled)
  selection = "random",
  iteration = 100,  # Number of iterations to stabilize folds
  verbose = TRUE
)
# Display the created folds
#plot(spatial_folds_old)

# The new function for creating spatial folds
spatial_folds <- cv_spatial(
  SCV_Data_sf, # My "x" argument providing the data and geographic info
  k = 3,
  #selection = 'predefined',
  selection = 'random',
  iteration = 10000, # switch
  #user_blocks = spatial_blocks,
  #folds_column = "Fold",
  hexagon = F,
  rows_cols = c(2,2),
  progress = T,
  report = T,
  plot = T,
)

## Step 3: The analysis ##

### THIS ONE WORKS ###

## This version uses just mean canopy dieback and sitewide mean NDVI
# Initialize lists to store results
mae_results <- numeric(spatial_folds$k)
mse_results <- numeric(spatial_folds$k)
r_squared_results <- numeric(spatial_folds$k)

## This version of the model uses linear modelling
# Loop over each fold to train and test the model
for (fold in 1:spatial_folds$k) {
  # Get train and test indices for the current fold
  train_indices <- unlist(spatial_folds$folds_list[[fold]][1])  # First list element is train indices
  test_indices <- unlist(spatial_folds$folds_list[[fold]][2])   # Second list element is test indices
  
  # Subset data
  train_data <- SCV_Data_sf[train_indices, ]
  test_data <- SCV_Data_sf[test_indices, ]
  
  # Fit the model
  model <- lm(Total_Mean_Dieback_Percent ~ Mean_NDVI, data = train_data)
  
  # Predict on the test data
  predictions <- predict(model, newdata = test_data)
  
  # Calculate errors
  mae_results[fold] <- mean(abs(predictions - test_data$Total_Mean_Dieback_Percent))
  mse_results[fold] <- mean((predictions - test_data$Total_Mean_Dieback_Percent)^2)
  
  # Calculate R-squared
  sst <- sum((test_data$Total_Mean_Dieback_Percent - mean(train_data$Total_Mean_Dieback_Percent))^2)
  sse <- sum((predictions - test_data$Total_Mean_Dieback_Percent)^2)
  r_squared <- 1 - (sse / sst)
  r_squared_results[fold] <- r_squared
}

# Summarize cross-validation results
mean_mae <- mean(mae_results)
mean_mse <- mean(mse_results)
mean_r_squared <- mean(r_squared_results)

# Output the results
cat("Mean Absolute Error (MAE):", mean_mae, "\n")
cat("Mean Squared Error (MSE):", mean_mse, "\n")
cat("Mean R-squared:", mean_r_squared, "\n")

## This version uses mixed effects modelling (the relationship between total dieback percent & Juniper proportion)
# Initialize lists to store results
mae_results <- numeric(spatial_folds$k)
mse_results <- numeric(spatial_folds$k)
r_squared_results <- numeric(spatial_folds$k)

## This version of the model uses linear modelling
# Loop over each fold to train and test the model
for (fold in 1:spatial_folds$k) {
  # Get train and test indices for the current fold
  train_indices <- unlist(spatial_folds$folds_list[[fold]][1])  # First list element is train indices
  test_indices <- unlist(spatial_folds$folds_list[[fold]][2])   # Second list element is test indices
  
  # Subset data
  train_data <- SCV_Data_sf[train_indices, ]
  test_data <- SCV_Data_sf[test_indices, ]
  
  # Fit the model
  model <- lm(Total_Mean_Dieback_Percent ~ Mean_NDVI*Juniper_Site_Percent, data = train_data)
  
  # Predict on the test data
  predictions <- predict(model, newdata = test_data)
  
  # Calculate errors
  mae_results[fold] <- mean(abs(predictions - test_data$Total_Mean_Dieback_Percent))
  mse_results[fold] <- mean((predictions - test_data$Total_Mean_Dieback_Percent)^2)
  
  # Calculate R-squared
  sst <- sum((test_data$Total_Mean_Dieback_Percent - mean(train_data$Total_Mean_Dieback_Percent))^2)
  sse <- sum((predictions - test_data$Total_Mean_Dieback_Percent)^2)
  r_squared <- 1 - (sse / sst)
  r_squared_results[fold] <- r_squared
}

# Summarize cross-validation results
mean_mae <- mean(mae_results)
mean_mse <- mean(mse_results)
mean_r_squared <- mean(r_squared_results)

# Output the results
cat("Mean Absolute Error (MAE):", mean_mae, "\n")
cat("Mean Squared Error (MSE):", mean_mse, "\n")
cat("Mean R-squared:", mean_r_squared, "\n")


### VVV Testing scripts (not used) ####

## This version uses the deprecated "spatial_folds_old" object created by "spatialBlock"

# Initialize vectors to store results
SCV_Predicted_Dieback <- numeric(nrow(SCV_Data_sf)) # Empty vector for predictions
SCV_Actual_Dieback <- SCV_Data_sf$Total_Mean_Dieback_Percent # Actual data for comparison

# Loop through each spatial block for cross-validation
for (fold in 1:length(spatial_folds_old$foldID)) {
  
  # Extract training and test indices for the current spatial block
  tmp.test_indices <- unlist(spatial_folds_old$foldID[fold])
  tmp.SCV_training_data <- SCV_Data_sf[-tmp.test_indices, ]
  tmp.SCV_test_data <- SCV_Data_sf[tmp.test_indices, ]
  
  # Create a data frame for the training NDVI values
  training_NDVI <- numeric(nrow(tmp.SCV_training_data))
  
  # Iterate through the training data to find the nearest NDVI dates
  for (j in 1:nrow(tmp.SCV_training_data)) {
    tmp.site <- tmp.SCV_training_data$Site[j]
    if (grepl("_May", tmp.site)) {
      tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site, end=-5)][1]
    } else {
      tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site, end=-5)][1]
    }
    
    tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]]),]
    tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),]
    training_NDVI[j] <- Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  }
  
  # Do the same for the test data
  test_NDVI <- numeric(nrow(tmp.SCV_test_data))
  
  for (j in 1:nrow(tmp.SCV_test_data)) {
    tmp.site <- tmp.SCV_test_data$Site[j]
    if (grepl("_May", tmp.site)) {
      tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site, end=-5)][1]
    } else {
      tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site, end=-5)][1]
    }
    
    tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]]),]
    tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),]
    test_NDVI[j] <- Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  }
  
  # Create a training data frame
  tmp.SCV_training_df <- data.frame(Total_Mean_Dieback = tmp.SCV_training_data$Total_Mean_Dieback_Percent/100, Mean_NDVI = training_NDVI)
  
  # Train the model using training data
  tmp.SCV_model <- betareg(Total_Mean_Dieback ~ Mean_NDVI, data = tmp.SCV_training_df)
  
  # Predict canopy dieback for the test set
  tmp.SCV_predictions <- predict(tmp.SCV_model, newdata = data.frame(Mean_NDVI = test_NDVI))
  
  # Store predictions
  SCV_Predicted_Dieback[tmp.test_indices] <- tmp.SCV_predictions
}

# Assess the performance
SCV_MAE <- mean(abs(SCV_Predicted_Dieback - SCV_Actual_Dieback)) # Mean Absolute Error
SCV_MSE <- mean((SCV_Predicted_Dieback - SCV_Actual_Dieback)^2)  # Mean Squared Error
SCV_RMSE <- sqrt(SCV_MSE)                                        # Root Mean Squared Error
SCV_R_Squared <- cor(SCV_Predicted_Dieback, SCV_Actual_Dieback)^2 # R-squared

# Print performance metrics
cat("Mean Absolute Error (MAE):", SCV_MAE, "\n")
cat("Mean Squared Error (MSE):", SCV_MSE, "\n")
cat("Root Mean Squared Error (RMSE):", SCV_RMSE, "\n")
cat("R-squared:", SCV_R_Squared, "\n")

# Plot predicted vs actual dieback with labels
plot(SCV_Actual_Dieback, SCV_Predicted_Dieback, xlab = "Actual Canopy Dieback", ylab = "Predicted Canopy Dieback",
     xlim = c(0, 110), ylim = c(0, 110), main = "Spatial Cross-Validation: Actual vs Predicted", pch = 16, col = "blue")
abline(0, 1, col = "red")  # y=x reference line

# Add labels (assuming KBerg_Stand_Composition_May.Oct$Site contains site names)
text(SCV_Actual_Dieback, SCV_Predicted_Dieback, labels = SCV_Data_sf$Site, pos = 2, cex = 0.8, col = "black")

## This version uses the new "spatial_folds" object created by "cv_spatial"

# Initialize vectors to store results
SCV_Predicted_Dieback <- numeric(nrow(SCV_Data_sf)) # Empty vector for predictions
SCV_Actual_Dieback <- SCV_Data_sf$Total_Mean_Dieback_Percent # Actual data for comparison

# Loop through each spatial block for cross-validation
for (fold in 1:length(spatial_folds$folds_ids)) {
  
  # Extract training and test indices for the current spatial block
  tmp.test_indices <- unlist(spatial_folds$folds_ids[fold])
  tmp.SCV_training_data <- SCV_Data_sf[-tmp.test_indices, ]
  tmp.SCV_test_data <- SCV_Data_sf[tmp.test_indices, ]
  
  # Create a data frame for the training NDVI values
  training_NDVI <- numeric(nrow(tmp.SCV_training_data))
  
  # Iterate through the training data to find the nearest NDVI dates
  for (j in 1:nrow(tmp.SCV_training_data)) {
    tmp.site <- tmp.SCV_training_data$Site[j]
    if (grepl("_May", tmp.site)) {
      tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site, end=-5)][1]
    } else {
      tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site, end=-5)][1]
    }
    
    tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]]),]
    tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),]
    training_NDVI[j] <- Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  }
  
  # Do the same for the test data
  test_NDVI <- numeric(nrow(tmp.SCV_test_data))
  
  for (j in 1:nrow(tmp.SCV_test_data)) {
    tmp.site <- tmp.SCV_test_data$Site[j]
    if (grepl("_May", tmp.site)) {
      tmp.survey.date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == str_sub(tmp.site, end=-5)][1]
    } else {
      tmp.survey.date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == str_sub(tmp.site, end=-5)][1]
    }
    
    tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]]),]
    tmp.nearest_date <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey.date))),]
    test_NDVI[j] <- Mean_Site_NDVI[[paste0(str_sub(tmp.site, end=-5), "_NDVI")]][Mean_Site_NDVI$Strip_ID == tmp.nearest_date$Strip_ID]
  }
  
  # Create a training data frame
  tmp.SCV_training_df <- data.frame(Total_Mean_Dieback = tmp.SCV_training_data$Total_Mean_Dieback_Percent/100, Mean_NDVI = training_NDVI)
  
  # Train the model using training data
  tmp.SCV_model <- betareg(Total_Mean_Dieback ~ Mean_NDVI, data = tmp.SCV_training_df)
  
  # Predict canopy dieback for the test set
  tmp.SCV_predictions <- predict(tmp.SCV_model, newdata = data.frame(Mean_NDVI = test_NDVI))
  
  # Store predictions
  SCV_Predicted_Dieback[tmp.test_indices] <- tmp.SCV_predictions
}

# Assess the performance
SCV_MAE <- mean(abs(SCV_Predicted_Dieback - SCV_Actual_Dieback)) # Mean Absolute Error
SCV_MSE <- mean((SCV_Predicted_Dieback - SCV_Actual_Dieback)^2)  # Mean Squared Error
SCV_RMSE <- sqrt(SCV_MSE)                                        # Root Mean Squared Error
SCV_R_Squared <- cor(SCV_Predicted_Dieback, SCV_Actual_Dieback)^2 # R-squared

# Print performance metrics
cat("Mean Absolute Error (MAE):", SCV_MAE, "\n")
cat("Mean Squared Error (MSE):", SCV_MSE, "\n")
cat("Root Mean Squared Error (RMSE):", SCV_RMSE, "\n")
cat("R-squared:", SCV_R_Squared, "\n")

# Plot predicted vs actual dieback with labels
plot(SCV_Actual_Dieback, SCV_Predicted_Dieback, xlab = "Actual Canopy Dieback", ylab = "Predicted Canopy Dieback",
     xlim = c(0, 110), ylim = c(0, 110), main = "Spatial Cross-Validation: Actual vs Predicted", pch = 16, col = "blue")
abline(0, 1, col = "red")  # y=x reference line

# Add labels (assuming KBerg_Stand_Composition_May.Oct$Site contains site names)
text(SCV_Actual_Dieback, SCV_Predicted_Dieback, labels = SCV_Data_sf$Site, pos = 2, cex = 0.8, col = "black")

### ^^^ End of test scripts ####
### Q4: Can these sites predict other sites in the forest? (Lateral Analysis Approach 3: Using Spatial Random Forests) ####
### Spatial Random Forests is a package Dirk recommended which can find the most suitable model for ir/regular spatial data
## It should be quite simple to implement, although I have to also incorporate Blue, Green, and Red as separate predictor variables for the function to work

## Step 1: Create separate R, G, B, & NIR data columns

# Create an object to store the data
SpatialRF_Upscaling_Data <- KBerg_Stand_Composition_May.Oct

## Step 2: Find the R, G, & B values, plus the other indices I'll test, following the same method as used to get NDVI
# This is slightly inelegant code because I'm going for simplicity and repurposing existing methods

## Starting with redness values only

# Initialize a vector to store Red values
tmp.site_Red <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available Red date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_Red")]]),] # Filter for available Red values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Extract the mean Red value from that row
    tmp.site_Red[i] <- tmp.nearest_row[[paste0(tmp.site, "_Red")]]
  } else {
    tmp.site_Red[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_Red <- tmp.site_Red

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_Red[12] <- Mean_Site_NDVI$MD4_Red[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Now greenness values

# Initialize a vector to store Green values
tmp.site_Green <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available Green date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_Green")]]),] # Filter for available Green values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Extract the mean Green value from that row
    tmp.site_Green[i] <- tmp.nearest_row[[paste0(tmp.site, "_Green")]]
  } else {
    tmp.site_Green[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_Green <- tmp.site_Green

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_Green[12] <- Mean_Site_NDVI$MD4_Green[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## And thereafter, blueness values

# Initialize a vector to store Blue values
tmp.site_Blue <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available Blue date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_Blue")]]),] # Filter for available Blue values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Extract the mean Blue value from that row
    tmp.site_Blue[i] <- tmp.nearest_row[[paste0(tmp.site, "_Blue")]]
  } else {
    tmp.site_Blue[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_Blue <- tmp.site_Blue

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_Blue[12] <- Mean_Site_NDVI$MD4_Blue[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Now Near Infra-red values

# Initialize a vector to store Blue values
tmp.site_NIR <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available NIR date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_NIR")]]),] # Filter for available NIR values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Extract the mean Blue value from that row
    tmp.site_NIR[i] <- tmp.nearest_row[[paste0(tmp.site, "_NIR")]]
  } else {
    tmp.site_NIR[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_NIR <- tmp.site_NIR

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_NIR[12] <- Mean_Site_NDVI$MD4_NIR[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Now adding indices! Starting with EVI

# Initialize a vector to store Blue values
tmp.site_EVI <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_EVI")]]),] # Filter for available Blue values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Calculate the mean EVI value from that row
    tmp.site_EVI[i] <- tmp.nearest_row[[paste0(tmp.site, "_EVI")]]
  } else {
    tmp.site_EVI[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_EVI <- tmp.site_EVI

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_EVI[12] <- Mean_Site_NDVI$MD4_EVI[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Soil Adjusted Vegetation Index (SAVI)

# Initialize a vector to store Blue values
tmp.site_SAVI <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_SAVI")]]),] # Filter for available Blue values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Calculate the mean EVI value from that row
    tmp.site_SAVI[i] <- tmp.nearest_row[[paste0(tmp.site, "_SAVI")]]
  } else {
    tmp.site_SAVI[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_SAVI <- tmp.site_SAVI

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_SAVI[12] <- Mean_Site_NDVI$MD4_SAVI[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Green Normalised Difference Vegetation Index (GNDVI)

# Initialize a vector to store Blue values
tmp.site_GNDVI <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_GNDVI")]]),] # Filter for available Blue values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Calculate the mean GNDVI value from that row
    tmp.site_GNDVI[i] <- tmp.nearest_row[[paste0(tmp.site, "_GNDVI")]]
  } else {
    tmp.site_GNDVI[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_GNDVI <- tmp.site_GNDVI

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_GNDVI[12] <- Mean_Site_NDVI$MD4_GNDVI[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## The final index: Multiple Soil Adjusted Vegetation Indices 2 (MSAVI2)

# Initialize a vector to store Blue values
tmp.site_MSAVI2 <- numeric(nrow(SpatialRF_Upscaling_Data))

# Loop through each row of SpatialRF_Upscaling_Data
for (i in 1:nrow(SpatialRF_Upscaling_Data)) {
  # Extract the site name without "_May" or "_Oct"
  tmp.site <- sub("_.*", "", SpatialRF_Upscaling_Data$Site[i])
  
  # Determine if it's a May or October site and get the correct survey date
  if (grepl("_May", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_May$`Survey date`[KBerg_Mortality_Data_May$Site == tmp.site][1] 
  } else if (grepl("_Oct", SpatialRF_Upscaling_Data$Site[i])) {
    tmp.survey_date <- KBerg_Mortality_Data_Oct$`Survey date`[KBerg_Mortality_Data_Oct$Site == tmp.site][1] 
  } else {
    next  # Skip if site name doesn't match expected format
  }
  
  # Find the nearest available date for the site
  tmp.available_dates <- Mean_Site_NDVI[!is.na(Mean_Site_NDVI[[paste0(tmp.site, "_MSAVI2")]]),] # Filter for available Blue values
  
  if (nrow(tmp.available_dates) > 0) {
    # Get the row with the closest date
    tmp.nearest_row <- tmp.available_dates[which.min(abs(as.Date(tmp.available_dates$Date) - as.Date(tmp.survey_date))),]
    
    # Calculate the mean MSAVI2 value from that row
    tmp.site_MSAVI2[i] <- tmp.nearest_row[[paste0(tmp.site, "_MSAVI2")]]
  } else {
    tmp.site_MSAVI2[i] <- NA
  }
}

# Append the new column to SpatialRF_Upscaling_Data
SpatialRF_Upscaling_Data$Mean_MSAVI2 <- tmp.site_MSAVI2

## The site MD4_May value is anomalous, so I'll manually overwrite the selected date with the next nearest normal selected date's NDVI
# The next nearest MD4 instance is Strip No.2319101 from 2019/04/27

SpatialRF_Upscaling_Data$Mean_MSAVI2[12] <- Mean_Site_NDVI$MD4_MSAVI2[which(Mean_Site_NDVI$Strip_ID == 2319101)]

## Step 3: Add location data

# I will again use co-ordinate data from my KBerg AoI dataframe
# I'm doing this simply. If I change the site order, then the co-ordinates being appended will be mismatched
SpatialRF_Upscaling_Data <- bind_cols(SpatialRF_Upscaling_Data, c(KBerg_AoI$Latitude, KBerg_AoI$Latitude)) # add latitude c(repeated because I need 2*12 occurrences)
colnames(SpatialRF_Upscaling_Data)[23] <- "Latitude"
SpatialRF_Upscaling_Data <- bind_cols(SpatialRF_Upscaling_Data, c(KBerg_AoI$Longitude, KBerg_AoI$Longitude)) # add longitude c(repeated because I need 2*12 occurrences)
colnames(SpatialRF_Upscaling_Data)[24] <- "Longitude"

# Now that my data has latitude and longitude, convert it to a spatial object
sp::coordinates(SpatialRF_Upscaling_Data) <- ~Longitude + Latitude  # Pointing out the column names for coordinates
SpatialRF_Upscaling_Data <- st_as_sf(SpatialRF_Upscaling_Data, coords = c("Longitude", "Latitude"), crs = crs(KBerg_AoI)) # CRS was matched to the CRS of KBerg_AoI

## Step 4: Create a spatialRF model

# This model focuses on NDVI - with its constituent coloured wavelength reflectance values thrown in mainly just to make the code work
SpatialRF_Upscaling_Results <- spatialRF::rf_spatial(
  data = SpatialRF_Upscaling_Data,
  dependent.variable.name = "Total_Mean_Dieback_Percent",
  predictor.variable.names = c("Mean_NDVI", "Mean_Red", "Mean_Green", "Mean_Blue", "Mean_NIR"),
  distance.matrix = fields::rdist(st_coordinates(SpatialRF_Upscaling_Data$geometry)),
  #xy = "geometry",  
  verbose = TRUE
)

## Step 5: Initial results!
# Print results
print(SpatialRF_Upscaling_Results)

## Step 6: Identify the best model and use it to make predictions

SpatialRF.param.grid <- expand.grid(
  num.trees = c(100, 200, 500, 2500),  # Different numbers of trees
  mtry = c(2, 3, 4, 5),             # Number of predictors sampled at each split
  min.node.size = c(1, 5, 10, 15)    # Minimum number of observations in terminal nodes
)

# Create an empty list to store results
SpatialRF_search_results <- list()

# Loop over all combinations of hyperparameters
for (i in 1:nrow(SpatialRF.param.grid)) {
  cat("Running model", i, "of", nrow(SpatialRF.param.grid), "\n")
  
  # Extract current parameter values
  tmp.num.trees <- SpatialRF.param.grid$num.trees[i]
  tmp.mtry <- SpatialRF.param.grid$mtry[i]
  tmp.min.node.size <- SpatialRF.param.grid$min.node.size[i]
  
  # Train the spatial RF model
  tmp.SpatialRF_Upscaling_Model <- spatialRF::rf_spatial(
    data = SpatialRF_Upscaling_Data,
    dependent.variable.name = "Total_Mean_Dieback_Percent",
    #predictor.variable.names = c("Mean_NDVI", "Mean_Red", "Mean_Green", "Mean_Blue", "Mean_NIR"), # These are the simple, mainly NDVI predictors
    #predictor.variable.names = c("Mean_NDVI", "Mean_EVI", "Mean_SAVI", "Mean_GNDVI", "Mean_MSAVI2"), # These are a variety of indices, including classic NDVI
    predictor.variable.names = c("Mean_NDVI", "Mean_EVI", "Mean_SAVI", "Mean_GNDVI", "Mean_MSAVI2", "Mean_Red", "Mean_Green", "Mean_Blue", "Mean_NIR"), # These are all my available predictors
    distance.matrix = fields::rdist(st_coordinates(SpatialRF_Upscaling_Data$geometry)),
    ranger.arguments = list(
      num.trees = tmp.num.trees,
      mtry = tmp.mtry,
      min.node.size = tmp.min.node.size),
    verbose = T
  )
  
  # Store the model and its R? value
  SpatialRF_search_results[[i]] <- list(
    params = SpatialRF.param.grid[i, ],
    model = tmp.SpatialRF_Upscaling_Model,
    r_squared = tmp.SpatialRF_Upscaling_Model[["r.squared"]]
  )
}

# Find the best model based on R?^2
#best_model_index <- which.max(sapply(SpatialRF_search_results, function(x) x$r_squared))
#best_model <- grid_search_results[[best_model_index]]$model
SpatialRF_Best_Model <- SpatialRF_search_results[[which.max(sapply(SpatialRF_search_results, function(x) x$r_squared))]]$model

# Print the best model's parameters
cat("Best Model Parameters:\n")
print(SpatialRF_search_results[[which.max(sapply(SpatialRF_search_results, function(x) x$r_squared))]]$params)

# Step 6: Import the rest of the forest as a new AoI for defoliation predictions to be made

#####################################################################################
### DEFUNCT CODE SECTION vvv (Doesn't compute indices)
## Step 6.1: Load the GeoTIFF and Compute NDVI
# Load the multispectral raster data (PlanetScope GeoTIFF)
Utah_Forest_Composite <- rast("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Demo Data Small/composite.tif")

# Check which raster band is which (PlanetScope tends to use B, G, R, + NIR)
names(Utah_Forest_Composite)  # Ensure you know which bands correspond to Red, NIR, etc.

# Assuming Red = band 3, NIR = band 4 (adjust based on your dataset)
red_band <- Utah_Forest_Composite[[3]]
nir_band <- Utah_Forest_Composite[[4]]

# Compute NDVI: (NIR - Red) / (NIR + Red)
Forest.Composite.NDVI <- (nir_band - red_band) / (nir_band + red_band)

# Assign a meaningful name to the NDVI raster
names(Forest.Composite.NDVI) <- "NDVI"

# Save NDVI raster to file (optional)
writeRaster(Forest.Composite.NDVI , 
            "~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Demo Data Small/NDVI_of_Composite.tif", 
            overwrite = TRUE)

####################################################################################
### ACTIVE CODE SECTION VVV (computes indices as well)

## Step 6.1: Load the GeoTIFF and compute the indices
# Load the multispectral raster data (PlanetScope GeoTIFF)
Utah_Forest_Composite <- rast("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Utah Large Core Woodlands/composite_Core_PJ_Woodlands_Cropped.tif")

# Rename bands for clarity (B1 = Blue, B2 = Green, B3 = Red, B4 = NIR)
names(Utah_Forest_Composite) <- c("Blue", "Green", "Red", "NIR")

# Extract bands
Utah.Forest.Blue <- Utah_Forest_Composite[["Blue"]]
Utah.Forest.Green <- Utah_Forest_Composite[["Green"]]
Utah.Forest.Red <- Utah_Forest_Composite[["Red"]]
Utah.Forest.NIR <- Utah_Forest_Composite[["NIR"]]

# Calculate vegetation indices
Utah.Forest.NDVI   <- (Utah.Forest.NIR - Utah.Forest.Red) / (Utah.Forest.NIR + Utah.Forest.Red)
Utah.Forest.EVI    <- 2.5 * (Utah.Forest.NIR - Utah.Forest.Red) / (Utah.Forest.NIR + 6 * Utah.Forest.Red - 7.5 * Utah.Forest.Blue + 1)
Utah.Forest.SAVI   <- ((Utah.Forest.NIR - Utah.Forest.Red) * (1 + 0.5)) / (Utah.Forest.NIR + Utah.Forest.Red + 0.5)
Utah.Forest.GNDVI  <- (Utah.Forest.NIR - Utah.Forest.Green) / (Utah.Forest.NIR + Utah.Forest.Green)
Utah.Forest.MSAVI2 <- (2 * Utah.Forest.NIR + 1 - sqrt((2 * Utah.Forest.NIR + 1)^2 - 8 * (Utah.Forest.NIR - Utah.Forest.Red))) / 2

# Assign names
names(Utah.Forest.NDVI)   <- "NDVI"
names(Utah.Forest.EVI)    <- "EVI"
names(Utah.Forest.SAVI)   <- "SAVI"
names(Utah.Forest.GNDVI)  <- "GNDVI"
names(Utah.Forest.MSAVI2) <- "MSAVI2"

# Combine all the layers into a single raster stack
Utah_Forest_Reflectances.Indices <- 
  c(Utah.Forest.Blue, Utah.Forest.Green, Utah.Forest.Red, Utah.Forest.NIR, Utah.Forest.NDVI, Utah.Forest.EVI, Utah.Forest.SAVI, Utah.Forest.GNDVI, Utah.Forest.MSAVI2)

# Write to file (optional, each index separately or as a multi-layer raster)
writeRaster(Utah_Forest_Reflectances.Indices,
            "~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Utah Large Core Woodlands/Indices_of_composite_Core_PJ_Woodlands_Cropped.tif",
            overwrite = TRUE)
# For some reason, the raster doesn't actually save as part of the .rdata file, so let's load it back in
Utah_Forest_Reflectances.Indices <- rast("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Utah Large Core Woodlands/Indices_of_composite_Core_PJ_Woodlands_Cropped.tif")

####################################################################################

## Step 6.2: Define Pseudo Field Sites (706m? each)
# Get raster extent and projection
Utah.Forest.Extent <- ext(Utah_Forest_Reflectances.Indices)
Utah.Forest.CRS <- crs(Utah_Forest_Reflectances.Indices)

# Define grid cell size (approx. 706m? assume 26.5m x 26.5m per site)
grid_size <- sqrt(706)  # ~26.5m

# Generate a grid over the raster extent
grid <- st_make_grid(st_as_sfc(st_bbox(Utah.Forest.Extent, crs = Utah.Forest.CRS)), 
                     cellsize = grid_size, 
                     what = "polygons",
                     square = F,
                     flat_topped = F)

# Convert grid to sf object and assign IDs
pseudo_sites <- st_sf(site_id = seq_len(length(grid)), geometry = grid)

## Step 6.3: Identify and Remove/Tag Sites Affected by Clouds
# Load the Usable Data Mask raster (assuming cloud-affected pixels have a value of 0)
cloud_mask <- rast("~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Data Downloads + Review/PlanetScope/Utah Maps/Spatial Upscaling/Utah Large Core Woodlands/composite_udm2.tif")

# Extract cloud coverage for each pseudo site
#cloud_coverage <- exact_extract(cloud_mask, pseudo_sites, "mean") ### PROBLEM

#pseudo_sites$cloud_cover <- exact_extract(cloud_mask, pseudo_sites, function(df) mean(df$value, na.rm = TRUE)) ### PROBLEM

# Add cloud coverage information to pseudo_sites
pseudo_sites$cloud_cover <- cloud_coverage

# Filter out sites with significant cloud coverage (e.g., if > 0.01% of pixels are cloudy)
pseudo_sites_clean <- pseudo_sites %>% filter(cloud_cover < 0.001)

## Step 6.4: Compute Mean Values of NDVI, Red, Green, Blue, and NIR per Site

# Extract mean values for each spectral band (assuming PlanetScope bands)
#Utah_Forest_Site_Means <- pseudo_sites_clean %>%
#  mutate(
#    Mean_NDVI = exact_extract(Forest.Composite.NDVI, pseudo_sites_clean, fun = mean, na.rm = TRUE, summarize_df = TRUE),
#    Mean_Red = exact_extract(red_band, pseudo_sites_clean, fun = mean, na.rm = TRUE, summarize_df = TRUE),
#    Mean_Green = exact_extract(Utah_Forest_Composite[[2]], pseudo_sites_clean, fun = mean, na.rm = TRUE, summarize_df = TRUE), # Assuming band 2 is green
#    Mean_Blue = exact_extract(Utah_Forest_Composite[[1]], pseudo_sites_clean, fun = mean, na.rm = TRUE, summarize_df = TRUE),  # Assuming band 1 is blue
#    Mean_NIR = exact_extract(nir_band, pseudo_sites_clean, fun = mean, na.rm = TRUE, summarize_df = TRUE)
#  )

Pinon.Juniper_SpatialUpscaling_Sites <- pseudo_sites %>%
  mutate(
    Mean_NDVI = exact_extract(Utah_Forest_Reflectances.Indices$NDVI, pseudo_sites, "mean"),
    Mean_GNDVI = exact_extract(Utah_Forest_Reflectances.Indices$GNDVI, pseudo_sites, "mean"),
    Mean_EVI = exact_extract(Utah_Forest_Reflectances.Indices$EVI, pseudo_sites, "mean"),
    Mean_MSAVI2 = exact_extract(Utah_Forest_Reflectances.Indices$MSAVI2, pseudo_sites, "mean"),
    Mean_SAVI = exact_extract(Utah_Forest_Reflectances.Indices$SAVI, pseudo_sites, "mean"),
    Mean_Red = exact_extract(Utah_Forest_Reflectances.Indices$Red, pseudo_sites, "mean"), # Wouldn't it be cool if instead of just R,G,B, & NDVI values, I had also GNDVI and other indices?
    Mean_Green = exact_extract(Utah_Forest_Reflectances.Indices$Green, pseudo_sites, "mean"), # Assuming band 2 is green
    Mean_Blue = exact_extract(Utah_Forest_Reflectances.Indices$Blue, pseudo_sites, "mean"),  # Assuming band 1 is blue
    Mean_NIR = exact_extract(Utah_Forest_Reflectances.Indices$NIR, pseudo_sites, "mean")
  )

# Convert to an sf object for compatibility with SpatialRF
# Pinon_Juniper_SpatialUpscaling_Sites <- st_as_sf(Utah_Forest_Site_Means) # THIS DOES NOTHING, both objects are exactly the same

## Step 6.5: Clean up the pseudo sites' data

# Save the prepared dataset
saveRDS(Pinon.Juniper_SpatialUpscaling_Sites, "Pinon_Juniper_SpatialUpscaling_Sites_V1.3.rds")

# Use the best model to make predictions
#SpatialRF_Upscaling_Predictions <- predict(SpatialRF_Best_Model, newdata = SpatialRF_Prediction_Data)  # Replace with your new sites data

Pinon.Juniper_SpatialUpscalingSites_Clean <- Pinon.Juniper_SpatialUpscaling_Sites |>
  dplyr::filter(
    !is.na(Mean_NDVI),
    !is.na(Mean_GNDVI),
    !is.na(Mean_EVI),
    !is.na(Mean_MSAVI2),
    !is.na(Mean_SAVI),
    !is.na(Mean_Red),
    !is.na(Mean_Green),
    !is.na(Mean_Blue),
    !is.na(Mean_NIR)
  )

SpatialUpscaling_Prediction_Matrix <- Pinon.Juniper_SpatialUpscalingSites_Clean |>
  dplyr::select(all_of(c("Mean_NDVI", "Mean_GNDVI", "Mean_EVI", "Mean_MSAVI2", "Mean_SAVI", "Mean_Red", "Mean_Green", "Mean_Blue", "Mean_NIR"))) |>
  as.data.frame()

## Step 7: Make predictions using the pseudo sites and best RF model

Pinon.Juniper_SpatialUpscaling_Results <- predict(
  object = SpatialRF_Best_Model,  # the best ranger model
  data = SpatialUpscaling_Prediction_Matrix,  # new data for prediction
  type = "response" # the default argument
)

# This can take forever to run, so it makes sense to save the output
saveRDS(Pinon.Juniper_SpatialUpscaling_Results, "~/Work work/Philipps Universit?t Marburg/Defoliation Timing and Tree Mortality/Code/Pinon_Juniper_SpatialUpscaling_Results_V1.3.rds")
Pinon.Juniper_SpatialUpscaling_Results <- readRDS("~/Work work/Philipps Universität Marburg/Defoliation Timing and Tree Mortality/Code/Pinon_Juniper_SpatialUpscaling_Results_V1.3.rds")

## Step 8: Export these data so I can make a map in QGIS

Pinon.Juniper_SpatialUpscalingSites_Clean$Predicted_Dieback <- Pinon.Juniper_SpatialUpscaling_Results[["predictions"]]

#pseudo_sites <- pseudo_sites |> 
#  dplyr::mutate(site_id = row_number()) |>  # Ensure a unique identifier
#  dplyr::left_join(
#    data.frame(site_id = 1:length(Pinon.Juniper_SpatialUpscaling_Results$predictions), 
#               predicted_dieback = Pinon.Juniper_SpatialUpscaling_Results$predictions),
#    by = "site_id"
#  )
st_write(Pinon.Juniper_SpatialUpscalingSites_Clean, "Pinon_Juniper_Spatial_Upscaling_Dieback_Predictions_V1.3.gpkg", append = F)

## Step 9: Make some other things!

## Distribution histograms
# Here's a basic histogram:
hist(Pinon.Juniper_SpatialUpscaling_Results[["predictions"]])

# Here's a more publication ready histogram (Figure 6):
ggplot(mapping = aes(x = Pinon.Juniper_SpatialUpscaling_Results[["predictions"]])) +
  geom_histogram(binwidth = 2.5, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(
    title = "Distribution of Predicted Sitewide Mean Relative Canopy Dieback",
    x = "Predicted Relative Canopy Dieback (%)",
    y = "Frequency of Pseudo-sites in Prediction Range"
  ) +
  scale_y_continuous(labels = scales::label_comma()) +  # Y axis in plain numbers
  scale_x_continuous(breaks = seq(0, 55, by = 5)) +  # X axis in 5% steps
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11))

# Now save it (MDPI proportions)!
ggsave("Figure6.png", width = 18.46, height = 11, units = "cm", dpi = 600)

## Mapping the predicted relative canopy dieback:
### SECTION UNUSED: this crashed my PC, so I did it in QGIS instead ###
# ggplot(Pinon.Juniper_SpatialUpscalingSites_Clean) +
#   geom_sf(aes(fill = Predicted_Dieback), color = NA) +  # Remove polygon borders for speed
#   scale_fill_viridis_c(
#     name = "Predicted Mean of Sitewide Canopy Dieback (%)",
#     option = "plasma",
#     direction = -1
#   ) +
#   labs(
#     title = "Spatial Upscaling of Predicted Relative Canopy Dieback",
#     subtitle = "Across Pi?on-Juniper Forests in Southeastern Utah"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5))
