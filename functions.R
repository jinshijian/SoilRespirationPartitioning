
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)
writ_file <- function(input, output) write.csv(input, file.path(OUT_DIR, output), row.names = FALSE)
'%!in%' <- function(x,y)!('%in%'(x,y))

#*****************************************************************************************************************
# functions for analysis
#*****************************************************************************************************************
get_evi <- function(sdata){
  evi_mean <- raster("Data/EVI_mean.tif")
  for (i in seq_len(nrow(sdata))){
    lat <- sdata$Latitude[i]
    lon <- sdata$Longitude[i]
    if (is.na(lat) | is.na(lon)) { next } 
    else {
      lat_low <- lat - 0.1
      lat_high <- lat + 0.1
      lon_low <- lon - 0.1
      lon_high <- lon + 0.1
      
      ex_cell <- raster::extent(c(lon_low, lon_high, lat_low, lat_high))
      evi_cell <- raster::extract(evi_mean, ex_cell, fun = mean, na.rm = TRUE)
      
      sdata$EVI_2[i] <- evi_cell
    }
    print(paste0("*****", i))
  }
  return(sdata)
}
