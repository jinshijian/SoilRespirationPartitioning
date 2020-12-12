library('raster')
library('ggplot2')
library('gridExtra')
library('grDevices')
library('RColorBrewer')
library("LSD")
library("rgdal")
library("maptools")
library("png") 
library("splancs")
library("plotrix")
library("phonR")
# install.packages('spatialEco')
library(spatialEco)
library(latticeExtra)
library(dplyr)


#### zonal mean of Ra for predicted and hashimoto Ra **********
### predicted Ra

# Bulk Density, EVI, and MAP

#Climate_type, MAP, MAT, Percentage_AM, Percentage_EcM, Percentage_ErM, Percentage_NM, 
#BMa, BMb, N_deposition, BD, Clay, SOC, EVI, RC, SD_RC, CV_RC

EVI <- raster("outputs/RC_factor_tif/EVI.tif")
BD <- raster("outputs/RC_factor_tif/BD.tif")
MAP <- raster("outputs/RC_factor_tif/MAP.tif")
RC <- raster("outputs/RC_factor_tif/RC.tif")

cor.BD <- rasterCorrelation(RC, BD, s = 3, na.rm = T, type = "pearson","cor.BD.tif",overwrite = T)
cor.EVI <- rasterCorrelation(RC, EVI, s = 3, na.rm = T, type = "pearson","cor.EVI.tif",overwrite = T)
cor.MAP <- rasterCorrelation(RC, MAP, s = 3, na.rm = T, type = "pearson","cor.MAP.tif",overwrite = T)

#### plot results out *****************************************
# Coordinates of the triangle
tri <- rbind(sin(0:2*2/3*pi), cos(0:2*2/3*pi))

# Function for calculating the color of a set of points `pt`
# in relation to the triangle
tricol <- function(pt, sharpness=1.6){
  #print("inner function")
  #print(pt)
  RGB <- sapply(1:3, function(i){
    a <- sweep(pt, 2, tri[,i])
    
    b <- apply(tri[,-i], 1, mean) - tri[,i]
    sharpness*((a %*% b) / sum(b^2))-sharpness+1
  })
  RGB[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
  
  do.call(rgb, unname(as.data.frame(pmin(pmax(RGB, 0), 1))))
}

# Plot


res <- 1000                       # Resolution
xi <- seq(-0.9, 0.9, length=res)        # Axis points
yi <- seq(-.8, 1, length=res)
x <- xi[1] + cumsum(diff(xi))       # Midpoints between axis points
y <- yi[1] + cumsum(diff(yi))
xy <- matrix(1:(length(x)*length(y)), length(x))

png(filename="outputs/dom_factor/legend.png",width = 600, height = 450)
op <- par(mfrow = c(1,1),mar=c(1,1,1,1)+0.1)
#png("test.png")
#image(xi, yi, xy, xlim=c(-0.9,0.9),ylim=c(-0.8,1.0),col=tricol(as.matrix(expand.grid(x,y))), useRaster=TRUE,xlab="",ylab="",axes=F)

image(xi, yi, xy, xlim=c(-1.2, 1.2),ylim=c(-1.2,1.5),col=tricol(as.matrix(expand.grid(x,y))), 
      useRaster=TRUE,xlab="",ylab="",axes=F)

lines(tri[1,c(1:3,1)], tri[2,c(1:3,1)], type="l")
#########

### red for MAP, green for BD, blue for EVI

# ,srt=-90  set the angle of the text

text( x = 0, y = -0.8, label = "MAP",cex=7) 
text( x = -.8, y = 0.4, label = "BD",cex=7)
text( x = .8, y = 0.4, label = "EVI",cex=7)

dev.off()


require(png)
rgb.legend <- readPNG("outputs/dom_factor/legend.png")


pdf("outputs/dom_factor/dominant.pdf",width = 8,height=5)


Cor.files <- list.files(pattern = "cor")

Cor.tif <-lapply(Cor.files,raster)

### red for MAP, green for BD, blue for EVI

rgbRaster.d <- stack(Cor.tif[[3]],Cor.tif[[1]],Cor.tif[[2]])

par(col.axis = "white", col.lab = "white", tck = 0)

plotRGB(rgbRaster.d, stretch = 'lin',axes = T, maxpixels = 100000)

box(col = "white")

rasterImage(rgb.legend,xleft = -170, ybottom = -65, xright = -80, ytop = 10)

dev.off()


