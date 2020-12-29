clc;
clear all;

%% MODIS Land Cocer
landcovers = imread('MODIS_LC_2010.tif');
fun = @(block_struct) mode(block_struct.data(:));
upscale_factor = [5 5];
landcovers = blockproc(landcovers, upscale_factor, fun);
masks = nan(size(landcovers));
masks((landcovers >0 & landcovers <13) | (landcovers == 14) ) = 1;

%% landmask
load('landmasks.mat');

BDs = double(imread('cor.BD.tif'));
EVIs = double(imread('cor.EVI.tif'));
MAPs = double(imread('cor.MAP.tif'));

%% Remove abnormal values and normalize to 0-1
BDs(BDs<-1 | BDs>1) = 1;
EVIs(EVIs<-1 | EVIs>1) = 1;
MAPs(MAPs<-1 | MAPs>1) = 1;

inmin = -1;
inmax = 1;
RGBimage = cat(3, rescale(BDs,'InputMin',inmin,'InputMax',inmax),...
    rescale(EVIs,'InputMin',inmin,'InputMax',inmax),...
    rescale(MAPs,'InputMin',inmin,'InputMax',inmax));

RGBimage = im2uint8(RGBimage);
%% lat lon
res_v = 0.5;
res_h = 0.5;
lon = (-180+res_h/2):res_h: (180-res_h/2);
lat = (90-res_v/2):-res_v: (-90+res_v/2);
[lons,lats]=meshgrid(lon,lat);

%% plot
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.38,0.3]);
set(gca, 'Position', [0.05 0 0.92 1])
plot_global_map(lat, lon, RGBimage);

%% save
print(gcf, '-dtiff', '-r600', 'dominant_effect.tif')
