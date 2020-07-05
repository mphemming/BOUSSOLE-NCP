%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% Compare_wind.m
%------------
% Compare wind reanalysis with BOUSSOLE meteorological buoy
%-----------------------------------------------------------------------------------------------------
% script created by MPH in Sydney, 05/07/2020
%

clear all
close all
clc

%% load data


options.directory = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts'
options.data_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/data/'
options.plot_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/Plots/'
addpath(genpath(options.directory));
load([options.data_dir,'prcdata.mat'],'prcdata');
load([options.data_dir,'METEO.mat']);  
load([options.data_dir,'BOUSSOLE.mat']);  
load([options.data_dir,'closeness.mat']);  
load([options.data_dir,'CTD.mat']);  
% reanalysis data
CCMP_201603 = get_mooring('CCMP_Wind_Analysis_201603_V02.0_L3.5_RSS.nc',0);
CCMP_201604 = get_mooring('CCMP_Wind_Analysis_201604_V02.0_L3.5_RSS.nc',0);
era5_201603 = get_mooring('Wind_WFDE5_CRU_201603_v1.0.nc',0);
era5_201604 = get_mooring('Wind_WFDE5_CRU_201604_v1.0.nc',0);

%% select data

[era5_201603.lon_mat era5_201603.lat_mat] = meshgrid(era5_201603.lon, era5_201603.lat);
era5_201603.time_mat = repmat(era5_201603.time,[1 360])';
check = era5_201603.lon_mat == 7.75 & era5_201603.lat_mat == 43.75;
hr = 0;
dy = 0;
for ts = 1:744
    hr = hr+1;
    if hr == 25
        hr = 1;
        dy = dy+1
    end
    era5_201603.selection.time(ts) = datenum(2016,03,dy,hr,0,0);
    era5_201603.selection.lon(ts) = era5_201603.lon_mat(check);    
    era5_201603.selection.lat(ts) = era5_201603.lat_mat(check); 
    wind_now = squeeze(era5_201603.Wind(:,:,ts));
    era5_201603.selection.wind(ts) = wind_now(check');
end

%% create comparison plot






