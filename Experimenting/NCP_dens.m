%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% NCP_dens.m
%------------
% Calculate NCP using a mass balance, including different physical terms
%-----------------------------------------------------------------------------------------------------
% script by MPH in Norwich, 14/11/2017
% script modified by MPH in Sydney, 21/06/2019
%

clear all
close all
clc

%% Addpath and load in data

options.directory = '/Users/Michael/Documents/Work/UEA/NCP_scripts'
addpath(genpath(options.directory));
load([options.directory,'/data/prcdata.mat'],'prcdata');
load([options.directory,'/data/METEO.mat']);  

%% define options for NCP calculation

% top layer of interest 
options.h = 1029; % kg m-3 isopycnal
% moving window range 
options.window = 2; % (e.g. +- number)
% specify timesteps
options.dayrange = [9:1:35]; % day range used
options.interval = 1; % need this for some scripts, day interval

disp(['Layer density = ',num2str(options.h), ' kg m-3'])
disp(['Moving window = +- ',num2str(options.window),' days'])

%% Organise variable arrays for NCP-related calculations

% Required parameters
vars.P = prcdata.timeseries.P;
vars.depth = prcdata.timeseries.depth;
vars.S = prcdata.timeseries.S;
vars.absS = prcdata.timeseries.absS;
vars.T = prcdata.timeseries.T;
vars.consT = prcdata.timeseries.T;
vars.t = prcdata.timeseries.t;
vars.lon = prcdata.timeseries.lon;
vars.lat = prcdata.timeseries.lat;
vars.sigma0 = prcdata.timeseries.sigma0;
vars.MLD = prcdata.timeseries.MLDO2; % MLD from oxygen concentration
vars.MLDO2 = prcdata.timeseries.MLDO2val; % oxygen concentration at MLD 
vars.oxycline = prcdata.timeseries.oxyclineD; % depth of oxycline
vars.oxyclineO2 = prcdata.timeseries.oxyclineDval; % oxygen concentration at oxycline
vars.O2 =  (prcdata.timeseries.sigma0/1000) .* prcdata.timeseries.O2; % calibrated O2, convert µmol kg -> mmol m-3
vars.Fl = prcdata.timeseries.Chl_Fl;
vars.Sc700 = prcdata.timeseries.Scatter_700;
vars.dive = prcdata.timeseries.dive_vector;
vars.DIC = ([prcdata.CO2SYS.DIC]./vars.S) * 38.3; % normalised DIC for S
vars.DIC = (prcdata.timeseries.sigma0/1000) .* vars.DIC; % convert µmol kg -> mmol m-3
vars.fCO213 = ([prcdata.CO2SYS.fCO2]).*exp(0.0423*(13-vars.T)); % normalised FCO2 for T
vars.fCO2 = [prcdata.CO2SYS.fCO2];
% add wind from meteo buoy,
% needs selecting and interpolating onto same time as glider data
wind_selection = METEO.time_wind_speed >= vars.t(1) & ...
    METEO.time_wind_speed <= vars.t(end);
l1 = linspace(1,length(METEO.time_wind_speed(wind_selection)),length(vars.t));
vars.Wind_time = interp1(METEO.time_wind_speed(wind_selection)',l1,'Linear');
vars.Wind = interp1(METEO.wind_speed(wind_selection)',l1,'Linear');
% same but with atmospheric pressure
press_selection = METEO.date >= vars.t(1) & METEO.date <= vars.t(end);
l1 = linspace(1,length(METEO.date(press_selection)),length(vars.t));
vars.sea_level_press_time = interp1(METEO.date(press_selection)',l1,'Linear');
vars.sea_level_press = interp1(METEO.sea_level_pressure(press_selection)',l1,'Linear');
% DACs
vars.DACs.DACu = [prcdata.hydrography.DAC_u];
vars.DACs.DACv = [prcdata.hydrography.DAC_v];
% DAC QC
vars.DACs.DACu(vars.DACs.DACu > 0.2) = NaN;
vars.DACs.DACv(vars.DACs.DACv > 0.2 | vars.DACs.DACv < -0.2) = NaN;
% get average/median parameters for each dive to match DACs in depth of
% layer (h) specified at start
for dive = 1:147
    check_dive = vars.sigma0 <= options.h & vars.dive == dive;
    vars.DACs.t(dive) = nanmedian([prcdata.hydrography(dive).time]);
    vars.DACs.lon(dive) = nanmedian([prcdata.hydrography(dive).lon]);
    vars.DACs.lat(dive) = nanmedian([prcdata.hydrography(dive).lat]);
    vars.DACs.absS_mean_50m(dive) = nanmean(vars.absS(check_dive));
    vars.DACs.consT_mean_50m(dive) = nanmean(vars.consT(check_dive));
    vars.DACs.O2_mean_50m(dive) = nanmean(vars.O2(check_dive));
    vars.DACs.wind(dive) = nanmean(vars.Wind(vars.dive == dive));
    vars.DACs.MLD(dive) = nanmean(vars.MLD(vars.dive == dive));
end
% satellite measurements
vars.GVsat_time =  [prcdata.GVsat.day];
vars.GVsat_U =  [prcdata.GVsat.abs_U_daymean_Glider];
vars.GVsat_V =  [prcdata.GVsat.abs_V_daymean_Glider];

clear l1 ans wind_selection check_dive dive 

%% run individual mass balance modules / terms
% calculate geopotential anomaly profiles
run NCP_dens_GPA.m;
% obtain plane-fits of geopotential anomalies and oxygen
run NCP_dens_planes.m
% get oxygen inventory change with time
run NCP_dens_inventory.m
% calculate entrainment
run NCP_entrainment.m
% calculate advection term
run NCP_advection.m
% calculate air-sea exchange
run NCP_airseaexchange.m

%% calculate NCP and associated errors
% concatenate advection
ADV = [O2_adv.adv];
NCP_est = O2_inv.inv(2:end-1)'  + ADV(2:end-1)  + O2_ase.ASE(2:end-1) - [O2_ent.ent];
% get errors
run NCP_error.m
