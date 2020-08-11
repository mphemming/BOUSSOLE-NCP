%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% NCP.m
%------------
% Calculate NCP using a mass balance, including different physical terms
%-----------------------------------------------------------------------------------------------------
% script by MPH in Norwich, 14/11/2017
% script modified by MPH in Sydney, 21/06/2020
%

clear all
close all
clc

%% Addpath and load in data

options.directory = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts'
options.data_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\data\'
options.plot_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\Plots\'
addpath(genpath(options.directory));
load([options.directory,'\data\prcdata.mat'],'prcdata');
load([options.directory,'\data\METEO.mat']);  

%% define options for NCP calculation

% top layer of interest 
options.h = 46; % metres
% moving window range 
options.window = 2; % (e.g. +- number)
% moving window range extender for advection planes
options.window_adv = 2; %e.g. total would be this multiplied by options.window  above
% specify timesteps
options.dayrange = [9:1:35]; % day range used
options.interval = 1; % need this for some scripts, day interval

disp(['Layer depth = ',num2str(options.h), ' metres'])
disp(['Moving window = +- ',num2str(options.window),' days'])

%% Organise variable arrays for NCP-related calculations

% Required parameters
vars.dive = prcdata.timeseries.dive_vector;
vars.P = prcdata.timeseries.P;
vars.depth = prcdata.timeseries.depth;
vars.S = prcdata.timeseries.S;
vars.absS = prcdata.timeseries.absS;
vars.T = prcdata.timeseries.T;
vars.consT = prcdata.timeseries.T;
vars.Chl = prcdata.timeseries.Chl_Fl;
vars.t = prcdata.timeseries.t;
vars.lon = prcdata.timeseries.lon;
vars.lat = prcdata.timeseries.lat;
vars.sigma0 = prcdata.timeseries.sigma0;
vars.MLD = prcdata.timeseries.MLDO2; % MLD from oxygen concentration
vars.MLDO2 = prcdata.timeseries.MLDO2val; % oxygen concentration at MLD 
vars.oxycline = prcdata.timeseries.oxyclineD; % depth of oxycline
vars.oxyclineO2 = prcdata.timeseries.oxyclineDval; % oxygen concentration at oxycline
%%%%%%%%%%%%%%%
vars.O2 =  (prcdata.timeseries.sigma0/1000) .* prcdata.timeseries.O2; % calibrated O2, convert µmol kg -> mmol m-3
% Some further QC for O2
vars.O2(vars.depth < 5 | vars.O2 > 300) = NaN; % remove top 5m where spikes
vars.O2_MLD = vars.O2;
for n_dives = 1:numel(prcdata.hydrography)
    % remove dives less than 20 m deep
    if prcdata.hydrography(n_dives).max_pressure < 20
        vars.O2(prcdata.timeseries.dive_vector == n_dives) = NaN; 
        vars.O2_MLD(prcdata.timeseries.dive_vector == n_dives) = NaN;         
    end
    % O2 only in MLD
    MLD = nanmean(vars.MLD(prcdata.timeseries.dive_vector == n_dives));
    vars.O2_MLD(prcdata.timeseries.dive_vector == n_dives & vars.depth > MLD) = NaN;
    % determine mean gradients in euphotic depth
     check_top = vars.depth <= 10 & vars.dive == n_dives;
     check_bot = vars.depth >= 40 & vars.depth <= 50 & vars.dive == n_dives;
     check_up = prcdata.timeseries.downup_vector == 2;
     grad_diff_up(n_dives) = nanmedian(vars.O2(check_top & check_up))-nanmedian(vars.O2(check_bot & check_up));
     grad_diff_down(n_dives) = nanmedian(vars.O2(check_top & ~check_up))-nanmedian(vars.O2(check_bot & ~check_up));
end
vars.O2_sat = o2satSTP(vars.T, vars.S, 1013); % should use changing pressure!
vars.O2_sat = (prcdata.timeseries.sigma0/1000) .* vars.O2_sat;
% use O2 profiles with low gradient for advection
for n_dives = 1:numel(prcdata.hydrography)
    check_down =  vars.dive == n_dives & prcdata.timeseries.downup_vector == 1;
    check_up =  vars.dive == n_dives & prcdata.timeseries.downup_vector == 2;
    % down
    if abs(grad_diff_down(n_dives)) < 17
        vars.O2_adv(check_down) =  vars.O2(check_down);
    else
        vars.O2_adv(check_down) =  ones(size(vars.O2(check_down)))*NaN;
    end
    % up
    if abs(grad_diff_down(n_dives)) < 17
        vars.O2_adv(check_up) =  vars.O2(check_up);
    else
        vars.O2_adv(check_up) =  ones(size(vars.O2(check_up)))*NaN;
    end    
end
vars.O2_adv = vars.O2; % this method not used
%%%%%%%%%%%%%%%
vars.compensation = NaN(size(vars.O2));
vars.compensation(vars.O2-vars.O2_sat > -1 & vars.O2-vars.O2_sat < 1) = vars.depth(vars.O2-vars.O2_sat > -1 & vars.O2-vars.O2_sat < 1);
vars.Fl = prcdata.timeseries.Chl_Fl;
vars.Sc700 = prcdata.timeseries.Scatter_700;
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
    check_dive = vars.P <= options.h & vars.dive == dive;
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

clear l1 ans wind_selection check* dive 

%% run individual mass balance modules / terms
% calculate geopotential anomaly profiles
run NCP_GPA.m
% obtain plane-fits of geopotential anomalies and oxygen
run NCP_planes_depthres.m
run NCP_planes_depthres_DIC.m
% run /Users/Michael/Documents/Work/UEA/NCP_Scripts/NCP_planes.m
% run /Users/Michael/Documents/Work/UEA/NCP_Scripts/NCP_planes_DIC.m
% get oxygen inventory change with time
run NCP_inventory.m
run NCP_inventory_DIC.m
% calculate entrainment
run NCP_entrainment.m
run NCP_entrainment_DIC.m
% calculate advection term
run NCP_advection_depthres.m
% run NCP_advection.m
run NCP_advection_depthres_DIC.m
% calculate air-sea exchange
run NCP_airseaexchange.m
run NCP_airseaexchange_DIC.m
% get diapycnal diffusion
run NCP_kz.m % exclude this, not required
%% calculate NCP and associated errors
% concatenate advection
ADV = [O2_adv.adv]; ADV_DIC = [DIC_adv.adv];
NCP_est = O2_inv.inv(2:end-1)'  + ADV(2:end-1)  + O2_ase.ASE(2:end-1) - [O2_ent.ent];
NCP_est_DIC = DIC_inv.inv(2:end-1)'  + ADV_DIC(2:end-1)  + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent];
NCP_est_kz = O2_inv.inv(2:end-1)'  + ADV(2:end-1)  + O2_ase.ASE(2:end-1) - [O2_ent.ent] - [kz(2:end-1).kz];
NCP_est_kz_no_adv = O2_inv.inv(2:end-1)' + O2_ase.ASE(2:end-1) - [O2_ent.ent] - [kz(2:end-1).kz];
NCP_est_kz_DIC = DIC_inv.inv(2:end-1)'  + ADV_DIC(2:end-1)  + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent] - [kz(2:end-1).kz_DIC];
NCP_est_kz_no_adv_DIC = DIC_inv.inv(2:end-1)' + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent] - [kz(2:end-1).kz_DIC];
NCP_est_no_adv = O2_inv.inv(2:end-1)' + O2_ase.ASE(2:end-1) - [O2_ent.ent];
NCP_est_no_adv_DIC = DIC_inv.inv(2:end-1)' + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent];
% get errors
% run NCP_error.m 
plot_days = options.dayrange(2:end-1);

%% run NCP_buoy and create table of vals
% get buoy vals
run NCP_buoy
% get table vales
run NCP_table.m

% NOTE: no errors added at the moment


%% save NCP
save([options.data_dir,'NCP_',num2str(options.h),'m.mat']);  
