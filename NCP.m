%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% NCP.m
%------------
% Calculate NCP using a budget
%-----------------------------------------------------------------------------------------------------
% script by MPH in Norwich, 14/11/2017
% script modified by MPH in Sydney, 05/10/2020
%

% clear all
% close all
clc

%% Addpath and load in data

options.directory = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts'
options.data_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\data\'
options.plot_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\Plots\'
addpath(genpath(options.directory));
load([options.directory,'\data\prcdata.mat'],'prcdata');
load([options.directory,'\data\METEO.mat']);  

%% define options for NCP calculation

% layer of interest 
options.h = 46; % metres
% moving window range 
options.window = 2; % (e.g. +- number)
% moving window range extender for advection planes
options.window_adv = 2; %e.g. total would be this multiplied by options.window  above
% specify timesteps
options.dayrange = [9:1:35]; % day range used
options.interval = 1; % need this for some scripts, day interval
options.vertical_grid = 1:5:46; % grid used for GPA and O2 gradients

disp(['Layer depth = ',num2str(options.h), ' metres'])
disp(['Moving window = +- ',num2str(options.window),' days'])
disp(['ADV moving window = +- ',num2str(options.window_adv*options.window),' days'])
%% Organise variable arrays for NCP-related calculations

run NCP_prepare_data.m
    
%% run individual mass balance modules / terms
% calculate geopotential anomaly profiles
run  'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_GPA.m'
% obtain plane-fits of geopotential anomalies and oxygen
run  'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_plane_fits.m'
% get oxygen inventory change with time
run  'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_inv.m'
% calculate entrainment
run 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_ent.m'
% calculate advection term
run 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_adv.m'
% calculate air-sea exchange
run 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_airsea.m'
% get diapycnal diffusion
run 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_kz.m' 
%% calculate NCP and associated errors
% concatenate advection
ADV = [O2_adv.adv]; ADV_DIC = [DIC_adv.adv];
NCP_est = O2_inv.inv(2:end-1)'  + ADV(2:end-1)  + O2_ase.ASE(2:end-1) - [O2_ent.ent];
NCP_est_DIC = DIC_inv.inv(2:end-1)'  + ADV_DIC(2:end-1)  + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent];
NCP_est_kz = O2_inv.inv(2:end-1)'  + ADV(2:end-1)  + O2_ase.ASE(2:end-1) - [O2_ent.ent] -  kz_(2:end-1);
NCP_est_kz_no_adv = O2_inv.inv(2:end-1)' + O2_ase.ASE(2:end-1) - [O2_ent.ent] - kz_(2:end-1);
NCP_est_kz_DIC = DIC_inv.inv(2:end-1)'  + ADV_DIC(2:end-1)  + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent] -  kz_DIC_(2:end-1);
NCP_est_kz_no_adv_DIC = DIC_inv.inv(2:end-1)' + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent] - kz_DIC_(2:end-1);
NCP_est_no_adv = O2_inv.inv(2:end-1)' + O2_ase.ASE(2:end-1) - [O2_ent.ent];
NCP_est_no_adv_DIC = DIC_inv.inv(2:end-1)' + DIC_ase.FDIC(2:end-1) - [DIC_ent.ent];
% get errors
run 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\NCP_error.m'
plot_days = options.dayrange(2:end-1);

%% run NCP_buoy and create table of vals
% get buoy vals
run NCP_buoy
% get table vales
run NCP_table.m


%% save NCP
save([options.data_dir,'NCP_',num2str(options.h),'m_',num2str(options.window),'_day_window_',num2str(options.window_adv*options.window),'_day_window_ADV.mat']);  
