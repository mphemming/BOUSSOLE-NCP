%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_airsea.m

% Script to calculate Air-sea gas exchange of O2

% created by MPH in Norwich, 15/11/2017
% modified by MPH in Sydney, 05/10/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get time window averaged variables for ASE

O2_ase.Salt = [means_struct.S_surf_profs]; % practical
O2_ase.Temp = [means_struct.T_surf_profs]; %degC
O2_ase.press = [means_struct.sea_level_pressure_atm];
O2_ase.wind10sq = [means_struct.wind_squared]; % mean(wind^2)
O2_ase.wind10 = [means_struct.wind]; % wind normal
% This is already wind speed at 10m, added
% 10% on during data processing - see mistrals.sedoo website

%% oxygen saturation
% Benson and Krause (1984) - equation 24 using tables 5 and 9 and
% atmospheric pressure (atm) from Meteo buoy
% O2 sat taking into account pressure
O2_ase.O2_saturation = o2satSTP(O2_ase.Temp, O2_ase.Salt, O2_ase.press*1013.25);
% convert to mmol m^-3
O2_ase.O2_saturation = ([means_struct.sig0_surf]/1000) .* O2_ase.O2_saturation;

%% Schmidt number
% % The Schmidt number is the kinematic viscosity of water divided by the
% % molecular diffusion coefficient of the gas, should = 568 at 20C in 35
% % salinity water
% % The kinematic viscosity for fresh water and seawater are from Sharqawy et al. (2010). 
% % The diffusion coefficients of gases are from the following: Ar, O2, N2, N2O, and CCl4 fit
% % using Wilke and Chang (1955) as adapted by Hayduk and Laudie (1974).
O2_ase.ScO2 = 1920.4 - (135.6 * O2_ase.Temp) + (5.2122 * O2_ase.Temp.^2) ...
    - (0.10939 * O2_ase.Temp.^3) + (0.0009377 * O2_ase.Temp.^4); 

%% Gas diffusivity
O2_ase.KO2 = 0.251 .* O2_ase.wind10.^2 .* ((O2_ase.ScO2/660).^-0.5); % gas diffusivity

%% Calculate airsea exchange
% gas diffusivity, bubble, and schmidt also calculated in funciton 'ASEflux'
[O2_ase.ASE O2_ase.ASE_uncertainty O2_ase.KO2 O2_ase.Sch O2_ase.bub] = ASEflux(O2_ase.Temp, O2_ase.wind10, O2_ase.wind10sq, ...
    [means_struct.O2_surf_profs],[means_struct.O2_surf_std_profs],O2_ase.O2_saturation,O2_ase.press,1,1);

%% correct ASE for mixing effects

Zlim = options.h;
MLDs = [means_struct.MLD_h_profs];
O2_ase.correction = Zlim ./ MLDs;
O2_ase.correction(O2_ase.correction > 1) = 1;

O2_ase.ASE = O2_ase.ASE.*O2_ase.correction;

disp('O2 air-sea exchange | calculated');


%% calculate air sea exchange for DIC

load([options.directory,'\data\METEO.mat']);  

xpCO2.time = [datenum(2016,03,04),datenum(2016,03,11),datenum(2016,03,18),datenum(2016,03,26),datenum(2016,04,1),datenum(2016,04,15)];
xpCO2.val = [408.23,407.24, 408.68, 409.12, 406.31,408.08];
xpCO2.val_int  = interp1(xpCO2.time,xpCO2.val,METEO.date,'Linear'); 
xpCO2.val_int_pCO2  = xpCO2.val_int .*smooth((METEO.sea_level_pressure/100000),5)';
[~,f] = unique(METEO.date);
DIC_ase.pCO2_atm = interp1(METEO.date(f),xpCO2.val_int_pCO2(f),datenum(2016,03,09):1:datenum(2016,04,04),'linear');

[DIC_ase.FDIC, ~,DIC_ase.Sc,DIC_ase.K,DIC_ase.a]=FCO2_updated([means_struct.fCO2_surf], DIC_ase.pCO2_atm,O2_ase.Temp,O2_ase.Salt,O2_ase.wind10,O2_ase.wind10sq); %  squared wind is used for oxygen, need to sqrt

DIC_ase.xpCO2 = xpCO2;

clear xpCO2 f METEO
disp('fCO2 (DIC) air-sea exchange | calculated');

