%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_airseaexchange_DIC.m

% Script to calculate Air-sea gas exchange of O2

% created by MPH in Norwich, 15/11/2017
% modified by MPH in Sydney, 03/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate air sea exchange for DIC

xpCO2.time = [datenum(2016,03,04),datenum(2016,03,11),datenum(2016,03,18),datenum(2016,03,26),datenum(2016,04,1),datenum(2016,04,15)];
xpCO2.val = [408.23,407.24, 408.68, 409.12, 406.31,408.08];
xpCO2.val_int  = interp1(xpCO2.time,xpCO2.val,METEO.date,'Linear'); 
xpCO2.val_int_pCO2  = xpCO2.val_int .*smooth((METEO.sea_level_pressure/100000),5)';
[~,f] = unique(METEO.date);
DIC_ase.pCO2_atm = interp1(METEO.date(f),xpCO2.val_int_pCO2(f),datenum(2016,03,09):1:datenum(2016,04,04),'linear');

[DIC_ase.FDIC, ~,DIC_ase.Sc,DIC_ase.K,DIC_ase.a]=FCO2_updated([means_struct.fCO2_surf], DIC_ase.pCO2_atm,O2_ase.Temp,O2_ase.Salt,O2_ase.wind10,O2_ase.wind10sq); %  squared wind is used for oxygen, need to sqrt

DIC_ase.xpCO2 = xpCO2;

clear xpCO2 f
disp('fCO2 (DIC) air-sea exchange | calculated');
