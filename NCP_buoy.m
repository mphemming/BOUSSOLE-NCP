%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_buoy.m

% Script to calculate NCP using buoy DIC and O2

% created by MPH in Norwich, 21/08/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in data and addpaths for functions
% 
% clearvars -except NCP* DIC* O2* errors planes_loop options
% close all
% clc

options.directory = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts'
options.data_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/data/'
options.plot_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/Plots/'
addpath(genpath(options.directory));
load([options.data_dir,'prcdata.mat'],'prcdata');
load([options.data_dir,'METEO.mat']);  
load([options.data_dir,'BOUSSOLE.mat']);  
load([options.data_dir,'closeness.mat']);  
load([options.data_dir,'CTD.mat']);  


%% Get variables for BOUSSOLE using 3 hr bins

BUOY.t = BOUSSOLE.time_CSYS;
BUOY.S = BOUSSOLE.sal_CSYS;
BUOY.ALK = BOUSSOLE.ALK_CSYS;
BUOY.T = BOUSSOLE.temp_CSYS;
BUOY.DIC = BOUSSOLE.DIC_CSYS_mmolm3;
BUOY.fCO2 = BOUSSOLE.fCO2_CSYS;
BUOY.O2 = BOUSSOLE.O2_raw_10m_Liliane;
check = isfinite(BOUSSOLE.time_CSYS) & isfinite(BOUSSOLE.rho_CSYS);
BUOY.dens = interp1(BOUSSOLE.time_CSYS(check),BOUSSOLE.rho_CSYS(check),BOUSSOLE.O2_raw_10m_Liliane_date,'Linear');
BUOY.O2 = (BUOY.dens/1000)' .* BUOY.O2;

BUOY.O2_date = BOUSSOLE.O2_raw_10m_Liliane_date;
BUOY.dv_O2 = datevec(BOUSSOLE.O2_raw_10m_Liliane_date);
BUOY.O2_hr = BUOY.dv_O2(:,4);
BUOY.DICnormS = BOUSSOLE.DIC_CSYS_Snorm_mmolm3;
BUOY.hr = datestr(BUOY.t(1:1080));
BUOY.hr = BUOY.hr(:,13:14);
BUOY.hr = str2num(BUOY.hr)';
BUOY.hr(1081) = NaN;

MET.wind = METEO.wind_speed;
MET.t = METEO.time_wind_speed;
MET.press = METEO.sea_level_pressure;
[C,ia,ic] = unique(METEO.date);
MET.press = interp1(METEO.date(ia),METEO.sea_level_pressure(ia),MET.t);

d = datenum(2016,03,10,00,00,00):datenum(00,00,00,03,00,00):datenum(2016,04,3,00,00,00);

for dl = 1:length(d);
    daynum = d(dl);
    bins(dl).Bt = daynum;
    bins(dl).Bt_str = datestr(daynum);
    bins(dl).BT = nanmedian(BUOY.T(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));
    bins(dl).BS = nanmedian(BUOY.S(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));
    bins(dl).BT_std = nanstd(BUOY.T(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));
    bins(dl).BS_std = nanstd(BUOY.S(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));
    bins(dl).BALK = nanmedian(BUOY.ALK(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));
    bins(dl).BDIC = nanmedian(BUOY.DIC(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));       
    bins(dl).BDICnormS = nanmedian(BUOY.DICnormS(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7));   
    bins(dl).BfCO2 = nanmedian(BUOY.fCO2(BUOY.t >= daynum-((1/48)*3) & BUOY.t <= daynum+((1/48)*3) & BUOY.hr > 4 & BUOY.hr <= 7)); 
     bins(dl).press_std = nanstd(MET.press(MET.t >= daynum-((1/48)*3) & MET.t <= daynum+((1/48)*3)));
    bins(dl).press = nanmedian(MET.press(MET.t >= daynum-((1/48)*3) & MET.t <= daynum+((1/48)*3)));
    bins(dl).wind = nanmedian(MET.wind(MET.t >= daynum-((1/48)*3) & MET.t <= daynum+((1/48)*3)));
    bins(dl).wind10 = bins(dl).wind.*((10/0.7)^0.108);
    bins(dl).O2 = nanmedian(BUOY.O2(BUOY.O2_date >= daynum-((1/48)*3) & BUOY.O2_date <= daynum+((1/48)*3) & BUOY.O2_hr > 4 & BUOY.O2_hr <= 7));
    bins(dl).dens = nanmedian(BUOY.dens(BUOY.O2_date >= daynum-((1/48)*3) & BUOY.O2_date <= daynum+((1/48)*3) & BUOY.O2_hr > 4 & BUOY.O2_hr <= 7));
end


clear d dl daynum BUOY MET 
save([options.directory,'/data/3hrbins'],'bins')

%% get variables ready for NCP calculation

t = [bins.Bt];
T = [bins.BT];
S = [bins.BS];
wind10 = [bins.wind10];
atmpress = [bins.press];
fCO2 = [bins.BfCO2];
DIC = [bins.BDICnormS];

%% DIC NCP
% get atmospheric pCO2 throughout the time period
NCP_DIC.xpCO2.time = [datenum(2016,03,04),datenum(2016,03,11),datenum(2016,03,18),datenum(2016,03,26),datenum(2016,04,1),datenum(2016,04,15)];
NCP_DIC.xpCO2.val = [408.23,407.24, 408.68, 409.12, 406.31,408.08];
NCP_DIC.xpCO2.val_int  = interp1(NCP_DIC.xpCO2.time,NCP_DIC.xpCO2.val,METEO.date,'Linear'); 
NCP_DIC.xpCO2.val_int_pCO2  = NCP_DIC.xpCO2.val_int .*smooth((METEO.sea_level_pressure/100000),5)';
%------------------------------------------------------------------------------------------------------------------------------------
% NCP 20 - 25th March
% get atm pCO2 for time period chosen
check = t > datenum(2016,03,20,00,00,00) & t < datenum(2016,03,25,00,00,00) & ~isnan(t) & ~isnan(DIC);
check_METEO = METEO.date > datenum(2016,03,20,00,00,00) & METEO.date < datenum(2016,03,25,00,00,00);
NCP_DIC.M2025.pCO2atm = interp1(METEO.date(check_METEO),NCP_DIC.xpCO2.val_int_pCO2(check_METEO),t(check),'Linear');
% calculate air-sea exchange
[NCP_DIC.M2025.F, NCP_DIC.M2025.dpCO2]=FCO2_updated(fCO2(check),NCP_DIC.M2025.pCO2atm,T(check),S(check),wind10(check),nanmean(wind10(check).^2));
[NCP_DIC.M2025.fit,NCP_DIC.M2025.gof,NCP_DIC.M2025.out] = fit(t(check)',DIC(check)','poly1');
NCP_DIC.M2025.t = t(check)';
NCP_DIC.M2025.DIC = DIC(check)';
% estimate NCP and advection
NCP_DIC.M2025.NCP = (1.029 * 30 * NCP_DIC.M2025.fit.p1 + nanmean(NCP_DIC.M2025.F)) ; % mmol m^-3 d-1
NCP_DIC.M2025.ADV = [DIC_adv.adv] .* ([planes_loop.date_num] >= datenum(2016,03,20) & [planes_loop.date_num] <= datenum(2016,03,25));
NCP_DIC.M2025.NCP_ADV = (NCP_DIC.M2025.NCP + nanmean(NCP_DIC.M2025.ADV(NCP_DIC.M2025.ADV ~= 0))); % mmol m^-3 d-1
NCP_DIC.M2025.check = check;
%------------------------------------------------------------------------------------------------------------------------------------
% NCP 29/03 - 01/04
% get atm pCO2 for time period chosen
check = t > datenum(2016,03,29,00,00,00) & t < datenum(2016,04,02,00,00,00) & ~isnan(t) & ~isnan(DIC);
check_METEO = METEO.date > datenum(2016,03,29,00,00,00) & METEO.date < datenum(2016,04,02,00,00,00);
tm = METEO.date(check_METEO);
valm = NCP_DIC.xpCO2.val_int_pCO2(check_METEO);
[f f1 f2] = unique(tm);
NCP_DIC.M2901.pCO2atm = interp1(tm(f1),valm(f1),t(check),'Linear');
% calculate air-sea exchange
[NCP_DIC.M2901.F, NCP_DIC.M2901.dpCO2]=FCO2_updated(fCO2(check),NCP_DIC.M2901.pCO2atm,T(check),S(check),wind10(check),nanmean(wind10(check).^2));
[NCP_DIC.M2901.fit,NCP_DIC.M2901.gof,NCP_DIC.M2901.out] = fit(t(check)',DIC(check)','poly1');
NCP_DIC.M2901.t = t(check)';
NCP_DIC.M2901.DIC =DIC(check)';
% estimate NCP and advection
NCP_DIC.M2901.NCP = (1.029 * 22 * NCP_DIC.M2901.fit.p1 + nanmean(NCP_DIC.M2901.F)) ; % mmol m^-3 d-1
NCP_DIC.M2901.ADV = [DIC_adv.adv] .* ([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,01));
NCP_DIC.M2901.NCP_ADV = NCP_DIC.M2901.NCP + nanmean(NCP_DIC.M2901.ADV(NCP_DIC.M2901.ADV ~= 0)); % mmol m^-3 d-1
NCP_DIC.M2901.check = check;
    
%% O2 NCP
t = [bins.Bt];
O2 = [bins.O2];
dens =[bins.dens];
T = [bins.BT];

% calculate O2 ASE

% oxygen saturation
% Benson and Krause (1984) - equation 24 using tables 5 and 9 and
% atmospheric pressure (atm) from Meteo buoy
% O2 sat taking into account pressure
O2_saturation = o2satSTP([bins.BT], [bins.BS], atmpress/100);
% convert to mmol m^-3
O2_saturation = (dens/1000) .* O2_saturation;
% schmidt
ScO2 = 1920.4 - (135.6 * T) + (5.2122 * T.^2) ...
    - (0.10939 * T.^3) + (0.0009377 * T.^4); 
% Gas diffusivity
KO2 = 0.251 .* wind10.^2 .* ((ScO2/660).^-0.5); % gas diffusivity
% Calculate airsea exchange
% gas diffusivity, bubble, and schmidt also calculated in funciton 'ASEflux'
[ASE ASE_uncertainty] = ASEflux(T, wind10, wind10.^2, ...
    O2,ones(size(O2)),O2_saturation,atmpress,1,1);

%------------------------------------------------------------------------------------------------------------------------------------
% NCP 20 - 25th March
check = t > datenum(2016,03,20,00,00,00) & t < datenum(2016,03,25,00,00,00) & ~isnan(t) & ~isnan(O2);
[NCP_O2.M2025.fit,NCP_O2.M2025.gof,NCP_O2.M2025.out] = fit(t(check)',O2(check)','poly1');
NCP_O2.M2025.t = t(check)';
NCP_O2.M2025.O2 =O2(check)';
% estimate NCP and advection
NCP_O2.M2025.ASE = nanmean(ASE(check));
% previously 46 m, can use option.h, or now trying mean MLD (options.h m)
NCP_O2.M2025.NCP = 1.029 * 30 * NCP_O2.M2025.fit.p1 + NCP_O2.M2025.ASE; % mmol m^-3 d-1
NCP_O2.M2025.ADV = [O2_adv.adv] .* ([planes_loop.date_num] > datenum(2016,03,20) & [planes_loop.date_num] < datenum(2016,03,25));
NCP_O2.M2025.NCP_ADV = NCP_O2.M2025.NCP + nanmean(NCP_O2.M2025.ADV(NCP_O2.M2025.ADV ~= 0)); % mmol m^-3 d-1
NCP_O2.M2025.check = check;    
%------------------------------------------------------------------------------------------------------------------------------------
% NCP 29/03 - 01/04
check = t > datenum(2016,03,29,00,00,00) & t < datenum(2016,04,02,00,00,00) & ~isnan(t) & ~isnan(O2);
[NCP_O2.M2901.fit,NCP_O2.M2901.gof,NCP_O2.M2901.out] = fit(t(check)',O2(check)','poly1');
NCP_O2.M2901.t = t(check)';
NCP_O2.M2901.O2 =O2(check)';
% estimate NCP and advection
NCP_O2.M2901.ASE = nanmean(ASE(check));
NCP_O2.M2901.NCP = 1.029 * 22 * NCP_O2.M2901.fit.p1 + NCP_O2.M2901.ASE; % mmol m^-3 d-1
NCP_O2.M2901.ADV = [O2_adv.adv] .* ([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,02));
NCP_O2.M2901.NCP_ADV = NCP_O2.M2901.NCP + nanmean(NCP_O2.M2901.ADV(NCP_O2.M2901.ADV ~= 0)); % mmol m^-3 d-1
NCP_O2.M2901.check = check;     
    
%% errors for buoy NCP

% get standard error for fits
% DIC

%------------------------------------------------------------------------------------------------------------------------------------
% buoy_error_DIC.SE_line_20_25.m = NCP_DIC.M2025.fit.p1;
% buoy_error_DIC.SE_line_20_25.t = NCP_DIC.M2025.t;
% buoy_error_DIC.SE_line_20_25.DIC = NCP_DIC.M2025.DIC;
% buoy_error_DIC.SE_line_20_25.DIC_pred = NCP_DIC.M2025.fit(buoy_error_DIC.SE_line_20_25.t);
% buoy_error_DIC.SE_line_20_25.residuals = buoy_error_DIC.SE_line_20_25.DIC_pred- buoy_error_DIC.SE_line_20_25.DIC;
% buoy_error_DIC.SE_line_20_25.t2 = buoy_error_DIC.SE_line_20_25.t.^2;
% buoy_error_DIC.SE_line_20_25.residuals2 = buoy_error_DIC.SE_line_20_25.residuals.^2;
% buoy_error_DIC.SE_line_20_25.sum_residuals2 = nansum(buoy_error_DIC.SE_line_20_25.residuals2);
% buoy_error_DIC.SE_line_20_25.sum_t2 = nansum(buoy_error_DIC.SE_line_20_25.t2);
% buoy_error_DIC.SE_line_20_25.sum_t2_n = buoy_error_DIC.SE_line_20_25.sum_t2+10;
% buoy_error_DIC.SE_line_20_25.sum_t = nansum(buoy_error_DIC.SE_line_20_25.t);
% buoy_error_DIC.SE_line_20_25.sum_t_2 = nansum(buoy_error_DIC.SE_line_20_25.t).^2;
% buoy_error_DIC.SE_line_20_25.sigma_m_2 = buoy_error_DIC.SE_line_20_25.sum_residuals2 / ...
%     (buoy_error_DIC.SE_line_20_25.sum_t2_n -buoy_error_DIC.SE_line_20_25.sum_t_2);
% buoy_error_DIC.SE_line_20_25.sigma = real(sqrt(buoy_error_DIC.SE_line_20_25.sigma_m_2));
%------------------------------------------------------------------------------------------------------------------------------------
% buoy_error_DIC.SE_line_29_01.m = NCP_DIC.M2901.fit.p1;
% buoy_error_DIC.SE_line_29_01.t = NCP_DIC.M2901.t;
% buoy_error_DIC.SE_line_29_01.DIC = NCP_DIC.M2901.DIC;
% buoy_error_DIC.SE_line_29_01.DIC_pred = NCP_DIC.M2901.fit(buoy_error_DIC.SE_line_29_01.t);
% buoy_error_DIC.SE_line_29_01.residuals = buoy_error_DIC.SE_line_29_01.DIC_pred- buoy_error_DIC.SE_line_29_01.DIC;
% buoy_error_DIC.SE_line_29_01.t2 = buoy_error_DIC.SE_line_29_01.t.^2;
% buoy_error_DIC.SE_line_29_01.residuals2 = buoy_error_DIC.SE_line_29_01.residuals.^2;
% buoy_error_DIC.SE_line_29_01.sum_residuals2 = nansum(buoy_error_DIC.SE_line_29_01.residuals2);
% buoy_error_DIC.SE_line_29_01.sum_t2 = nansum(buoy_error_DIC.SE_line_29_01.t2);
% buoy_error_DIC.SE_line_29_01.sum_t2_n = buoy_error_DIC.SE_line_29_01.sum_t2+10;
% buoy_error_DIC.SE_line_29_01.sum_t = nansum(buoy_error_DIC.SE_line_29_01.t);
% buoy_error_DIC.SE_line_29_01.sum_t_2 = nansum(buoy_error_DIC.SE_line_29_01.t).^2;
% buoy_error_DIC.SE_line_29_01.sigma_m_2 = buoy_error_DIC.SE_line_29_01.sum_residuals2 / ...
%     (buoy_error_DIC.SE_line_29_01.sum_t2_n -buoy_error_DIC.SE_line_29_01.sum_t_2);
% buoy_error_DIC.SE_line_29_01.sigma = real(sqrt(buoy_error_DIC.SE_line_29_01.sigma_m_2));

% buoy_error_DIC.SE_line_20_25 = NCP_DIC.M2025.gof.rmse;
% buoy_error_DIC.SE_line_29_01 = NCP_DIC.M2901.gof.rmse;
% buoy_error_DIC.ase_20_25 = nanmean(...
%     errors.ASE.errors_ASE_DIC([planes_loop.date_num] > datenum(2016,03,20) & [planes_loop.date_num] < datenum(2016,03,25)))
% buoy_error_DIC.ase_29_01 = nanmean(...
%     errors.ASE.errors_ASE_DIC([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,02)))
% EADV = [errors.ADV_DIC.errors_adv].*([planes_loop.date_num] > datenum(2016,03,20) & [planes_loop.date_num] < datenum(2016,03,25));
% EADV(EADV == 0) = [];
% buoy_error_DIC.adv_20_25 = nanmean(EADV);
% EADV = [errors.ADV_DIC.errors_adv].*([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,01));
% EADV(EADV == 0) = [];
% buoy_error_DIC.adv_29_01 = nanmean(EADV);
% 
% buoy_error_DIC.error_20_25_adv = sqrt ( ((buoy_error_DIC.SE_line_20_25).^2) + ((buoy_error_DIC.ase_20_25).^2) + ((buoy_error_DIC.adv_20_25).^2) );
% buoy_error_DIC.error_29_01_adv = sqrt ( ((buoy_error_DIC.SE_line_29_01).^2) + ((buoy_error_DIC.ase_29_01).^2) + ((buoy_error_DIC.adv_29_01).^2) );
% buoy_error_DIC.error_20_25 = sqrt ( ((buoy_error_DIC.SE_line_20_25).^2) + ((buoy_error_DIC.ase_20_25).^2) );
% buoy_error_DIC.error_29_01 = sqrt (  ((buoy_error_DIC.SE_line_29_01).^2) + ((buoy_error_DIC.ase_29_01).^2) );

% O2
% 
% buoy_error_O2.SE_line_20_25 = NCP_O2.M2025.gof.rmse;
% buoy_error_O2.SE_line_29_01 = NCP_O2.M2901.gof.rmse;
% % ASE errors
% buoy_error_O2.o2sat = o2satSTP([bins.BT], [bins.BS], atmpress/100) - ...
%     o2satSTP([bins.BT]+[bins.BT_std], [bins.BS]+[bins.BS_std], (atmpress+[bins.press_std])/100);
% buoy_error_O2.o2sat = (dens/1000) .* buoy_error_O2.o2sat; 
% buoy_error_O2.KO2 = ASE_uncertainty.ASEAlkireKO2val/100*24;
% buoy_error_O2.ase_20_25= sqrt((buoy_error_O2.KO2 .* buoy_error_O2.o2sat).^2);  
% buoy_error_O2.ase_20_25 = nanmean(buoy_error_O2.ase_20_25(t > datenum(2016,03,20,00,00,00) & t < datenum(2016,03,25,00,00,00)));
% buoy_error_O2.ase_29_01= sqrt((buoy_error_O2.KO2 .* buoy_error_O2.o2sat).^2);  
% buoy_error_O2.ase_29_01 = nanmean(buoy_error_O2.ase_29_01(t > datenum(2016,03,29,00,00,00) & t < datenum(2016,04,02,00,00,00)));
% 
% % buoy_error_O2.ase_20_25 = nanmean(...
% %     errors.ASE.errors_ASE([planes_loop.date_num] > datenum(2016,03,20) & [planes_loop.date_num] < datenum(2016,03,25)));
% % buoy_error_O2.ase_29_01 = nanmean(...
% %     errors.ASE.errors_ASE([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,02)));
% EADV = [errors.ADV.errors_adv].*([planes_loop.date_num] > datenum(2016,03,20) & [planes_loop.date_num] < datenum(2016,03,25));
% EADV(EADV == 0) = [];
% buoy_error_O2.adv_20_25 = nanmean(EADV);
% EADV = [errors.ADV.errors_adv].*([planes_loop.date_num] > datenum(2016,03,29) & [planes_loop.date_num] < datenum(2016,04,01));
% EADV(EADV == 0) = [];
% buoy_error_O2.adv_29_01 = nanmean(EADV);
%     
% buoy_error_O2.error_20_25_adv = sqrt ( ((buoy_error_O2.SE_line_20_25).^2) + ((buoy_error_O2.ase_20_25).^2) + ((buoy_error_O2.adv_20_25).^2) );
% buoy_error_O2.error_29_01_adv = sqrt ( ((buoy_error_O2.SE_line_29_01).^2) + ((buoy_error_O2.ase_29_01).^2) + ((buoy_error_O2.adv_29_01).^2) );
% buoy_error_O2.error_20_25 = sqrt ( ((buoy_error_O2.SE_line_20_25).^2) + ((buoy_error_O2.ase_20_25).^2) );
% buoy_error_O2.error_29_01 = sqrt (  ((buoy_error_O2.SE_line_29_01).^2) + ((buoy_error_O2.ase_29_01).^2) );
    
%% clear uneccessary variables

clearvars -except NCP* DIC* O2* errors planes_loop options buoy_error* NCP_O2 NCP_DIC vars means_struct

    
    
    



