%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% NCP_prepare_data.m
%------------
% Script to prepare data for use in NCP calculation
%-----------------------------------------------------------------------------------------------------
% script modified by MPH in Sydney, 05/10/2020
%
%% get variables

vars.dive = prcdata.timeseries.dive_vector;
vars.downup = prcdata.timeseries.downup_vector;
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
vars.oxycline = prcdata.timeseries.oxyclineD; % depth of oxycline
vars.oxyclineO2 = prcdata.timeseries.oxyclineDval; % oxygen concentration at oxycline
vars.O2 =  (prcdata.timeseries.sigma0/1000) .* prcdata.timeseries.O2; % calibrated O2, convert µmol kg -> mmol m-3
% Some further QC for O2
vars.O2(vars.depth < 5 | vars.O2 > 300) = NaN; % remove top 5m where spikes
vars.O2_MLD = vars.O2;
vars.O2_sat = o2satSTP(vars.T, vars.S, 1013); % should use changing pressure!
vars.O2_sat = (prcdata.timeseries.sigma0/1000) .* vars.O2_sat; % convert to mmol
vars.Fl = prcdata.timeseries.Chl_Fl;
vars.Sc700 = prcdata.timeseries.Scatter_700;
% vars.DIC = ([prcdata.CO2SYS.DIC]./vars.S) * 38.3; % normalised DIC for S
vars.DIC = prcdata.timeseries.DIC_recalibrated;
vars.DIC = (prcdata.timeseries.sigma0/1000) .* vars.DIC; % convert µmol kg -> mmol m-3
vars.fCO213 = ([prcdata.CO2SYS.fCO2]).*exp(0.0423*(13-vars.T)); % normalised FCO2 for T
vars.fCO2 = [prcdata.CO2SYS.fCO2];
%% add wind from meteo buoy,
% needs selecting and interpolating onto same time as glider data
wind_selection = METEO.time_wind_speed >= vars.t(1) & ...
    METEO.time_wind_speed <= vars.t(end);
l1 = linspace(1,length(METEO.time_wind_speed(wind_selection)),length(vars.t));
vars.Wind_time = interp1(METEO.time_wind_speed(wind_selection)',l1,'Linear');
vars.Wind = interp1(METEO.wind_speed(wind_selection)',l1,'Linear');
%% atm pressure
press_selection = METEO.date >= vars.t(1) & METEO.date <= vars.t(end);
l1 = linspace(1,length(METEO.date(press_selection)),length(vars.t));
vars.sea_level_press_time = interp1(METEO.date(press_selection)',l1,'Linear');
vars.sea_level_press = interp1(METEO.sea_level_pressure(press_selection)',l1,'Linear');
%% DACs
vars.DACs.DACu = [prcdata.hydrography.DAC_u];
vars.DACs.DACv = [prcdata.hydrography.DAC_v];
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
end
% DAC QC
vars.DACs.DACu(vars.DACs.DACu > 0.2) = NaN;
vars.DACs.DACv(vars.DACs.DACv > 0.2 | vars.DACs.DACv < -0.2) = NaN;
%% satellite measurements

% vars.GVsat_time =  [prcdata.GVsat.day];
% vars.GVsat_U =  [prcdata.GVsat.abs_U_daymean_Glider];
% vars.GVsat_V =  [prcdata.GVsat.abs_V_daymean_Glider];

% AVISO
% obtain area average
for n = 1:numel(prcdata.GVsat)
    AVISO.time(n) = prcdata.GVsat(n).day;
    % get domain average
    U = prcdata.GVsat(n).abs_U';
    V = prcdata.GVsat(n).abs_V';
    f_lon = find(prcdata.GVsat(n).Lon >= 7.4 & prcdata.GVsat(n).Lon <= 8.2);
    f_lat = find(prcdata.GVsat(n).Lat >= 43 & prcdata.GVsat(n).Lat <= 43.6);
    AVISO.absU_ave(n) = nanmean(nanmean(U(f_lat,f_lon)));
    AVISO.absV_ave(n) = nanmean(nanmean(V(f_lat,f_lon)));  
    m_lon = find(prcdata.GVsat(n).Lon == 7.8125);
    m_lat = find(prcdata.GVsat(n).Lat == 43.3125);    
    AVISO.absU_middle(n) = U(m_lat,m_lon);
    AVISO.absV_middle(n) = V(m_lat,m_lon);
end

% Look at JPL Oscar for comparison
%https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplOscar.html      for attributes

load('C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\data\jplOscar.mat.mat');
Oscar.time = (jplOscar.time/60/60/24) + datenum(1970,01,01);
f_lon = find(jplOscar.longitude-360 >= 7.2 & jplOscar.longitude-360 <= 8.4);
f_lat = find(jplOscar.latitude >= 43 & jplOscar.latitude <= 43.7);
m_lon = find(jplOscar.longitude-360 >= 7.4 & jplOscar.longitude-360 <= 8.3);
m_lat = find(jplOscar.latitude >= 43.3 & jplOscar.latitude <= 43.4);
for n = 1:numel(Oscar.time)
    u = squeeze(jplOscar.u(n,1,:,:));
    v = squeeze(jplOscar.v(n,1,:,:));
    % averages
    Oscar.U_ave(n) = nanmean(nanmean(u(f_lat,f_lon)));
    Oscar.V_ave(n) = nanmean(nanmean(v(f_lat,f_lon)));
    % middle
    Oscar.U_middle(n) = nanmean(nanmean(u(m_lat,m_lon)));
    Oscar.V_middle(n) = nanmean(nanmean(v(m_lat,m_lon)));    
end

%% get MLD parameters
MLD = prcdata.timeseries.MLDO2; % MLD from oxygen concentration
load([options.data_dir,'MLD_O2_eye.mat'],'xi','yi')
MLD_O2_eye = yi;
MLD_O2_eye(MLD_O2_eye <= 5) = NaN;
load([options.data_dir,'MLD_T_eye.mat'],'xi','yi')
MLD_T_eye = yi;
MLD_T_eye(MLD_T_eye <= 5) = NaN;

% profile n for use in obtaining profiles
[vars ~] = get_profile_means(vars,options.h,numel(vars.O2));

for n = 1:294
    check = vars.profile_n == n;
    vars.MLD_combined(check) = nanmean([unique(MLD(check)), MLD_T_eye(n),MLD_O2_eye(n)]);
end

for n_profs = 1:294
    check = vars.profile_n == n_profs;
    % remove profs less than 20 m deep
    if nanmax(vars.depth(check)) < 20
        vars.O2(check) = NaN; 
        vars.O2_MLD(check) = NaN;         
    end
    % O2 only in MLD
    MLD = nanmean(vars.MLD_combined(check));
    vars.O2_MLD(check & vars.depth > MLD) = NaN;
end

%% get profile means for layer from data 

[~, vars_profile_means] = get_profile_means(vars,options.h,numel(vars.O2));

%% clear unimportant variables
clearvars -except options prcdata vars vars_profile_means AVISO Oscar