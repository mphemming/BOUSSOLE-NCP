%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% pHcalibration.m
%------------
% calibrate pH using offset correction using CTD profiles
%-----------------------------------------------------------------------------------------------------
% script by MPH in Sydney, 01\11\2020
% script modified by MPH in Sydney, 01\11\2020
%

clear all
close all
clc

%% load data

options.directory = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts'
options.data_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\data\'
options.plot_dir = 'C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_scripts\Plots\'
load C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\closeness.mat;
load C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\BOUSSOLE_data\BOUSSOLE.mat;
load C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\Meteo_Buoy_data\METEO.mat;  
load('C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\tseries');
load('C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\distBOUSSOLE');
load C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\BOUSSOLE_data\CTD.mat;
load C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\DYFAMED_data\DYFAMED.mat;
load([options.directory,'\data\prcdata.mat'],'prcdata');
load('C:\Users\mphem\Documents\Work\UEA\UEA_work\BOUSSOLE_data\BOUSSOLE_data\BOUSSOLE.mat','BOUSSOLE');

%% pH cal notes

% pH cal steps before this point:
% o   TRIS and AMP correction for temperature
% o   Calibrate pH using TRIS and AMP pH




%% depth offset correction using ship CTD profiles

% interpolate to get pH at 1000m before and after deployment

March_pH = interp1(CTD.pH_pCO2_depth,CTD.pH,950,'Linear');
April_pH = interp1(CTD.pH_pCO2_depth_A,CTD.pH_A,950,'Linear');

March_pH_int_prof = interp1(CTD.pH_pCO2_depth,CTD.pH,1:1:1000,'Linear');
April_pH_int_prof = interp1(CTD.pH_pCO2_depth_A,CTD.pH_A,1:1:1000,'Linear');

% use CTD profiles either side of deployment to determine offset

pH_vals = [March_pH, April_pH];
dates = [CTD.time-0.0833, datenum(2016,04,16,12,11,00)]; %changed March time to UTC again, not sure if glider data UTC or local

for dive = 1:147
    d(dive) = nanmedian(prcdata.hydrographyCl(dive).time);
    offset_pH_950m(dive) = interp1(dates,pH_vals,d(dive),'Linear');
    
    % get unique pH vals
    D = prcdata.hydrographyCl(dive).depth;
    pH = prcdata.pH_info_cal(dive).pH_driftcorr;
    [~, un_pH] = unique(pH);  
    pH = pH(un_pH);
    D = D(un_pH);
    % get unique depths
    [~, un_D] = unique(D);  
    pH = pH(un_D);
    D = D(un_D);   
    % ascending order
    [~, index] = sort(D);
    D = D(index);
    pH = pH(index);
    
    c = isfinite(pH);
    
    if max(D(c) >= 900)
        pH_glider_950m(dive) = interp1(D(c),pH(c),950,'Linear');
        prcdata.pH_info_cal(dive).pHdrift_offset_CTDs = ones(size(prcdata.hydrographyCl(dive).depth)).*(offset_pH_950m(dive)-pH_glider_950m(dive));
        prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs = ...
        prcdata.pH_info_cal(dive).pH_driftcorr + prcdata.pH_info_cal(dive).pHdrift_offset_CTDs;
    else
        prcdata.pH_info_cal(dive).pHdrift_offset_CTDs = 0;
        prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs = ...
        prcdata.pH_info_cal(dive).pH_driftcorr + prcdata.pH_info_cal(dive).pHdrift_offset_CTDs;        
    end
    
end

prcdata.pH_info_cal(147).pH_driftcorrFull_CTDs = NaN(size(prcdata.pH_info_cal(147).pHdrift_offset_CTDs));

%% get interpolated ship profiles 

time_ship = [ones(size(March_pH_int_prof))*dates(1), ones(size(April_pH_int_prof))*dates(2)];
pH_ship = [March_pH_int_prof, April_pH_int_prof];
D_ship = [1:1000,1:1000];

[X, Y] = meshgrid(d,1:1:1000);
ship_profs = griddata(time_ship',D_ship',pH_ship',X,Y);

%% pH QC
depth_grid = 10:10:1000;

for dive = 2:146
    % to calculate differences between last up and next down, and at the apogee
        % top diff
        top_t(dive) = nanmedian(prcdata.hydrography(dive).time(prcdata.hydrography(dive).downup_vector == 2 & prcdata.hydrography(dive).depth <=5));
        top_diff(dive) = nanmedian(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(...
            prcdata.hydrography(dive).downup_vector == 2 & prcdata.hydrography(dive).depth <=5)) ...
            - nanmedian(prcdata.pH_info_cal(dive+1).pH_driftcorrFull_CTDs(...
            prcdata.hydrography(dive+1).downup_vector == 1 & prcdata.hydrography(dive+1).depth <=5));
        % bottom diff
        maxD(dive) = nanmax(prcdata.hydrography(dive).depth);
        bott_diff(dive) = nanmedian(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(...
            prcdata.hydrography(dive).downup_vector == 1 & prcdata.hydrography(dive).depth > maxD(dive)-5)) ...
            - nanmedian(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(...
            prcdata.hydrography(dive).downup_vector == 2 & prcdata.hydrography(dive).depth > maxD(dive)-5));   
        % calculate difference profiles
        for n = 1:numel(depth_grid)
            c_down = prcdata.hydrography(dive).downup_vector == 1 & ...
                prcdata.hydrography(dive).depth > depth_grid(n)-5 & ...
                prcdata.hydrography(dive).depth <= depth_grid(n)+5;
            c_up = prcdata.hydrography(dive).downup_vector == 2 & ...
                prcdata.hydrography(dive).depth > depth_grid(n)-5 & ...
                prcdata.hydrography(dive).depth <= depth_grid(n)+5;                
            diff_profiles(dive).depth(n) = depth_grid(n);
            diff_profiles(dive).diff(n) = nanmedian(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(c_down)) - ...
                nanmedian(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(c_up));
        end
end
% mean difference
diff_arr = [diff_profiles.diff];
depth_arr = [diff_profiles.depth];
for n = 1:numel(depth_grid)
    c = depth_arr == depth_grid(n);
    mean_diff(n) = nanmean(diff_arr(c));
    median_diff(n) = nanmedian(diff_arr(c));
    std_diff(n) = nanstd(diff_arr(c));
end


% for dive = 1:147
%     check = prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs < 7.95 | ...
%         prcdata.hydrographyCl(dive).depth > 980 | prcdata.hydrographyCl(dive).depth < 4;
%     prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(check) = NaN;
%     % use downcast only
%     prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(prcdata.hydrography(dive).downup_vector == 2) == NaN;    
% end
depth_grid = 0:20:1000;
% use median profiles instead
for dive = 1:147
    check = prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs < 7.95 | ...
        prcdata.hydrographyCl(dive).depth > 980 | prcdata.hydrographyCl(dive).depth < 4;
    prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(check) = NaN;
    % get median profile
    for n = 1:numel(depth_grid)
        c = prcdata.hydrography(dive).depth > depth_grid(n)-10 & ...
            prcdata.hydrography(dive).depth <= depth_grid(n)+10;    
        pH_bin(n) = nanmean(prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(c));
    end
    c = prcdata.hydrography(dive).downup_vector == 1;
    prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(c) = interp1(depth_grid,pH_bin,prcdata.hydrographyCl(dive).depth(c),'Linear');
    % flag upcast
    c = prcdata.hydrography(dive).downup_vector == 2;
    prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs(c) = NaN;    
end



%% Multiple linear regression for temperature, salinity, and pressure

% get delta pH
for dive = 1:147
        d(dive) = nanmedian(prcdata.hydrographyCl(dive).time); 
        % depth-binning
        T = prcdata.hydrographyCl(dive).temp;
        S = prcdata.hydrographyCl(dive).salinity;
        D = prcdata.hydrographyCl(dive).depth;
        P = prcdata.hydrographyCl(dive).pressure;
        pH = prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs;
        for n = 1:1000
            check = D > n-0.5 & D <= n+0.5;
            bin(dive).T(n) = nanmedian(T(check));
            bin(dive).S(n) = nanmedian(S(check));
            bin(dive).D(n) = n;
            bin(dive).P(n) = nanmedian(P(check));   
            bin(dive).pH(n) = nanmedian(pH(check));   
            bin(dive).t(n) = d(dive);
        end
        delta_pH(dive).vals = bin(dive).pH' - ship_profs(:,dive);
end

t_b = [bin.t]; t_b = t_b(:);
D_b = [bin.D]; D_b = D_b(:);
T_b = [bin.T]; T_b = T_b(:);
S_b = [bin.S]; S_b = S_b(:);
P_b = [bin.P]; P_b = P_b(:);
delta_pH = [delta_pH.vals]; delta_pH = delta_pH(:);
delta_pH(delta_pH < -0.1 | delta_pH > 0.19) = NaN;

% % regression of temperature using data > 350m only
% 
% c = isfinite(delta_pH) & isfinite(T_b) & isfinite(P_b) & isfinite(S_b) & D_b > 350;
% % [fit_T,gof_T,out_T] = fit(T_b(c),delta_pH(c),'poly1');
% 
% [b,bint,r,rint,stats] = regress(delta_pH(c),[ones(size(T_b(c))),T_b(c),P_b(c), S_b(c)]);
% 
% % correct pH for temp, press, sal
% for dive = 1:147
%     prcdata.pH_info_cal(dive).pH_driftcorrFull_TP_CTDs = ...
%         prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs - (b(1) + b(2)*prcdata.hydrographyCl(dive).temp +...
%         b(3)*prcdata.hydrographyCl(dive).pressure + b(4)*prcdata.hydrographyCl(dive).salinity);
% end


%% try using delta pH first week vs CTD and and last week vs CTD

% first week

c = isfinite(delta_pH) & isfinite(T_b) & isfinite(P_b) & isfinite(S_b) & t_b < datenum(2016,03,14);
[first.b,first.bint,first.r,first.rint,first.stats] = regress(delta_pH(c),[ones(size(T_b(c))),T_b(c),P_b(c), S_b(c)]);

% last week
c = isfinite(delta_pH) & isfinite(T_b) & isfinite(P_b) & isfinite(S_b) & t_b > datenum(2016,03,29);
[last.b,last.bint,last.r,last.rint,last.stats] = regress(delta_pH(c),[ones(size(T_b(c))),T_b(c),P_b(c), S_b(c)]);

% calculate delta pH for each dive using coefficients at start and end of
% deployment

% grid coefficients

b_1 = interp1([datenum(2016,03,07) datenum(2016,04,16)],[first.b(1) last.b(1)],d);
b_2 = interp1([datenum(2016,03,07) datenum(2016,04,16)],[first.b(2) last.b(2)],d);
b_3 = interp1([datenum(2016,03,07) datenum(2016,04,16)],[first.b(3) last.b(3)],d);
b_4 = interp1([datenum(2016,03,07) datenum(2016,04,16)],[first.b(4) last.b(4)],d);

for dive = 1:147
    prcdata.pH_info_cal(dive).pH_driftcorrFull_TP_CTDs = ...
        prcdata.pH_info_cal(dive).pH_driftcorrFull_CTDs - (b_1(dive) + b_2(dive)*prcdata.hydrographyCl(dive).temp +...
        b_3(dive)*prcdata.hydrographyCl(dive).pressure + b_4(dive)*prcdata.hydrographyCl(dive).salinity);
end

% remove that one errorneous profile
prcdata.timeseries.pH_recalibrated = [prcdata.pH_info_cal.pH_driftcorrFull_TP_CTDs];
prcdata.timeseries.pH_recalibrated(prcdata.timeseries.dive_vector == 147) = NaN;
prcdata.timeseries.pH_recalibrated(prcdata.timeseries.dive_vector == 146) = NaN;

%% NEED TO DO FURTHER QC

% DIFFERENCE BETWEEN UP/DOWN CASTS
% OUTLIERS AT VARIOUS DEPTHS

% Down seems more reliable than up in the top 45 m - plot time vs pH for up
% and down for 10-20m
% only use down casts for everything?
% upcasts also have drift issue as function of pressure towards bottom

% or use median for each dive to get profiles
% try below:

% for dive = 1:147
%     for n = 10:10:1000
%         check = prcdata.hydrographyCl(dive).depth > n-5 &  prcdata.hydrographyCl(dive).depth <= n+5;
%         medvals(dive).pH(n) = nanmean(prcdata.pH_info_cal(dive).pH_driftcorrFull_TP_CTDs(check));
%         medvals(dive).D(n) = n;
%         medvals(dive).t(n) = d(dive);
%     end  
%     medvals(dive).pH(medvals(dive).pH == 0) = NaN;
%     medvals(dive).D(medvals(dive).D == 0) = NaN;
% end
% 
% mvpH = [medvals.pH]; mvpH = mvpH(:);
% mvD = [medvals.D]; mvD = mvD(:);
% mvt = [medvals.t]; mvt = mvt(:);


%% Derive Alkalinity and DIC using CO2SYS

addpath('C:\Users\mphem\Documents\Work\UEA\UEA_work\NCP_Scripts\Utilities')
[RESULT,~,~]=CO2SYS(prcdata.timeseries.Alk,prcdata.timeseries.pH_recalibrated,1,3,...
    prcdata.timeseries.S,prcdata.timeseries.T,0,prcdata.timeseries.P,0,0,0,1,4,1);

prcdata.timeseries.DIC_recalibrated = RESULT(:,2)';
prcdata.timeseries.pCO2_recalibrated = RESULT(:,4)';
prcdata.timeseries.fCO2_recalibrated = RESULT(:,5)';

save([options.directory,'\data\prcdata.mat'],'prcdata');

%% create comparison between BOUSSOLE buoy

t_grid = datenum(2016,03,07):1:datenum(2016,04,05);

for n = 1:numel(t_grid)
    c = prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5 & prcdata.timeseries.depth > 5 & ...
        prcdata.timeseries.depth < 15 & prcdata.timeseries.close_vector <= 5;
    binned.DIC(n) = nanmedian(prcdata.timeseries.DIC(c));
    binned.DIC_recalibrated(n) = nanmedian(prcdata.timeseries.DIC_recalibrated(c));
end

%% pH calibration plot

figure('units','normalized','position',[0 0 .8 .9]);

% Drift 
drifted_pH = [prcdata.pH_info_cal.pH_driftcorr];
c10 = prcdata.timeseries.depth > 5 & prcdata.timeseries.depth < 15;
c100 = prcdata.timeseries.depth > 95 & prcdata.timeseries.depth < 105;
c250 = prcdata.timeseries.depth > 245 & prcdata.timeseries.depth < 255;
c500 = prcdata.timeseries.depth > 495 & prcdata.timeseries.depth < 505;
c750 = prcdata.timeseries.depth > 745 & prcdata.timeseries.depth < 755;
c1000 = prcdata.timeseries.depth > 995 & prcdata.timeseries.depth < 1005;
% daily resolution
for n = 1:numel(t_grid)
    check = c10 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH10(n) = nanmedian(drifted_pH(check));
    check = c100 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH100(n) = nanmedian(drifted_pH(check));    
    check = c250 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH250(n) = nanmedian(drifted_pH(check));    
    check = c500 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH500(n) = nanmedian(drifted_pH(check)); 
    check = c750 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH750(n) = nanmedian(drifted_pH(check)); 
    check = c1000 & prcdata.timeseries.t > t_grid(n) - 0.5 & prcdata.timeseries.t <= t_grid(n)+0.5;
    binned.pH1000(n) = nanmedian(drifted_pH(check));    
end

% get historical data

fields = fieldnames(DYFAMED.historical);
np = 0;
for n = 1:numel(fields)
    fd = fields(n);
    if ~isempty(strfind(cell2mat(fd),'pH')) & isempty(strfind(cell2mat(fd),'25'))
        insert = cell2mat(fd);
        np = np+1;
        historical(np).pH = eval(['DYFAMED.historical.',insert,';']); 
        historical(np).P = eval(['DYFAMED.historical.',insert(1:6),'P;']);   
        historical(np).fields = insert(1:6);
    end
end

% Create axes
axes('Parent',gcf,...
    'Position',[0.08056640625 0.589248971193416 0.383463541666667 0.335751028806585]);

p1 = plot(t_grid,binned.pH10,'LineWidth',2)
hold on
p2 = plot(t_grid,binned.pH100,'LineWidth',2)
p3 = plot(t_grid,binned.pH250,'LineWidth',2)
p4 = plot(t_grid,binned.pH500,'LineWidth',2)
p5 = plot(t_grid,binned.pH750,'LineWidth',2)
p6 = plot(t_grid,binned.pH1000,'LineWidth',2)

set(gca,'FontSize',16,'LineWidth',2,'Box','On');
ylim([7.5 7.95]); xlim([datenum(2016,03,06) datenum(2016,04,06)]);
datetick('x','dd/mm','KeepLimits'); xlabel('Date'); ylabel('pH_{T,drift}');
leg = legend([p1 p2 p3 p4 p5 p6], '10','100','250','500','750','1000');
set(leg,'Location','SouthWest','Orientation','Horizontal','Box','Off','FontSize',11);

annotation(gcf,'textbox',...
    [0.09084375 0.878600823045268 0.04978125 0.0411522633744856],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off');

% Create axes
axes('Parent',gcf,...
    'Position',[0.08056640625 0.152006172839506 0.384114583333333 0.351851851851852]);

hold on
% raw pH
raw_pH = vertcat(prcdata.pH_info_cal.pH_before);
[N,x] = hist(raw_pH-nanmean(raw_pH),100);
N = (N-nanmin(N)) / (nanmax(N) - nanmin(N));
a = area(x,N);
set(a,'FaceAlpha',0.8);
% pH drift
pH_drift = [prcdata.pH_info_cal.pH_driftcorr];
[N,x] = hist(pH_drift-nanmean(pH_drift),100);
N = (N-nanmin(N)) / (nanmax(N) - nanmin(N));
a = area(x,N);
set(a,'FaceAlpha',0.7);

% pH drift offset
pH_drift_offset = [prcdata.pH_info_cal.pH_driftcorrFull_CTDs];
[N,x] = hist(pH_drift_offset-nanmean(pH_drift_offset),100);
N = (N-nanmin(N)) / (nanmax(N) - nanmin(N));
a = area(x,N);
set(a,'FaceAlpha',0.7);

% pH correct
pH_corr = prcdata.timeseries.pH_recalibrated;
[N,x] = hist(pH_corr-nanmean(pH_corr),100);
N = (N-nanmin(N)) / (nanmax(N) - nanmin(N));
a = area(x,N);
set(a,'FaceAlpha',0.7);

% pH historical
pH_hist = vertcat(historical.pH);
[N,x] = hist(pH_hist-nanmean(pH_hist),100);
N = (N-nanmin(N)) / (nanmax(N) - nanmin(N));
a = area(x,N,'FaceColor','k');
set(a,'FaceAlpha',0.3);

set(gca,'FontSize',16,'LineWidth',2,'Box','On');
xlim([-0.5 0.7]);
leg = legend('pH_T','pH_{T,drift}','pH_{T,drift,offset}','pH_{T,drift,offset,T,P,S}','pH_{T,historical}');
set(leg,'FontSize',10,'Box','Off');
xlabel('pH_T anomaly');
ylabel('Normalised Count');



% Create axes
axes('Parent',gcf,...
    'Position',[0.54931640625 0.152006172839506 0.414713541666667 0.815]);

hold on
s0 = scatter(raw_pH,prcdata.timeseries.depth,3,'filled');
s1 = scatter(pH_drift,prcdata.timeseries.depth,3,'filled');
s2 = scatter(pH_drift_offset,prcdata.timeseries.depth,3,'filled');
s3 = scatter(pH_corr,prcdata.timeseries.depth,3,'filled');
s4 = scatter(pH_hist,vertcat(historical.P),10,'k','filled');
s5 = scatter(vertcat(CTD.pH,CTD.pH_A),vertcat(CTD.pH_pCO2_depth,CTD.pH_pCO2_depth_A),20,'MarkerFaceColor','W','MarkerEdgeColor','k');

set(gca,'YDir','Reverse','LineWidth',2,'Box','On','FontSize',16);
ylim([0 980]); xlim([7.55 8.8]);
xlabel('pH_T'); ylabel('Depth [m]');
leg = legend('pH_T','pH_{T,drift}','pH_{T,drift,offset}','pH_{T,drift,offset,T,P,S}','pH_{T,historical}','pH_{T,ship}');
set(leg,'FontSize',10,'Location','SouthEast');

% Create textbox
annotation(gcf,'textbox',...
    [0.0787994791666667 0.459362139917696 0.0497812500000001 0.0411522633744856],...
    'String','(b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.550153645833333 0.911008230452675 0.0497812500000001 0.0411522633744856],...
    'String','(c)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off')

print(gcf, '-dpng','-r400', [options.plot_dir,'Fig_pHcal.png'])

%% Figure comparing glider O2 and pH with ship profiles

figure('units','normalized','position',[0 0 .8 .9]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.16 0.414108072916667 0.789454732510288]);

vars.O2 =  (prcdata.timeseries.sigma0/1000) .* prcdata.timeseries.O2; % calibrated O2, convert µmol kg -> mmol m-3
% Some further QC for O2
vars.O2(prcdata.timeseries.depth < 5 | vars.O2 > 300) = NaN; % remove top 5m where spikes
O2_down = ((CTD.sigmatheta_do+1000)/1000).*CTD.O2_do;
sig_A = interp1(CTD.P_April_bottle, CTD.sig0_April_bottle,CTD.depth_do_Apr,'Linear')+1000;
O2_down_A = (sig_A /1000)'.*CTD.O2_do_Apr;
bott_O2 = ((CTD.bottle_sig+1000)/1000) .* CTD.bottle_O2;
bott_O2_A = ((CTD.sig0_April_bottle+1000)/1000) .* CTD.O2_April_bottle;

s1 = scatter(vars.O2,prcdata.timeseries.depth,5,'MarkerFaceColor',[.6 .8 .8],'MarkerEdgeColor',[.6 .8 .8])
hold on
% scatter(O2_down,CTD.depth_do,10,'MarkerFaceColor',[0 .4 .4],'MarkerEdgeColor',[0 .4 .4])
% scatter(O2_down_A,CTD.depth_do_Apr,10,'MarkerFaceColor',[.8 .6 .6],'MarkerEdgeColor',[.8 .6 .6])

c = prcdata.timeseries.t < datenum(2016,03,09);
s2 = scatter(vars.O2(c),prcdata.timeseries.depth(c),10,'MarkerFaceColor',[0 .4 .4],'MarkerEdgeColor',[0 .4 .4])
c = prcdata.timeseries.t > datenum(2016,04,04);
s3 = scatter(vars.O2(c),prcdata.timeseries.depth(c),10,'MarkerFaceColor',[.8 .6 .6],'MarkerEdgeColor',[.8 .6 .6])

scatter(bott_O2,CTD.bottle_P,60,'MarkerFaceColor','k','MarkerEdgeColor','k');
s4 = scatter(bott_O2,CTD.bottle_P,30,'MarkerFaceColor',[0 .4 .4],'MarkerEdgeColor','k');
scatter(bott_O2_A,CTD.P_April_bottle,60,'MarkerFaceColor','k','MarkerEdgeColor','k');
s5 = scatter(bott_O2_A,CTD.P_April_bottle,30,'MarkerFaceColor',[.8 .6 .6],'MarkerEdgeColor','k');

set(gca,'YDir','Reverse','LineWidth',2,'Box','On','FontSize',16,'YLim',[-5 1020],'XLim',[170 285]);
grid on;
xlabel('\itc\rm(O_2) [mmol m^{-3}]');
ylabel('Depth [m]');

leg = legend([s1 s2 s3 s4 s5],'Glider 7 March - April 5','Glider 7 - 8 March','Glider 4 - 5 April','Ship 7 March','Ship 16 April');
set(leg,'Location','SouthEast','Box','On');

% Create axes
axes('Parent',gcf,...
    'Position',[0.570340909090909 0.16 0.406709872159091 0.788425925925926]);

scatter(prcdata.timeseries.pH_recalibrated,prcdata.timeseries.depth,5,'MarkerFaceColor',[.6 .8 .8],'MarkerEdgeColor',[.6 .8 .8])
hold on

c = prcdata.timeseries.t < datenum(2016,03,9)
s2 = scatter(prcdata.timeseries.pH_recalibrated(c),prcdata.timeseries.depth(c),10,'MarkerFaceColor',[0 .4 .4],'MarkerEdgeColor',[0 .4 .4])
c = prcdata.timeseries.t > datenum(2016,04,04);
s3 = scatter(prcdata.timeseries.pH_recalibrated(c),prcdata.timeseries.depth(c),10,'MarkerFaceColor',[.8 .6 .6],'MarkerEdgeColor',[.8 .6 .6])

scatter(CTD.pH,CTD.pH_pCO2_depth,60,'MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(CTD.pH,CTD.pH_pCO2_depth,30,'MarkerFaceColor',[0 .4 .4],'MarkerEdgeColor','k');
scatter(CTD.pH_A,CTD.pH_pCO2_depth_A,60,'MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(CTD.pH_A,CTD.pH_pCO2_depth_A,30,'MarkerFaceColor',[.8 .6 .6],'MarkerEdgeColor','k');

set(gca,'YDir','Reverse','LineWidth',2,'Box','On','FontSize',16,'YLim',[-5 1020],'XLim',[8.03 8.19],'YTickLabels','');
grid on;
xlabel('pH_T');

annotation(gcf,'textbox',...
    [0.568684895833334 0.941901234567902 0.001 0.001],'String','(b)',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.133463541666667 0.939329218106996 0.000999999999999973 0.001],...
    'String',{'(a)'},...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','off');

print(gcf, '-dpng','-r400', [options.plot_dir,'Fig_TScal'])

