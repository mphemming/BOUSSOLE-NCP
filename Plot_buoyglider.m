
%% Load data

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

%% prepare data for figure


for dive = 1:147
    closeness(dive).vectorBOUSSOLE = ones(size(prcdata.hydrography(dive).time))*closeness(dive).BOUSSOLEdist;
end

close1 = [closeness.vectorBOUSSOLE];
P = [prcdata.timeseries.P];
t = [prcdata.timeseries.t];
S = [prcdata.timeseries.S];
T = [prcdata.timeseries.T];
sig0 = [prcdata.timeseries.sigma0];
Sc700 = [prcdata.timeseries.Scatter_700];
O2 = [prcdata.timeseries.O2]; O2 = (sig0/1000) .* O2; % to mmol m-3
pH = [prcdata.timeseries.pH_re];

% select depth and close to BUSSOLE data

checks = (close1 <= 5 & P > 5 & P < 15);

tc = t(checks);
Pc = P(checks);
Sc = S(checks);
Tc = T(checks);
Sc700c = Sc700(checks);
O2c = O2(checks);
pHc = pH(checks);
sig0c = sig0(checks);

% bin data in time 

% daily bins

xrange = datenum(2016,03,06,00,00,00):1:datenum(2016,04,06,00,00,00);

for xr = 1:length(xrange)

    xtime = xrange(xr);
    
    checks = tc > xrange(xr) & tc < xrange(xr)+1;

    bins.t(xr) = xtime+0.5;
    bins.P(xr) = nanmean(Pc(checks));
    bins.S(xr) = nanmean(Sc(checks));
    bins.T(xr) = nanmean(Tc(checks));
    bins.Sc700(xr) = nanmean(Sc700c(checks));
    bins.O2(xr) = nanmean(O2c(checks));
    bins.pH(xr) = nanmean(pHc(checks));
    bins.sig0(xr) = nanmean(sig0c(checks));
    
    bins.P_std(xr) = nanstd(Pc(checks));
    bins.S_std(xr) = nanstd(Sc(checks));
    bins.T_std(xr) = nanstd(Tc(checks));
    bins.Sc700_std(xr) = nanstd(Sc700c(checks));
    bins.O2_std(xr) = nanstd(O2c(checks));
    bins.pH_std(xr) = nanstd(pHc(checks));   

% bin BOUSSOLE data in same way

    checks = BOUSSOLE.time_3m > xrange(xr)-(1/8) & BOUSSOLE.time_3m < xrange(xr)+(1/8); 

    bins.fCO213(xr) = nanmedian(BOUSSOLE.fCO213_3m(checks));    
    bins.fCO213_std(xr) = nanstd(BOUSSOLE.fCO213_3m(checks));        
    
    
end

% save('bins_near_boussole','bins')

%% create figure

 figure('units','normalized','position',[0 0 1 .9]);

axes4 = axes('Parent',gcf,'Position',[0.13 0.636484687083888 0.635432098765432 0.288515312916113]);
p4 = area(METEO.time_wind_speed,METEO.wind_speed,'FaceColor',[0.494117647409439 0.494117647409439 0.494117647409439],...
    'EdgeColor','none')

hold on

xlim([datenum(2016,03,07,00,00,00) datenum(2016,04,05,00,00,00)]);
set(axes4,'FontSize',22,'XTickLabel',[''],'LineWidth',2,'XGrid','On');
set(axes4,'XTick',[datenum(2016,03,13,00,00,00) datenum(2016,03,20,00,00,00) datenum(2016,03,27,00,00,00) datenum(2016,04,03,00,00,00)],...
    'XTickLabel','');
ylabel('[m s^{-1}]')

axes1 = axes('Parent',gcf,'Position',[0.13 0.636484687083888 0.635432098765432 0.288515312916113]);


p3 = plot(ones(size(320:1:450))*datenum(2016,03,18,18,00,00),320:1:450,...
    'Color',[0 0.498039215803146 0],...
    'LineStyle',':','LineWidth',6)

hold on

xpCO2.time = [datenum(2016,03,04),datenum(2016,03,11),datenum(2016,03,18),datenum(2016,03,26),datenum(2016,04,1),datenum(2016,04,15)];
xpCO2.val = [408.23,407.24, 408.68, 409.12, 406.31,408.08];
xpCO2.val_int  = interp1(xpCO2.time,xpCO2.val,METEO.date,'Linear'); 
xpCO2.val_int_pCO2  = xpCO2.val_int .*smooth((METEO.sea_level_pressure/100000),5)';

p1 = plot(BOUSSOLE.time_3m,BOUSSOLE.fCO213_3m,'LineWidth',4,'Color',[0 0.447058826684952 0.74117648601532])
xlim([datenum(2016,03,07,00,00,00) datenum(2016,04,05,00,00,00)]);
ylim([325 450])
p2 = plot(METEO.date,xpCO2.val_int_pCO2 ,'LineStyle',':','LineWidth',2,'Color',[0 0.447058826684952 0.74117648601532]);

set(axes1,'FontSize',14,'XTickLabel',[''],'LineWidth',2,'Visible','Off');

axes5 = axes('Parent',gcf,...
    'Position',[0.0673809523809524 0.632490013315579 0.00198412698412698 0.292509986684423]);

ylim([325 450])
set(axes5,'FontSize',22,'XTickLabel',[''],'LineWidth',2,'YColor',[0 0.447058826684952 0.74117648601532]);
ylabel('[µatm]')

axes2 = axes('Parent',gcf,'Position',[0.13 0.636484687083888 0.635432098765432 0.288515312916113]);

p2 = plot(BOUSSOLE.time_3m,BOUSSOLE.temp_3m,'LineWidth',4,'Color','r')
xlim([datenum(2016,03,07,00,00,00) datenum(2016,04,05,00,00,00)]);
ylim([13 14.5])
hold on;

l = legend([p1 p2],'f(CO_2)_{13}','SST','Location','NorthWest','Orientation','Horizontal');
lpos = get(l,'Position');
lpos(1) = 0.17;
set(l,'FontSize',24,'Box','Off','Position',lpos);


set(axes2,'Visible','Off')

axes3 = axes('Parent',gcf,'Position',[0.831901408450704 0.637816245006658 0.00176056338028174 0.285852197070575]);
ylim([13 14.5])
set(axes3,'FontSize',22,'YColor','r','LineWidth',2)
y = ylabel('[^\circC]','Rotation',270)
ypos = get(y,'Position')
ypos(1) = 5;
set(y,'Position',ypos);

% axes('Position',[0.831901408450704 0.637816245006658 0.00176056338028174 0.285852197070575])
% 
% l2 = legend([p3 p4],'Bloom','Wind Speed','Location','NorthEast','Orientation','Horizontal');
% lpos = get(l2,'Position');
% lpos(1) = 0.17;
% set(l2,'FontSize',24,'Box','Off','Position',lpos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



axes6 = axes('Parent',gcf,'Position',[0.13 0.11 0.635432098765432 0.499853528628495]);

p3 = plot(ones(size(-1000:1:1000))*datenum(2016,03,18,18,00,00),-1000:1:1000,...
    'Color',[0 0.498039215803146 0],...
    'LineStyle',':','LineWidth',6)
xlim([datenum(2016,03,7,00,00,00) datenum(2016,04,05,00,00,00)]);

set(axes6,'LineWidth',2,'FontSize',22,'Color',[0.972549019607843 0.972549019607843 0.972549019607843],'YTick',...
    [1.5 3 4.5 6],'XGrid','On');
hold on

pSc = plot(bins.t(~isnan(bins.t) & ~isnan(bins.Sc700)),bins.Sc700(~isnan(bins.t) & ~isnan(bins.Sc700))*10000,'k','LineWidth',3)
scatter(bins.t,bins.Sc700*10000,100,'filled','k')
ylim([1.5 6])
ylabel('x 10^{-4} [700 nm]');

axes7 = axes('Parent',gcf,'Position',[0.13 0.11 0.635432098765432 0.499853528628495]);

pS = plot(bins.t(~isnan(bins.t) & ~isnan(bins.S)),bins.S(~isnan(bins.t) & ~isnan(bins.S)),'Color',[0.8 0.4 0],'LineWidth',3)
hold on
scatter(bins.t,bins.S,100,'filled','MarkerFaceColor',[0.8 0.4 0.0])
ylim([38.22 38.45])


xlim([datenum(2016,03,7,00,00,00) datenum(2016,04,05,00,00,00)]);

set(axes7,'LineWidth',2,'FontSize',14,'Visible','Off');


axes8 = axes('Parent',gcf,'Position',[0.0734126984126984 0.11 0.00297619047619048 0.498521970705726]);
ylim([38.22 38.45])
set(axes8,'LineWidth',2,'FontSize',22,'YColor',[0.8 0.4 0]);



axes9 = axes('Parent',gcf,'Position',[0.13 0.11 0.635432098765432 0.499853528628495]);

pO2 = plot(bins.t(~isnan(bins.t) & ~isnan(bins.S)),bins.O2(~isnan(bins.t) & ~isnan(bins.S)),'Color',[0.4 0.6 0.6],'LineWidth',3)
hold on
scatter(bins.t,bins.O2,100,'filled','MarkerFaceColor',[0.4 0.6 0.6])
ylim([230 280])

O2satt = o2satv2b(bins.S,bins.T);
O2satt = O2satt .* (bins.sig0/1000);
O2satt2 = inpaint_nans(O2satt',4);
O2satt2(end-2:end) = NaN;
O2satt2(1:3) = NaN;



s = scatter(bins.t,O2satt,100,'filled','Sq','MarkerFaceColor',[0.4 0.6 0.6])
hold on
plot(bins.t,O2satt2,':','LineWidth',3,'Color',[0.4 0.6 0.6]);


xlim([datenum(2016,03,7,00,00,00) datenum(2016,04,05,00,00,00)]);

set(axes9,'LineWidth',2,'FontSize',14,'Visible','Off');

axes10 = axes('Parent',gcf,'Position',[0.832460317460318 0.11 0.00198412698412698 0.498521970705726]);
ylim([230 270])
set(axes10,'LineWidth',2,'FontSize',22,'YColor',[0.4 0.6 0.6],'YTick',[230 240 250 260 270]);
y = ylabel('[mmol m^{-3}]','Rotation',270)
ypos = get(y,'Position');
ypos(1) = 5;
set(y,'Position',ypos)



axes11 = axes('Parent',gcf,'Position',[0.13 0.11 0.635432098765432 0.499853528628495]);

c = isfinite(bins.t) & isfinite(bins.S) & isfinite(bins.pH);

ppH = plot(bins.t(c),bins.pH(c),'Color',[0.8 0.2 0.4],'LineWidth',3)
hold on
scatter(bins.t,bins.pH,100,'filled','MarkerFaceColor',[0.8 0.2 0.4])
ylim([8.04 8.17])


xlim([datenum(2016,03,7,00,00,00) datenum(2016,04,05,00,00,00)]);

set(axes11,'LineWidth',2,'FontSize',14,'Visible','Off');


axes12 = axes('Parent',gcf,'Position',[0.90503968253968 0.11051930758988 0.00198412698412698 0.496671105193077]);
ylim([8.04 8.24])
set(axes12,'LineWidth',2,'FontSize',22,'YColor',[0.8 0.2 0.4],'YTick',[8.04 8.08 8.12 8.16 8.2 8.24]);

datetick(axes6,'x','dd/mm','KeepLimits')


l = legend([pSc pS pO2 ppH],'Back scatter','Salinity','\itc\rm(O_2)','pH','Location','NorthWest')
set(l,'FontSize',22)
set(gcf,'Color','w')
set(l,'Color','w','Box','Off')
lpos = get(l,'Position');
lpos(2) = 0.445;
lpos(1) = 0.15;
set(l,'Position',lpos)

xlabel(axes6,'Date [2016]','FontSize',22)

annotation(gcf,'textbox',...
    [0.251520833333333 0.436213991769547 0.0802499999999993 0.0771604938271604],...
    'String','\itc_{\rmsat}\rm(O_2)',...
    'LineStyle','none',...
    'FontSize',22,...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.445010416666666 0.848251028806585 0.0802499999999992 0.0771604938271604],...
    'String','\chi\rmatm(CO_2)',...
    'LineStyle','none',...
    'FontSize',22,...
    'FitBoxToText','off');

annotation(gcf,'line',[0.414583333333333 0.441666666666667],...
    [0.890946502057613 0.890946502057613],...
    'Color',[0 0.447058823529412 0.741176470588235],...
    'LineWidth',4,...
    'LineStyle',':');

annotation(gcf,'line',[0.223697916666666 0.250781249999999],...
    [0.476851851851853 0.476851851851852],...
    'Color',[0.4 0.6 0.6],...
    'LineWidth',4,...
    'LineStyle',':');

% print(gcf, '-dpng','-r400', [options.plot_dir,'Fig_buoyglider'])
