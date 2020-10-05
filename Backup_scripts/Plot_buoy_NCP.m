
load([options.directory,'/data/BOUSSOLE.mat']);  
load([options.directory,'/data/3hrbins'],'bins')

%% plot DIC and morning maxim between March 16th and March 26th plus temperature and alkalinity


figure('units','normalized','position',[.01 .3 .9 .75]) 

% Create axes

axes('Parent',gcf,...
    'Position',[0.13 0.539259259259259 0.775 0.441481481481481]);
[handle] = add_bathymetry_to_plot([datenum(2016,03,20) datenum(2016,03,25) datenum(2016,03,25) datenum(2016,03,20)],[1 1 0 0]);
hold on 
plot(BOUSSOLE.time_CSYS,BOUSSOLE.sal_CSYS,'Color','k','LineWidth',5)    
a1 = plot(BOUSSOLE.time_CSYS,BOUSSOLE.sal_CSYS,'Color',[0.85 0.325 0.098],'LineWidth',3)        
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'XGrid','On','FontSize',24,'LineWidth',3,'XTickLabels','','Box','On','YTick','','XTick',[datenum(2016,03,10,00,00,00):1:datenum(2016,04,03,00,00,00)])
ylim([38.2 38.5])

pos = get(gca,'Position');
pos(1) = 0.12; pos(3) = 0;
axes('Position',pos)
ylim([38.2 38.5])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0.85 0.325 0.098])
ylabel('Salinity');

axes('Parent',gcf,...
    'Position',[0.13 0.539259259259259 0.775 0.441481481481481]);
plot(BOUSSOLE.time_3m,BOUSSOLE.temp_3m,'Color','k','LineWidth',5)  
hold on 
a2 = plot(BOUSSOLE.time_3m,BOUSSOLE.temp_3m,'Color',[0.929 0.694 0.125],'LineWidth',3)    
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'Visible','Off')
ylim([13.1 14.5])

pos = get(gca,'Position');
pos(1) = 0.95; pos(3) = 0;
axes('Position',pos)
ylim([13.1 14.5])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0.929 0.694 0.125])
y = ylabel('SST [^\circC]');
ylp = get(y, 'Position');
ext=get(y,'Extent');
set(y, 'Rotation',270, 'Position',ylp+[ext(3)+30 0 0]);

axes('Parent',gcf,...
    'Position',[0.13 0.539259259259259 0.775 0.441481481481481]);
[handle] = add_bathymetry_to_plot([datenum(2016,03,20) datenum(2016,03,25) datenum(2016,03,25) datenum(2016,03,20)],[1 1 0 0]);

set(gca,'Visible','Off')
set(handle,'FaceAlpha',0.05,'LineStyle','None')
hold on
[handle] = add_bathymetry_to_plot([datenum(2016,03,29) datenum(2016,04,02) datenum(2016,04,02) datenum(2016,03,29)],[1 1 0 0]);
set(handle,'FaceAlpha',0.05,'LineStyle','None')




% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.423333333333333]);

hold on
plot(BOUSSOLE.time_CSYS,BOUSSOLE.DIC_CSYS_Snorm_mmolm3,'LineWidth',5,'Color','k')    
plot(BOUSSOLE.time_CSYS,BOUSSOLE.DIC_CSYS_Snorm_mmolm3,'LineWidth',3,'Color',[0 0.447 0.741])
%     plot([bins.Bt],[bins.BDICnormS],'LineWidth',3,'Color','k','LineStyle','--')  
scatter([bins.Bt],[bins.BDICnormS],160,'k','filled')   
scatter([bins.Bt],[bins.BDICnormS],70,'MarkerFaceColor',[0 0.447 0.741])    
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
ylim([2335 2395])

t = [bins.Bt];
check = t > datenum(2016,03,20) & t < datenum(2016,03,25);
plot(t(check),NCP_DIC.M2025.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_DIC.M2025.fit(t(check)),'LineWidth',3,'Color',[0.2 0.8 0.4]);

check = t > datenum(2016,03,29) & t < datenum(2016,04,02);
plot(t(check),NCP_DIC.M2901.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_DIC.M2901.fit(t(check)),'LineWidth',3,'Color',[0.2 0.8 0.4]);

set(gca,'XGrid','On','FontSize',24,'LineWidth',3,'XTickLabels','','Box','On','YTick','','XTick',[datenum(2016,03,10,00,00,00):1:datenum(2016,04,03,00,00,00)])

pos = get(gca,'Position');
pos(1) = 0.12; pos(3) = 0;
axes('Position',pos)
ylim([2335 2395])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0 0.447 0.741])
ylabel('\itc\rm(DIC) [mmol m^{-3}]');

% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.423333333333333]);

BUOY.O2 = BOUSSOLE.O2_raw_10m_Liliane;
check = isfinite(BOUSSOLE.time_CSYS) & isfinite(BOUSSOLE.rho_CSYS);
BUOY.dens = interp1(BOUSSOLE.time_CSYS(check),BOUSSOLE.rho_CSYS(check),BOUSSOLE.O2_raw_10m_Liliane_date,'Linear');
BUOY.O2 = (BUOY.dens/1000)' .* BUOY.O2;
BUOY.O2_date = BOUSSOLE.O2_raw_10m_Liliane_date;

hold on
plot(BUOY.O2_date,BUOY.O2,'LineWidth',8,'Color','k')    
plot(BUOY.O2_date,BUOY.O2,'LineWidth',3,'Color',[1 1 1])
%     plot([bins.Bt],[bins.BDICnormS],'LineWidth',3,'Color','k','LineStyle','--')  
scatter([bins.Bt],[bins.O2],160,'k','filled')   
scatter([bins.Bt],[bins.O2],70,'MarkerFaceColor',[1 1 1])    
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'Visible','Off')
ylim([220 290])

check = t > datenum(2016,03,20) & t < datenum(2016,03,25);
plot(t(check),NCP_O2.M2025.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_O2.M2025.fit(t(check)),'LineWidth',3,'Color',[1 0 .6]);

check = t > datenum(2016,03,29) & t < datenum(2016,04,02);
plot(t(check),NCP_O2.M2901.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_O2.M2901.fit(t(check)),'LineWidth',3,'Color',[1 0 .6]);


axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.423333333333333]);
[handle] = add_bathymetry_to_plot([datenum(2016,03,20) datenum(2016,03,25) datenum(2016,03,25) datenum(2016,03,20)],[1 1 0 0]);
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'Visible','Off')
set(handle,'FaceAlpha',0.05,'LineStyle','None')
hold on
[handle] = add_bathymetry_to_plot([datenum(2016,03,29) datenum(2016,04,02) datenum(2016,04,02) datenum(2016,03,29)],[1 1 0 0]);
set(handle,'FaceAlpha',0.05,'LineStyle','None')


axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.423333333333333]);
pos = get(gca,'Position');
pos(1) = 0.95; pos(3) = 0;
set(gca,'FontSize',24,'LineWidth',5,'YColor',[0 0 0],'Position',pos)
ylim([220 290])
y = ylabel('\itc\rm(O_2) [mmol m^{-3}]');
ylp = get(y, 'Position');
ext=get(y,'Extent');
set(y, 'Rotation',270, 'Position',ylp+[ext(3) 0 0]);


axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.423333333333333]);
pos = get(gca,'Position');
pos(1) = 0.95; pos(3) = 0;
set(gca,'FontSize',24,'LineWidth',2,'YColor',[1 1 1],'Position',pos,'YTickLabels','')
ylim([220 290])

% Create axes
axes12 = axes('Parent',gcf,...
    'Position',[0.13 0.11 0.775 0.00296296296296296]);
set(gca,'FontSize',24,'LineWidth',2')
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
datetick('x','dd/mm','KeepLimits');
xlabel('Date')

set(gcf,'Color','W')


h=gcf;
set(h,'position',[.01 .3 .9 .75])
print(gcf, '-dtiff','-r400',[options.plot_dir,'DIC_SST_AT.tiff'])   
