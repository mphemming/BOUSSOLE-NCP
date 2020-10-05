
%% plot DIC and morning maxim between March 16th and March 26th plus temperature and alkalinity

figure('units','normalized','position',[.01 .3 .9 .75]) 
hold on 
plot(BOUSSOLE.time_CSYS,BOUSSOLE.sal_CSYS,'Color','k','LineWidth',4)    
a1 = plot(BOUSSOLE.time_CSYS,BOUSSOLE.sal_CSYS,'Color',[0.85 0.325 0.098],'LineWidth',2)        
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
datetick('x','dd/mm','KeepLimits')    
set(gca,'Visible','Off')
ylim([38 38.6])


pos = get(gca,'Position');
pos(1) = 0.2;
pos(3) = 0.6;
set(gca,'XGrid','On','Position',pos,'FontSize',24,'LineWidth',3)
datetick('x','dd/mm','KeepLimits')  


axes('Position',[0.854197530864197 0.11 0.00231481481481521 0.815])
ylim([38 38.6])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0.85 0.325 0.098])

y = ylabel('S')
ypos = get(y,'Position'); 
ypos(1) = 12;
set(y,'Position',ypos);

axes('Position',pos)
hold on
plot(BOUSSOLE.time_3m,BOUSSOLE.temp_3m,'Color','k','LineWidth',4)  
a2 = plot(BOUSSOLE.time_3m,BOUSSOLE.temp_3m,'Color',[0.929 0.694 0.125],'LineWidth',2)    
   
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
datetick('x','dd/mm','KeepLimits')  
set(gca,'Visible','Off')
ylim([13.1 14.5])

axes('Position',[0.1 0.11 0.00231481481481521 0.815])
ylim([13.1 14.5])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0.929 0.694 0.125])

y = ylabel('SST [^\circC]')


%    axes('Position',pos)   
%    plot([bins.Bt],[bins.BALK],'LineWidth',3,'Color','k','LineStyle','--')
%    hold on
%    scatter([bins.Bt],[bins.BALK],100,[0.85 0.325 0.098],'filled','MarkerEdgeColor','k')    
%    xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
%    datetick('x','KeepLimits')    
%    set(gca,'Visible','Off')
%    ylim([2566 2586])
%    


axes('Position',pos)
hold on
plot(BOUSSOLE.time_CSYS,BOUSSOLE.DIC_CSYS_Snorm_mmolm3,'LineWidth',4,'Color','k')    
plot(BOUSSOLE.time_CSYS,BOUSSOLE.DIC_CSYS_Snorm_mmolm3,'LineWidth',2,'Color',[0 0.447 0.741])
%     plot([bins.Bt],[bins.BDICnormS],'LineWidth',3,'Color','k','LineStyle','--')  
scatter([bins.Bt],[bins.BDICnormS],160,'k','filled')   
scatter([bins.Bt],[bins.BDICnormS],70,'MarkerFaceColor',[0 0.447 0.741])    
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'Visible','Off')
ylim([2335 2390])

t = [bins.Bt];
check = t > datenum(2016,03,20) & t < datenum(2016,03,25);
plot(t(check),NCP_DIC.M2025.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_DIC.M2025.fit(t(check)),'LineWidth',3,'Color',[0.2 0.8 0.4]);

check = t > datenum(2016,03,29) & t < datenum(2016,04,02);
plot(t(check),NCP_DIC.M2901.fit(t(check)),'LineWidth',8,'Color','k');
plot(t(check),NCP_DIC.M2901.fit(t(check)),'LineWidth',3,'Color',[0.2 0.8 0.4]);

BUOY.O2 = BOUSSOLE.O2_raw_10m_Liliane;
check = isfinite(BOUSSOLE.time_CSYS) & isfinite(BOUSSOLE.rho_CSYS);
BUOY.dens = interp1(BOUSSOLE.time_CSYS(check),BOUSSOLE.rho_CSYS(check),BOUSSOLE.O2_raw_10m_Liliane_date,'Linear');
BUOY.O2 = (BUOY.dens/1000)' .* BUOY.O2;
BUOY.O2_date = BOUSSOLE.O2_raw_10m_Liliane_date;

axes('Position',pos)
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


%     
%     plot(t(fits.check12),fits.B12_fit(t(fits.check12)),'LineWidth',8,'Color','k');
%     plot(t(fits.check12),fits.B12_fit(t(fits.check12)),'LineWidth',3,'Color',[0.2 0.8 0.4]);        
%     
%     plot(t(fits.check17),fits.B17_fit(t(fits.check17)),'LineWidth',8,'Color','k');
%     plot(t(fits.check17),fits.B17_fit(t(fits.check17)),'LineWidth',3,'Color',[0.2 0.8 0.4]);      
%     
axes('Position',[0.19 0.11 0.00231481481481521 0.815])
ylim([2335 2390])
set(gca,'FontSize',24,'LineWidth',3,'YColor',[0 0.447 0.741])
ylabel('\itc\rm(DIC) [mmol m^{-3}]')

axes('Position',[0.92 0.11 0.00231481481481521 0.815])
set(gca,'FontSize',24,'LineWidth',10,'YColor','k');
ylim([220 300])
axes('Position',[0.92 0.11 0.00231481481481521 0.815])
set(gca,'FontSize',24,'LineWidth',5,'YColor','w','YTickLabels','');
ylim([220 300])
y = ylabel('\itc\rm(O_2) [mmol m^{-3}]','Color','k');
ypos = get(y,'Position'); 
ypos(1) = 20;
set(y,'Position',ypos);

axes( 'Position',[0.203 0.082962962962963 0.595796296296296 0.0074074074074074]);
xlim([datenum(2016,03,10,00,00,00) datenum(2016,04,03,00,00,00)]);
set(gca,'FontSize',24,'LineWidth',3)    
datetick('x','dd/mm','KeepLimits')    
x1 = xlabel('Date')
x1pos = get(x1,'Position');
x1pos(2) = -3;
x1pos(1) =datenum(2016,03,23);
set(x1,'Position',x1pos);
set(gca,'YTickLabel','')

set(gcf,'Color','W')


h=gcf;
set(h,'position',[.01 .3 .9 .75])
print(gcf, '-dtiff','-r400',[options.plot_dir,'DIC_SST_AT.tiff'])   
