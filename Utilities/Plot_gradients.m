%% plot for paper combining two figures O2 gradient and speed

% plot O2 x and y gradients

for n = 1:27
    oxy_x(n) = nanmean(O2_adv(n).oxy_x)*1000;
    oxy_y(n) = nanmean(O2_adv(n).oxy_y*1000);
    grad_x(n) = [planes_loop(n+8).across_O2.gradient_x_across];
    grad_y(n) = [planes_loop(n+8).across_O2.gradient_y_across];
    U(n) = nanmean(O2_adv(n).U);
    V(n) = nanmean(O2_adv(n).V);
    U_surf(n) = O2_adv(n).U(2);    
    V_surf(n) = O2_adv(n).V(2);    
    AVISO_U(n) = planes_loop(n+8).means.GVsat_U;
    AVISO_V(n) = planes_loop(n+8).means.GVsat_V;
    Oscar_U(n) = planes_loop(n+8).means.GVsat_Oscar_U;
    Oscar_V(n) = planes_loop(n+8).means.GVsat_Oscar_V;    
end
% NOTES: Need to replicate across transects using gridded data
%%
figure('units','normalized','position',[0 0 .6 .9]);

%----------------------------------------------------------------------
% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.800677966101695 0.775 0.170258247890075]);

plot(datenum(2016,03,options.dayrange),oxy_x,'-o','LineWidth',3)
% hold on
% plot(datenum(2016,03,options.dayrange),grad_x,'-o','LineWidth',3)
set(gca,'LineWidth',3,'FontSize',16,'XTickLabels','')
annotation(gcf,'textbox',...
    [0.915930555555557 0.86 0.0667083333333336 0.0606995884773663],...
    'String',{'\delta/\deltax'},...
    'FitBoxToText','off','LineStyle','None','FontSize',22);
ylabel('         \itc\rm(O_2)  \newline \newline [mmol m^{-3} km^{-1}]','FontSize',14)
xlim([datenum(2016,03,8) datenum(2016,04,05)])
ylim([-0.6 0.6]); add_zero

%----------------------------------------------------------------------
% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.628008474576271 0.775 0.16597300690521]);

plot(datenum(2016,03,options.dayrange),oxy_y,'-o','LineWidth',3)
% hold on
% plot(datenum(2016,03,options.dayrange),grad_y,'-o','LineWidth',3)
set(gca,'LineWidth',3,'FontSize',18,'XTickLabels','')
annotation(gcf,'textbox',...
    [0.915930555555557 0.68 0.0667083333333336 0.0606995884773663],...
    'String',{'\delta/\deltay'},...
    'FitBoxToText','off','LineStyle','None','FontSize',22);
y = ylabel('         \itc\rm(O_2)  \newline \newline [mmol m^{-3} km^{-1}]','FontSize',14)
xlim([datenum(2016,03,8) datenum(2016,04,05)])
ylim([-0.6 0.6])
add_zero

%----------------------------------------------------------------------
% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.455338983050847 0.775 0.164774185673433]);

plot(datenum(2016,03,options.dayrange),U_surf*100,'-o','linewidth',3)
hold on
plot(datenum(2016,03,options.dayrange),AVISO_U*100,'-o','linewidth',3)
plot(datenum(2016,03,options.dayrange),Oscar_U*100,'-o','linewidth',3)
box on;
set(gca,'LineWidth',3,'FontSize',18,'XTickLabels','')
ylabel('[cm s^{-1}]','FontSize',14)
datetick('x','KeepLimits')
xlim([datenum(2016,03,08,00,00,00) datenum(2016,04,05,00,00,00)]);
ylim([-20 20])
add_zero
annotation(gcf,'textbox',...
    [0.915930555555557 0.51 0.0667083333333336 0.0606995884773663],...
    'String',{'\itU_{\rm{abs}}'},...
    'FitBoxToText','off','LineStyle','None','FontSize',22);
set(gca,'YTick',[-10 0 10]);

%----------------------------------------------------------------------
% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.282669491525424 0.775 0.167690590779103]);

plot(datenum(2016,03,options.dayrange),V_surf*100,'-o','linewidth',3)
hold on
plot(datenum(2016,03,options.dayrange),AVISO_V*100,'-o','linewidth',3)
plot(datenum(2016,03,options.dayrange),Oscar_V*100,'-o','linewidth',3)
box on;
set(gca,'LineWidth',3,'FontSize',18,'XTickLabels','')
ylabel('[cm s^{-1}]','FontSize',14)
xlim([datenum(2016,03,08,00,00,00) datenum(2016,04,05,00,00,00)]);
ylim([-20 20])
add_zero
annotation(gcf,'textbox',...
    [0.915930555555557 0.34 0.0667083333333336 0.0606995884773663],...
    'String',{'\itV_{\rm{abs}}'},...
    'FitBoxToText','off','LineStyle','None','FontSize',22);
set(gca,'YTick',[-10 0 10]);

%----------------------------------------------------------------------
% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.166491769547325]);

plot(datenum(2016,03,options.dayrange),sqrt(U_surf.^2 +V_surf.^2)*100,'-o','linewidth',3)
hold on
plot(datenum(2016,03,options.dayrange),sqrt(AVISO_U.^2 +AVISO_V.^2)*100,'-o','linewidth',3)
plot(datenum(2016,03,options.dayrange),sqrt(Oscar_U.^2 +Oscar_V.^2)*100,'-o','linewidth',3)

set(gca,'FontSize',18,'LineWidth',3,'XTick',[datenum(2016,03,12) datenum(2016,03,19) datenum(2016,03,26) datenum(2016,04,02)])
box on;
ylabel('[cm s^{-1}]','FontSize',14)
datetick('x','KeepLimits')
xlim([datenum(2016,03,08,00,00,00) datenum(2016,04,05,00,00,00)]);
xlabel('Date')
datetick('x','dd/mm','KeepLimits')
ylim([0 30])
set(gca,'YTick',[0 10 20]);

l = legend('Glider_{abs,0-10m}','AVISO_{abs,surf}','OSCAR_{abs,15m}','Location','NorthWest','Orientation','Horizontal')
set(l,'Box','Off')
annotation(gcf,'textbox',...
    [0.915930555555557 0.17 0.0667083333333336 0.0606995884773663],...
    'String',{'Speed'},...
    'FitBoxToText','off','LineStyle','None','FontSize',22);


print(gcf, '-dtiff','-r400', [options.plot_dir,'GVsat_compare_combined.tiff'])



xlabel('Date')