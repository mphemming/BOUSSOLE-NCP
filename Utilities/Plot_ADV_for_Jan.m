
% script to plot advection to show Jan in meeting

%% Figure
% Advection: Could you include time series of u, v, dc/dx, dc/dy, for each of the depth bins (0 to 45/46 m)? 
% Maybe 4 plots with 10 different lines (one for each of the depth bins)?

% day = 19

% concatenate all the data
for n_day = options.dayrange-8
    for n_depth = 1:10
        data_ADV(n_day,n_depth).U = O2_adv(n_day).DU_abs(n_depth);
        data_ADV(n_day,n_depth).U = O2_adv(n_day).DU_errors;
        data_ADV(n_day,n_depth).V = O2_adv(n_day).DV_abs(n_depth); 
        if n_depth < numel(O2_adv(n_day).oxy_x) 
            data_ADV(n_day,n_depth).oxy_x = O2_adv(n_day).oxy_x(n_depth);
            data_ADV(n_day,n_depth).oxy_x_err = O2_adv(n_day).oxy_x_errors(n_depth);
            data_ADV(n_day,n_depth).oxy_y = O2_adv(n_day).oxy_y(n_depth);
            data_ADV(n_day,n_depth).oxy_y_err = O2_adv(n_day).oxy_y_errors(n_depth);
        else
            data_ADV(n_day,n_depth).oxy_x = NaN;
            data_ADV(n_day,n_depth).oxy_x_err = NaN;
            data_ADV(n_day,n_depth).oxy_y = NaN;
            data_ADV(n_day,n_depth).oxy_y_err = NaN;
        end
    end
end

figure('units','normalized','position',[0 0 1 .9]);

cm = colormap(cbrewer('qual','Set1', 10));

subplot(2,2,1)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).U],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\itU \rm[m s^{-1}]')


subplot(2,2,2)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).V],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\itV \rm[m s^{-1}]')

subplot(2,2,3)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).oxy_x],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\partial\itc\rm(O_2) / \partialx [mmol m^{-3} m^{-1}]')

subplot(2,2,4)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).oxy_y],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\partial\itc\rm(O_2) / \partialy [mmol m^{-3} m^{-1}]')

leg = legend('5m','10m','15m','20m','25m','30m','35m','40m','45m');

%% Figure - same but errors

figure('units','normalized','position',[0 0 1 .9]);

cm = colormap(cbrewer('qual','Set1', 10));

subplot(2,2,1)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).U],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\itU \rm[m s^{-1}]')


subplot(2,2,2)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).V],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\itV \rm[m s^{-1}]')

subplot(2,2,3)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).oxy_x],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\partial\itc\rm(O_2) / \partialx [mmol m^{-3} m^{-1}]')

subplot(2,2,4)

for n_depth = 2:10
    plot(datenum(2016,03,options.dayrange),[data_ADV(:,n_depth).oxy_y],'LineWidth',2,'Color',cm(n_depth,:));
    hold on
end
datetick
xlim([datenum(2016,03,7) datenum(2016,04,7)])
set(gca,'FontSize',16,'LineWidth',2,'XGrid','On');
ylabel('\partial\itc\rm(O_2) / \partialy [mmol m^{-3} m^{-1}]')

leg = legend('5m','10m','15m','20m','25m','30m','35m','40m','45m');

%% create figure to be opened for gradients

%% 18/03
figure('units','normalized','position',[0 0 1 .9]);
subplot(1,2,1);
plot(planes_loop(18).oxygen_planes.gridded_O2fit(2).vals, ...
    planes_loop(18).oxygen_planes.gridded_O2fit(2).gridded_plane_matrix(:,1:2), ...
    planes_loop(18).oxygen_planes.gridded_O2fit(2).O2_vals_fit)
xlabel('Longitude distance [km]')
ylabel('Lattude distance [km]')
zlabel('\itc\rm(O_2) [mmol m^{-3}]');
title('18/03 between 2.5-7.5m');
annotation(gcf,'textbox',...
    [0.02 0.884615384615384 0.110458333333333 0.0737179487179488],...
    'String',{['m_x =', num2str(planes_loop(18).oxygen_planes.gridded_O2fit(2).vals.p10)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.83 0.110458333333333 0.0737179487179488],...
    'String',{['m_y =', num2str(planes_loop(18).oxygen_planes.gridded_O2fit(2).vals.p01)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.78 0.110458333333333 0.0737179487179488],...
    'String',{['R^2 =', num2str(planes_loop(18).oxygen_planes.gridded_gof(2).vals.rsquare)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.73 0.110458333333333 0.0737179487179488],...
    'String',{['RMSE =', num2str(planes_loop(10).oxygen_planes.gridded_gof(2).vals.rmse)]},...
    'FitBoxToText','off','LineStyle','None');
set(gca,'Box','On','FontSize',16);

subplot(1,2,2);
scatter(planes_loop(18).oxygen_planes.O2_adv,planes_loop(18).oxygen_planes.P_adv,'filled')
add_l(2.5,1);
add_l(7.5,1);
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');


%% 28/03
figure('units','normalized','position',[0 0 1 .9]);
subplot(1,2,1);
plot(planes_loop(28).oxygen_planes.gridded_O2fit(8).vals, ...
    planes_loop(28).oxygen_planes.gridded_O2fit(8).gridded_plane_matrix(:,1:2), ...
    planes_loop(28).oxygen_planes.gridded_O2fit(8).O2_vals_fit)
xlabel('Longitude distance [km]')
ylabel('Lattude distance [km]')
zlabel('\itc\rm(O_2) [mmol m^{-3}]');
title('28/03 between 32.5-37.5m');
annotation(gcf,'textbox',...
    [0.02 0.884615384615384 0.110458333333333 0.0737179487179488],...
    'String',{['m_x =', num2str(planes_loop(28).oxygen_planes.gridded_O2fit(8).vals.p10)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.83 0.110458333333333 0.0737179487179488],...
    'String',{['m_y =', num2str(planes_loop(28).oxygen_planes.gridded_O2fit(8).vals.p01)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.78 0.110458333333333 0.0737179487179488],...
    'String',{['R^2 =', num2str(planes_loop(28).oxygen_planes.gridded_gof(8).vals.rsquare)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.73 0.110458333333333 0.0737179487179488],...
    'String',{['RMSE =', num2str(planes_loop(10).oxygen_planes.gridded_gof(8).vals.rmse)]},...
    'FitBoxToText','off','LineStyle','None');
set(gca,'Box','On','FontSize',16);

subplot(1,2,2);
scatter(planes_loop(28).oxygen_planes.O2_adv,planes_loop(28).oxygen_planes.P_adv,'filled')
add_l(32.5,1);
add_l(37.5,1);
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');


%% 29/03
figure('units','normalized','position',[0 0 1 .9]);
subplot(1,2,1);
plot(planes_loop(29).oxygen_planes.gridded_O2fit(2).vals, ...
    planes_loop(29).oxygen_planes.gridded_O2fit(2).gridded_plane_matrix(:,1:2), ...
    planes_loop(29).oxygen_planes.gridded_O2fit(2).O2_vals_fit)
xlabel('Longitude distance [km]')
ylabel('Lattude distance [km]')
zlabel('\itc\rm(O_2) [mmol m^{-3}]');
title('29/03 between 32.5-37.5m');
annotation(gcf,'textbox',...
    [0.02 0.884615384615384 0.110458333333333 0.0737179487179488],...
    'String',{['m_x =', num2str(planes_loop(29).oxygen_planes.gridded_O2fit(2).vals.p10)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.83 0.110458333333333 0.0737179487179488],...
    'String',{['m_y =', num2str(planes_loop(29).oxygen_planes.gridded_O2fit(2).vals.p01)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.78 0.110458333333333 0.0737179487179488],...
    'String',{['R^2 =', num2str(planes_loop(29).oxygen_planes.gridded_gof(2).vals.rsquare)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.73 0.110458333333333 0.0737179487179488],...
    'String',{['RMSE =', num2str(planes_loop(10).oxygen_planes.gridded_gof(2).vals.rmse)]},...
    'FitBoxToText','off','LineStyle','None');
set(gca,'Box','On','FontSize',16);

subplot(1,2,2);
scatter(planes_loop(29).oxygen_planes.O2_adv,planes_loop(29).oxygen_planes.P_adv,'filled')
add_l(32.5,1);
add_l(37.5,1);
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');



%% 03/04
figure('units','normalized','position',[0 0 1 .9]);
subplot(1,2,1);
plot(planes_loop(34).oxygen_planes.gridded_O2fit(5).vals, ...
    planes_loop(34).oxygen_planes.gridded_O2fit(5).gridded_plane_matrix(:,1:2), ...
    planes_loop(34).oxygen_planes.gridded_O2fit(5).O2_vals_fit)
xlabel('Longitude distance [km]')
ylabel('Lattude distance [km]')
zlabel('\itc\rm(O_2) [mmol m^{-3}]');
title('03/04 between 17.5 - 22.5 m');
annotation(gcf,'textbox',...
    [0.02 0.884615384615384 0.110458333333333 0.0737179487179488],...
    'String',{['m_x =', num2str(planes_loop(34).oxygen_planes.gridded_O2fit(5).vals.p10)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.83 0.110458333333333 0.0737179487179488],...
    'String',{['m_y =', num2str(planes_loop(34).oxygen_planes.gridded_O2fit(5).vals.p01)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.78 0.110458333333333 0.0737179487179488],...
    'String',{['R^2 =', num2str(planes_loop(34).oxygen_planes.gridded_gof(5).vals.rsquare)]},...
    'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',...
    [0.02 0.73 0.110458333333333 0.0737179487179488],...
    'String',{['RMSE =', num2str(planes_loop(10).oxygen_planes.gridded_gof(5).vals.rmse)]},...
    'FitBoxToText','off','LineStyle','None');
set(gca,'Box','On','FontSize',16);

subplot(1,2,2);
scatter(planes_loop(34).oxygen_planes.O2_adv,planes_loop(34).oxygen_planes.P_adv,'filled')
add_l(17.5,1);
add_l(22.5,1);
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');
