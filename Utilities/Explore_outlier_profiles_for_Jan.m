% Script to investigate utlier profiles on certain days when calculating
% advection

%% 18th March

check = vars.t > datenum(2016,03,14) & vars.t < datenum(2016,03,22) & vars.P <= 46;
check_low = vars.O2 < 230;

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(1,4,1);
scatter(vars.O2(check),vars.P(check),'filled')
hold on
scatter(vars.O2(check & check_low),vars.P(check & check_low),'filled')
ylim([4 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');

subplot(1,4,2);
scatter(vars.O2_sat(check),vars.P(check),'filled')
hold on
scatter(vars.O2_sat(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm_{sat}(O_2) [mmol m^{-3}]');

subplot(1,4,3);
scatter(vars.S(check),vars.P(check),'filled')
hold on
scatter(vars.S(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Salinity');

subplot(1,4,4);
scatter(vars.T(check),vars.P(check),'filled')
hold on
scatter(vars.T(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Temperature [^\circC]');
set(gcf,'Color','w')

% create long-lat figure

figure;
scatter(vars.lon(check),vars.lat(check),'filled')
hold on
scatter(vars.lon(check & check_low),vars.lat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Latitude [^\circ N]');
xlabel('Longitude [^\circ E]');
set(gcf,'Color','w')

% create O2 in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.O2(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc\rm(O_2) [mmol m^{-3}]');
datetick
set(gcf,'Color','w')

% create MLD in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.MLD(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.MLD(check & check_low),'filled')
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
datetick
set(gcf,'Color','w')

% create O2sat, T and S in time plot

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(3,1,1);
scatter(vars.t(check),vars.O2_sat(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2_sat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc_{sat}\rm(O_2) [mmol m^{-3}]');
datetick

subplot(3,1,2);
scatter(vars.t(check),vars.T(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.T(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Temp. [^\circC]');
datetick

subplot(3,1,3);
scatter(vars.t(check),vars.S(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.S(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Salinity');
datetick

set(gcf,'Color','w')

%% 28th March

check = vars.t > datenum(2016,03,24) & vars.t < datenum(2016,04,01) & vars.P <= 46;
check_low_1 = vars.O2 < 233 & vars.P > 28;
check_low_2 = vars.O2 < 243 & vars.P < 28;
check_low_3 = vars.O2 > 265;
check_low = check_low_1 + check_low_2 + check_low_3;
check_low(check_low > 1) = 1;

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(1,4,1);
scatter(vars.O2(check),vars.P(check),'filled')
hold on
scatter(vars.O2(check & check_low),vars.P(check & check_low),'filled')
ylim([4 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');

subplot(1,4,2);
scatter(vars.O2_sat(check),vars.P(check),'filled')
hold on
scatter(vars.O2_sat(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm_{sat}(O_2) [mmol m^{-3}]');

subplot(1,4,3);
scatter(vars.S(check),vars.P(check),'filled')
hold on
scatter(vars.S(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Salinity');

subplot(1,4,4);
scatter(vars.T(check),vars.P(check),'filled')
hold on
scatter(vars.T(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Temperature [^\circC]');
set(gcf,'Color','w')

% create long-lat figure

figure;
scatter(vars.lon(check),vars.lat(check),'filled')
hold on
scatter(vars.lon(check & check_low),vars.lat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Latitude [^\circ N]');
xlabel('Longitude [^\circ E]');
set(gcf,'Color','w')

% create O2 in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.O2(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc\rm(O_2) [mmol m^{-3}]');
datetick
set(gcf,'Color','w')

% create MLD in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.MLD(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.MLD(check & check_low),'filled')
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
datetick
set(gcf,'Color','w')


% create O2sat, T and S in time plot

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(3,1,1);
scatter(vars.t(check),vars.O2_sat(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2_sat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc_{sat}\rm(O_2) [mmol m^{-3}]');
datetick

subplot(3,1,2);
scatter(vars.t(check),vars.T(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.T(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Temp. [^\circC]');
datetick

subplot(3,1,3);
scatter(vars.t(check),vars.S(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.S(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Salinity');
datetick

set(gcf,'Color','w')

%% 3rd May

check = vars.t > datenum(2016,04,03)-4 & vars.t < datenum(2016,04,03)+4 & vars.P <= 46;
check_low_1 = vars.O2 < 245;
check_low_2 = vars.O2 > 270;
check_low = check_low_1 + check_low_2;
check_low(check_low > 1) = 1;

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(1,4,1);
scatter(vars.O2(check),vars.P(check),'filled')
hold on
scatter(vars.O2(check & check_low),vars.P(check & check_low),'filled')
ylim([4 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm(O_2) [mmol m^{-3}]');

subplot(1,4,2);
scatter(vars.O2_sat(check),vars.P(check),'filled')
hold on
scatter(vars.O2_sat(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('\itc\rm_{sat}(O_2) [mmol m^{-3}]');

subplot(1,4,3);
scatter(vars.S(check),vars.P(check),'filled')
hold on
scatter(vars.S(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Salinity');

subplot(1,4,4);
scatter(vars.T(check),vars.P(check),'filled')
hold on
scatter(vars.T(check & check_low),vars.P(check & check_low),'filled')
ylim([0 46])
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
xlabel('Temperature [^\circC]');
set(gcf,'Color','w')

% create long-lat figure

figure;
scatter(vars.lon(check),vars.lat(check),'filled')
hold on
scatter(vars.lon(check & check_low),vars.lat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Latitude [^\circ N]');
xlabel('Longitude [^\circ E]');
set(gcf,'Color','w')

% create O2 in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.O2(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc\rm(O_2) [mmol m^{-3}]');
datetick
set(gcf,'Color','w')

% create MLD in time plot

figure('units','normalized','position',[0 0.1 1 .5]);
scatter(vars.t(check),vars.MLD(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.MLD(check & check_low),'filled')
set(gca,'YDir','Reverse','Box','On','FontSize',16);
ylabel('Depth [m]');
datetick
set(gcf,'Color','w')


% create O2sat, T and S in time plot

figure('units','normalized','position',[0 0.1 1 .8]);

subplot(3,1,1);
scatter(vars.t(check),vars.O2_sat(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.O2_sat(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('\itc_{sat}\rm(O_2) [mmol m^{-3}]');
datetick

subplot(3,1,2);
scatter(vars.t(check),vars.T(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.T(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Temp. [^\circC]');
datetick

subplot(3,1,3);
scatter(vars.t(check),vars.S(check),'filled')
hold on
scatter(vars.t(check & check_low),vars.S(check & check_low),'filled')
set(gca,'Box','On','FontSize',16);
ylabel('Salinity');
datetick

set(gcf,'Color','w')