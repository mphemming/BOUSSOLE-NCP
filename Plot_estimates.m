%--------------------------------------------------------------------
%--------------------------------------------------------------------
%------------
% Plot_estimates.m
%------------
%-----------------------------------------------------------------------------------------------------
% script by MPH in Norwich, 05/07/2020
% script modified by MPH in Sydney, 05/07/2020
%

clear all
close all
clc

%-----------------------------------------------------------------------------------------------------
%% load data

options.directory = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts'
options.data_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/data/'
options.plot_dir = '/Volumes/Hemming_backup/UNSW_work/UEA/NCP_scripts/Plots/'
addpath(genpath(options.directory));

NCP_46 = load([options.data_dir,'NCP_46m.mat']);
NCP_30 = load([options.data_dir,'NCP_30m.mat']);
NCP_22 = load([options.data_dir,'NCP_22m.mat']);
%-----------------------------------------------------------------------------------------------------
%% create figure

figure('units','normalized','position',[.1 .1 .7 .8]);
%-----------------------------------------------------------------------------------------------------
% NCP with ADV
subplot(2,2,1)
hold on;
% March 20-25
% buoy
scatter(1,NCP_46.N_table.O2_buoy.M2025.NCP_ADV_46,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.O2_buoy.M2025.NCP_ADV_46,250,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(1,NCP_46.N_table.DIC_buoy.M2025.NCP_ADV_46*-1,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.DIC_buoy.M2025.NCP_ADV_46*-1,250,'d','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% glider
scatter(1,NCP_46.N_table.O2_glider.M2025.NCP,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.O2_glider.M2025.NCP,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(1,NCP_46.N_table.DIC_glider.M2025.NCP*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.DIC_glider.M2025.NCP*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% March 29-April 01
% buoy
scatter(2,NCP_46.N_table.O2_buoy.M2901.NCP_ADV_46,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.O2_buoy.M2901.NCP_ADV_46,250,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% scatter(2,NCP_46.N_table.DIC_buoy.M2901.NCP_ADV*-1,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
% scatter(2,NCP_46.N_table.DIC_buoy.M2901.NCP_ADV*-1,250,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% glider
scatter(2,NCP_46.N_table.O2_glider.M2901.NCP,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.O2_glider.M2901.NCP,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(2,NCP_46.N_table.DIC_glider.M2901.NCP*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.DIC_glider.M2901.NCP*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% Whole period
% glider
scatter(3,NCP_46.N_table.O2_glider.M1025.NCP,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(3,NCP_46.N_table.O2_glider.M1025.NCP,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(3,NCP_46.N_table.DIC_glider.M1025.NCP*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(3,NCP_46.N_table.DIC_glider.M1025.NCP*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);

title('\itN \rm(ADV)');
xlim([0.5 3.5]); ylim([-50 300]);
add_zero
set(gca,'XTick',[1 2 3],'XTickLabels',[{'M20 - M25'} {'M29 - A01'} {'M10 - M25'}], ...
    'FontSize',18,'LineWidth',2,'Box','On');
%-----------------------------------------------------------------------------------------------------
% NCP without ADV
subplot(2,2,2)
hold on;
% March 20-25
% buoy
scatter(1,NCP_46.N_table.O2_buoy.M2025.NCP_46,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.O2_buoy.M2025.NCP_46,250,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(1,NCP_46.N_table.DIC_buoy.M2025.NCP_46*-1,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.DIC_buoy.M2025.NCP_46*-1,250,'d','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% glider
scatter(1,NCP_46.N_table.O2_glider.M2025.NCP_no_adv,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.O2_glider.M2025.NCP_no_adv,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(1,NCP_46.N_table.DIC_glider.M2025.NCP_no_adv*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.DIC_glider.M2025.NCP_no_adv*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% March 29-April 01
% buoy
scatter(2,NCP_46.N_table.O2_buoy.M2901.NCP_46,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.O2_buoy.M2901.NCP_46,250,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% scatter(2,NCP_46.N_table.DIC_buoy.M2901.NCP_46*-1,400,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
% scatter(2,NCP_46.N_table.DIC_buoy.M2901.NCP_46*-1,250,'d','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% glider
scatter(2,NCP_46.N_table.O2_glider.M2901.NCP_no_adv,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.O2_glider.M2901.NCP_no_adv,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(2,NCP_46.N_table.DIC_glider.M2901.NCP_no_adv*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(2,NCP_46.N_table.DIC_glider.M2901.NCP_no_adv*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% Whole period
% glider
scatter(3,NCP_46.N_table.O2_glider.M1025.NCP_no_adv,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(3,NCP_46.N_table.O2_glider.M1025.NCP_no_adv,250,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(3,NCP_46.N_table.DIC_glider.M1025.NCP_no_adv*-1,400,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(3,NCP_46.N_table.DIC_glider.M1025.NCP_no_adv*-1,250,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);

title('\itN \rm(No ADV)');
xlim([0.5 3.5]); ylim([-50 300]);
add_zero
set(gca,'XTick',[1 2 3],'XTickLabels',[{'M20 - M25'} {'M29 - A01'} {'M10 - M25'}], ...
    'FontSize',18,'LineWidth',2,'Box','On');

% create legend
s1 = scatter(-1,NCP_46.N_table.O2_glider.M1025.NCP_no_adv,600,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
s2 = scatter(-1,NCP_46.N_table.O2_buoy.M2901.NCP_46,600,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
s3 = scatter(-1,NCP_46.N_table.O2_buoy.M2901.NCP_46,600,'k','Sq');
s4 = scatter(-1,NCP_46.N_table.O2_buoy.M2901.NCP_46,600,'k','d');

leg = legend([s1 s2 s3 s4],'\itN_{\rmDIC}','\itN_{\rmO_2}','Glider','Buoy');
clear leg_pos
leg_pos = leg.Position;
leg_pos(1) = 0.45; leg_pos(2) = 0.44; 
set(leg,'FontSize',22,'Orientation','Horizontal','Position',leg_pos,'Box','Off');
%-----------------------------------------------------------------------------------------------------
% Create NCP patches
subplot(2,2,3:4)
hold on;

% NCP period 1
p1(:,1) = [datenum(2016,03,20) datenum(2016,03,25) datenum(2016,03,25) datenum(2016,03,20)];
p1(:,2) = [500 500 -400 -400];
% NCP period 2
p2(:,1) = [datenum(2016,03,29) datenum(2016,04,02) datenum(2016,04,02) datenum(2016,03,29)];
p2(:,2) = [500 500 -400 -400];

hold on;
pp1 = patch(p1(:,1),p1(:,2),[.6 .6 .6]);
set(pp1,'FaceAlpha',0.1,'LineStyle','None');
pp2 = patch(p2(:,1),p2(:,2),[.6 .6 .6]);
set(pp2,'FaceAlpha',0.1,'LineStyle','None');

add_l(datenum(2016,03,20),2)
add_l(datenum(2016,03,25),2)
add_l(datenum(2016,03,29),2)
add_l(datenum(2016,04,02),2)
add_zero
set(gca,'YLim',[-200 400],'XLim',[datenum(2016,03,09) datenum(2016,04,04)])
%-----------------------------------------------------------------------------------------------------
% NCP time comparison
% glider with time
plot(datenum(2016,03,NCP_22.options.dayrange(2:end-1)),NCP_46.NCP_est_kz,'LineWidth',6,'Color','k');
plot(datenum(2016,03,NCP_22.options.dayrange(2:end-1)),NCP_46.NCP_est_kz,'LineWidth',3,'Color',[.0 .6 .6]);
plot(datenum(2016,03,NCP_22.options.dayrange(2:17)),NCP_46.NCP_est_kz_DIC(2:17)*-1,'LineWidth',6,'Color','k');
plot(datenum(2016,03,NCP_22.options.dayrange(2:17)),NCP_46.NCP_est_kz_DIC(2:17)*-1,'LineWidth',3,'Color',[.8 .8 .4]);
% NCP estimates
% March 20-25
% buoy
scatter(datenum(2016,03,23),NCP_46.N_table.O2_buoy.M2025.NCP_ADV_46,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(1,NCP_46.N_table.O2_buoy.M2025.NCP_ADV_46,50,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% scatter(datenum(2016,03,23),NCP_46.N_table.DIC_buoy.M2025.NCP_ADV*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
% scatter(datenum(2016,03,23),NCP_46.N_table.DIC_buoy.M2025.NCP_ADV*-1,50,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% glider
scatter(datenum(2016,03,23),NCP_46.N_table.O2_glider.M2025.NCP,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.O2_glider.M2025.NCP,50,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_glider.M2025.NCP*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_glider.M2025.NCP*-1,50,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% buoy
scatter(datenum(2016,03,23),NCP_46.N_table.O2_buoy.M2025.NCP_46,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.O2_buoy.M2025.NCP_46,50,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_buoy.M2025.NCP_46*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_buoy.M2025.NCP_46*-1,50,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% glider
scatter(datenum(2016,03,23),NCP_46.N_table.O2_glider.M2025.NCP_no_adv,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.O2_glider.M2025.NCP_no_adv,50,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_glider.M2025.NCP_no_adv*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,23),NCP_46.N_table.DIC_glider.M2025.NCP_no_adv*-1,50,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% March 29-April 01
% buoy
scatter(datenum(2016,03,31),NCP_46.N_table.O2_buoy.M2901.NCP_ADV_46,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.O2_buoy.M2901.NCP_ADV_46,50,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% scatter(datenum(2016,03,31),NCP_46.N_table.DIC_buoy.M2901.NCP_ADV*-1,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
% scatter(datenum(2016,03,31),NCP_46.N_table.DIC_buoy.M2901.NCP_ADV*-1,50,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
% glider
scatter(datenum(2016,03,31),NCP_46.N_table.O2_glider.M2901.NCP,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.O2_glider.M2901.NCP,50,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_glider.M2901.NCP*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_glider.M2901.NCP*-1,50,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% buoy
scatter(datenum(2016,03,31),NCP_46.N_table.O2_buoy.M2901.NCP_46,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.O2_buoy.M2901.NCP_46,50,'d','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_buoy.M2901.NCP_46*-1,100,'d','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_buoy.M2901.NCP_46*-1,50,'d','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);
% glider
scatter(datenum(2016,03,31),NCP_46.N_table.O2_glider.M2901.NCP_no_adv,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.O2_glider.M2901.NCP_no_adv,50,'Sq','MarkerFaceColor',[0 .6 .6],'MarkerEdgeColor',[0 .6 .6]);
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_glider.M2901.NCP_no_adv*-1,100,'Sq','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(datenum(2016,03,31),NCP_46.N_table.DIC_glider.M2901.NCP_no_adv*-1,50,'Sq','MarkerFaceColor',[.8 .8 .4],'MarkerEdgeColor',[.8 .8 .4]);

set(gca,'YLim',[-200 450],'XLim',[datenum(2016,03,09) datenum(2016,04,04)],'Box','On','LineWidth',2,'FontSize',18);
datetick('x','KeepLimits');
add_zero

% save figure
% print(gcf,'-dpng','-r400',[options.plot_dir,'NCP_estimates.png']);

%% Plot PQ values

figure('units','normalized','position',[.1 .4 .7 1]);
plot(datenum(2016,03,NCP_22.options.dayrange(2:end-1)),NCP_46.NCP_est_kz./NCP_46.NCP_est_kz_DIC)

% NOTES: how to compare like for like (depths buoy vs glider)?
% wind comparisons
% display PQ








