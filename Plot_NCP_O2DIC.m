
figure('units','normalized','position',[.1 0 .6 .9]); 

%% O2 first
% Create axes
a_pos = [0.13 0.517601043024772 0.775 0.407398956975228];
axes('Parent',gcf,...
    'Position',[0.13 0.517601043024772 0.775 0.407398956975228]);
hold on;

ADV = [O2_adv.adv];
y = [O2_inv.inv(2:end-1), ADV(2:end-1)', O2_ase.ASE(2:end-1)', ([O2_ent.ent])'.*-1];

yneg = y;
yneg(yneg>0) = 0;
ypos = y;
ypos(ypos<0) = 0;
hold on
b1 = bar(yneg,'stack')
b2 = bar(ypos,'stack')


[b4 b4a] = boundedline(1:25,NCP_est_kz,errors.error_NCP_kz)
set(b4a,'FaceAlpha',0.1,'FaceColor',[.2 .6 .2],'LineStyle',':','EdgeColor',[.2 .6 .2],'LineWidth',2)
set(b4, 'LineWidth',2,'Color',[.2 .6 .2])
box on;
hold on;
[b5 b5a] = boundedline(1:25,NCP_est_kz_no_adv,errors.error_NCP_kz_no_adv)
set(b5a,'FaceAlpha',0.1,'FaceColor',[1 .2 .4],'LineStyle',':','EdgeColor',[1 .2 .4],'LineWidth',2)
set(b5, 'LineWidth',2,'Color',[1 .2 .4])

p5 = plot(NCP_est_kz,'LineWidth',6,'Color','k');
p5 = plot(NCP_est_kz,'LineWidth',3,'Color',[.0 .6 .0]);

p6 = plot(NCP_est_kz_no_adv,'LineWidth',6,'Color','k');
p6 = plot(NCP_est_kz_no_adv,'LineWidth',3,'Color',[1 .2 .4]);

b1(1).FaceColor = [0 0.2 0.6];
b2(1).FaceColor = [0 0.2 0.6];

b1(2).FaceColor = [0.2 0.6 0.6];
b2(2).FaceColor = [0.2 0.6 0.6];

b1(3).FaceColor = [0.8 0.8 0.8];
b2(3).FaceColor = [0.8 0.8 0.8];

b1(4).FaceColor = [0.8 0.8 0.4];
b2(4).FaceColor = [0.8 0.8 0.4];

% b1(5).FaceColor = [0.8 0.4 0.6];
% b2(5).FaceColor = [0.8 0.4 0.6];

set(gca,'FontSize',18,'XTick',[1:length(options.dayrange(2:end-1))],...
    'LineWidth',3,'XTickLabels',[num2cell(options.dayrange(2:end-1))])
y1 = ylabel('Flux [mmol O_2 m^{-2} d^{-1}]')
box on
xlim([0 length(options.dayrange(2:end-1))+1])
ylim([-200 350])

set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}],...
    'YTick',[-300 -200 -100 0 100 200 300],'XGrid','On')

add_zero

l = legend('^{\partial I}/_{\partial t}','\itF_{\rmADV}','\itF_{\rmASE}','\itF_{\rmENT}','Location','NorthWest','Orientation','horizontal');
set(l,'FontSize',16,'Box','Off');


a2 = axes;
set(a2,'Visible','Off','XTickLabels','','YTickLabels','')

l2 = legend(a2,[p5],'\itN','Location','NorthWest','Orientation','Horizontal')
pos = get(l2,'Position');
pos(1) = 0.2;
pos(2) = 0.82;
set(l2,'FontSize',16,'Box','Off','Position',pos);

a3 = axes;
set(a3,'Visible','Off','XTickLabels','','YTickLabels','')

l3 = legend(a3,[p6],'\itN_{no ADV}','Location','NorthWest','Orientation','Horizontal')
pos = get(l3,'Position');
pos(1) = 0.30;
pos(2) = 0.82;
set(l3,'FontSize',16,'Box','Off','Position',pos);

set(gcf,'Color','W')

%% DIC
% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.775 0.398474576271186]);
a_pos = get(gca,'Position');
hold on;
ADV = [DIC_adv.adv];
y = [DIC_inv.inv(2:16), ADV(2:16)', DIC_ase.FDIC(2:16)', ([DIC_ent([1:16]).ent])'.*-1];

yneg = y;
yneg(yneg>0) = 0;
ypos = y;
ypos(ypos<0) = 0;
hold on
b1 = bar(yneg,'stack')
b2 = bar(ypos,'stack')


[b4 b4a] = boundedline(1:15,NCP_est_kz_DIC(1:15),errors.error_NCP_DIC_kz(1:15))
set(b4a,'FaceAlpha',0.1,'FaceColor',[.2 .6 .2],'LineStyle',':','EdgeColor',[.2 .6 .2],'LineWidth',2)
set(b4, 'LineWidth',2,'Color',[.2 .6 .2])
box on;
hold on;
[b5 b5a] = boundedline(1:15,NCP_est_kz_no_adv_DIC(1:15),errors.error_NCP_DIC_kz_no_adv(1:15))
set(b5a,'FaceAlpha',0.1,'FaceColor',[1 .2 .4],'LineStyle',':','EdgeColor',[1 .2 .4],'LineWidth',2)
set(b5, 'LineWidth',2,'Color',[1 .2 .4])

p5 = plot(1:15,NCP_est_kz_DIC(1:15),'LineWidth',6,'Color','k');
p5 = plot(1:15,NCP_est_kz_DIC(1:15),'LineWidth',3,'Color',[.0 .6 .0]);

p6 = plot(1:15,NCP_est_kz_no_adv_DIC(1:15),'LineWidth',6,'Color','k');
p6 = plot(1:15,NCP_est_kz_no_adv_DIC(1:15),'LineWidth',3,'Color',[1 .2 .4]);

b1(1).FaceColor = [0 0.2 0.6];
b2(1).FaceColor = [0 0.2 0.6];

b1(2).FaceColor = [0.2 0.6 0.6];
b2(2).FaceColor = [0.2 0.6 0.6];

b1(3).FaceColor = [0.8 0.8 0.8];
b2(3).FaceColor = [0.8 0.8 0.8];

b1(4).FaceColor = [0.8 0.8 0.4];
b2(4).FaceColor = [0.8 0.8 0.4];

% b1(5).FaceColor = [0.8 0.4 0.6];
% b2(5).FaceColor = [0.8 0.4 0.6];


set(gca,'FontSize',18,'XTick',[1:length(options.dayrange(2:end-1))],...
    'LineWidth',3,'XTickLabels',[num2cell(options.dayrange(2:end-1))])
y1 = ylabel('Flux [mmol C m^{-2} d^{-1}]')
box on
xlim([0 length(options.dayrange(2:end-1))+1])
ylim([-500 400])

set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}],...
    'YTick',[-500 -400 -300 -200 -100 0 100 200 300],'XGrid','On')

add_zero

%% save figure

print(gcf, '-dpng','-r400', [options.plot_dir,'NCP_O2_DIC'])

clear a2 ans b1 b2 l l2 p5 pos y yneg ypos b4 b4a ADV




