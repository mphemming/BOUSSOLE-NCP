figure('units','normalized','position',[.1 .3 .9 .7]); hold on;
y = [DIC_inv.inv(2:17), ADV_DIC(2:17)', DIC_ase.FDIC(2:17)', ([DIC_ent(2:17).ent])'.*-1];

yneg = y;
yneg(yneg>0) = 0;
ypos = y;
ypos(ypos<0) = 0;
hold on
b1 = bar(yneg,'stack')
b2 = bar(ypos,'stack')

p5 = plot(NCP_est_DIC(1:16),'LineWidth',6,'Color','k');
p5 = plot(NCP_est_DIC(1:16),'LineWidth',3,'Color',[.0 .6 .0]);

p6 = plot(NCP_est_no_adv_DIC(1:16),'LineWidth',6,'Color','k');
p6 = plot(NCP_est_no_adv_DIC(1:16),'LineWidth',3,'Color',[1 .2 .4]);

b1(1).FaceColor = [0 0.2 0.6];
b2(1).FaceColor = [0 0.2 0.6];

b1(2).FaceColor = [0.2 0.6 0.6];
b2(2).FaceColor = [0.2 0.6 0.6];

b1(3).FaceColor = [0.8 0.8 0.8];
b2(3).FaceColor = [0.8 0.8 0.8];

b1(4).FaceColor = [0.8 0.8 0.4];
b2(4).FaceColor = [0.8 0.8 0.4];

set(gca,'FontSize',22,'XTick',[1:length(options.dayrange(2:1))],...
    'LineWidth',3,'XTickLabels',[num2cell(options.dayrange(2:end-1))])
y1 = ylabel('Flux [mmol C m^{-2} d^{-1}]')
y1_pos = get(y1,'Position');
y1_pos(2) = -200; set(y1,'Position',y1_pos);
xlabel('March [2016]')
box on
xlim([0 17])
ylim([-800 200])

set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}],...
    'YTick',[-600 -500 -400 -300 -200 -100 0 100 200 300])

l = legend('^{\partial I}/_{\partial t}','\itF_{\rmADV}','\itF_{\rmASE}','\itF_{\rmENT}','Location','NorthWest','Orientation','horizontal');
leg_pos = get(l,'Position');
leg_pos(1) = 0.45;
set(l,'FontSize',24,'Box','Off','Position',leg_pos);


a2 = axes;
set(a2,'Visible','Off','XTickLabels','','YTickLabels','')

l2 = legend(a2,[p5],'\itN','Location','NorthWest','Orientation','Horizontal')
pos = get(l2,'Position');
pos(1) = 0.74;
pos(2) = 0.86;
set(l2,'FontSize',24,'Box','Off','Position',pos);


a3 = axes;
set(a3,'Visible','Off','XTickLabels','','YTickLabels','')

l3 = legend(a3,[p6],'\itN_{no ADV}','Location','NorthWest','Orientation','Horizontal')
pos = get(l3,'Position');
pos(1) = 0.8;
pos(2) = 0.855;
set(l3,'FontSize',24,'Box','Off','Position',pos);

set(gcf,'Color','W')

axes('Parent',gcf,...
    'Position',[0.2 0.21 0.324845679012346 0.3]);

[b4 b4a] = boundedline(1:16,NCP_est_DIC(1:16),errors.error_NCP_DIC(1:16))
set(b4a,'FaceAlpha',0.1,'FaceColor',[.2 .6 .2],'LineStyle',':','EdgeColor',[.2 .6 .2],'LineWidth',2)
set(b4, 'LineWidth',2,'Color',[.2 .6 .2])
box on;
hold on;
[b5 b5a] = boundedline(1:16,NCP_est_no_adv_DIC(1:16),errors.error_NCP_no_adv_DIC(1:16))
set(b5a,'FaceAlpha',0.1,'FaceColor',[1 .2 .4],'LineStyle',':','EdgeColor',[1 .2 .4],'LineWidth',2)
set(b5, 'LineWidth',2,'Color',[1 .2 .4])

ylim([-1650 1000]);
set(gca,'FontSize',16,'LineWidth',2)
set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}])
xlim([0.5 16.5])
ylabel('Flux [mmol C m^{-2} d^{-1}]')
xlabel('March [2016]')
hold on;
plot(1:16,zeros(1,16),'k','LineStyle','--')

clear a2 ans b1 b2 l l2 p5 pos y yneg ypos b4 b4a ADV

 
 print(gcf, '-dpng','-r400', [options.plot_dir,'NCP_DIC'])
