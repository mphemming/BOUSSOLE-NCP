figure('units','normalized','position',[.1 .3 .9 .7]); hold on;
ADV = [O2_adv.adv];
y = [O2_inv.inv(2:end-1), ADV(2:end-1)', O2_ase.ASE(2:end-1)', ([O2_ent.ent])'.*-1,[kz(2:end-1).kz]'];

yneg = y;
yneg(yneg>0) = 0;
ypos = y;
ypos(ypos<0) = 0;
hold on
b1 = bar(yneg,'stack')
b2 = bar(ypos,'stack')

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

b1(5).FaceColor = [0.6 0.4 0.4];
b2(5).FaceColor = [0.6 0.4 0.4];

set(gca,'FontSize',22,'XTick',[1:length(options.dayrange(2:end-1))],...
    'LineWidth',3,'XTickLabels',[num2cell(options.dayrange(2:end-1))])
ylabel('Flux [mmol O_2 m^{-2} d^{-1}]')
xlabel('March - April [2016]')
box on
xlim([0 length(options.dayrange(2:end-1))+1])
ylim([-300 600])

set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}])

l = legend('^{\partial I}/_{\partial t}','\itF_{\rmADV}',...
    '\itF_{\rmASE}','\itF_{\rmENT}','\itF_{\rmkz}',...
    'Location','NorthWest','Orientation','horizontal');
set(l,'FontSize',24,'Box','Off');

a2 = axes;
set(a2,'Visible','Off','XTickLabels','','YTickLabels','')

l2 = legend(a2,[p5],'\itN','Location','NorthWest','Orientation','Horizontal')
pos = get(l2,'Position');
pos(1) = 0.25;
pos(2) = 0.765;
set(l2,'FontSize',24,'Box','Off','Position',pos);


a3 = axes;
set(a3,'Visible','Off','XTickLabels','','YTickLabels','')

l3 = legend(a3,[p6],'\itN_{no ADV}','Location','NorthWest','Orientation','Horizontal')
pos = get(l3,'Position');
pos(1) = 0.32;
pos(2) = 0.76;
set(l3,'FontSize',24,'Box','Off','Position',pos);


set(gcf,'Color','W')

axes('Parent',gcf,...
    'Position',[0.543086419753086 0.58 0.324845679012346 0.3]);

[b4 b4a] = boundedline(1:25,NCP_est_kz,errors.error_NCP_kz)
set(b4a,'FaceAlpha',0.1,'FaceColor',[.2 .6 .2],'LineStyle',':','EdgeColor',[.2 .6 .2],'LineWidth',2)
set(b4, 'LineWidth',2,'Color',[.2 .6 .2])
box on;
hold on
[b5 b5a] = boundedline(1:25,NCP_est_kz_no_adv,errors.error_NCP_kz_no_adv)
set(b5a,'FaceAlpha',0.1,'FaceColor',[1 .2 .4],'LineStyle',':','EdgeColor',[1 .2 .4],'LineWidth',2)
set(b5, 'LineWidth',2,'Color',[1 .2 .4])

set(gca,'FontSize',16,'LineWidth',2)
set(gca,'XTick',[1:3:25],'XTickLabels',[{'10'} {'13'} {'16'} {'19'} {'22'} {'25'} {'28'} {'31'} {'3'}])
xlim([0 26]); ylim([-600 600]);
ylabel('Flux [mmol O_2 m^{-2} d^{-1}]')
xlabel('March - April [2016]')
hold on;
plot(1:25,zeros(1,25),'k','LineStyle','--')

 clear a2 a3 ans b1 b2 l l2 l3 p5 pos y yneg ypos b4 b4a b5 b5a ADV
 
 print(gcf, '-dpng','-r400', [options.plot_dir,'NCP_with_kz'])
 
 