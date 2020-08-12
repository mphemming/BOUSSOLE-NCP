function plot_rand_error(dates,term,error,title_string)

    
%% create plot

figure('units','normalized','position',[0 0.1 .8 .6]);

for n = 1:numel(dates)
    stds(n) = nanstd(error(n,:));
    upper(n) = term(n) + stds(n);
    lower(n) = term(n) - stds(n);
    err = term(n) + (error(n,:) - nanmean(error(n,:)));
    
    hold on
    s = scatter(dates(n)*ones(size(err)),err,10,'k','filled')
end
p1 = plot(dates,term,'g','LineWidth',2);
p2 = plot(dates,upper,'r','LineStyle','--');
plot(dates,lower,'r','LineStyle','--');
add_zero

set(gca,'FontSize',16,'Box','On','XGrid','On')
datetick('x','dd/mm','KeepLimits');
ylabel('Error [mmol m^{-3}]');
title(title_string)

leg = legend([s p2 p1],'Error range','Standard Deviation','Term')
set(leg,'Location','NorthWest','Box','Off');

annotation(gcf,'textbox',[0.136 0.4 0.2 0.292727272727272], ...
    'String',{['mean = ',num2str(nanmean(stds)),' mmol m^{-3}']},'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',[0.136 0.35 0.2 0.292727272727272], ...
    'String',{['min = ',num2str(nanmin(stds)),' mmol m^{-3}']},'FitBoxToText','off','LineStyle','None');
annotation(gcf,'textbox',[0.136 0.3 0.2 0.292727272727272], ...
    'String',{['max = ',num2str(nanmax(stds)),' mmol m^{-3}']},'FitBoxToText','off','LineStyle','None');

end