



for n = 50:80

     figure('units','normalized','position',[0 0.1 .35 .8]);
    check = vars.profile_n == n;
    plot(vars.T(check),vars.depth(check),'LineWidth',2);
    add_l(profs.MLD(n),1);
    title(['Prof ',num2str(n)]);
    ylim([0 200]);
    set(gca,'YDir','Reverse','FontSize',20,'LineWidth',2);
    [xi(n),yi(n)] = getpts(gcf)
    
end