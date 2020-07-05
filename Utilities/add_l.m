function add_l(number,HorV)
% number = constant number for line
% HorV = horizontal (1) or vertical (2)?

if HorV == 1
    properties = gca
    hold on
    x = properties.XLim(1):0.1:properties.XLim(end);
    plot(x,ones(size(x))*number,'k:','LineWidth',2)
else
    properties = gca
    hold on
    y = properties.YLim(1):0.1:properties.YLim(end);
    plot(ones(size(y))*number,y,'k:','LineWidth',2)    
end
end