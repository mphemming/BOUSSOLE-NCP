function add_zero
properties = gca
hold on
x = properties.XLim(1):0.1:properties.XLim(end);
plot(x,zeros(size(x)),'k:','LineWidth',2)

end