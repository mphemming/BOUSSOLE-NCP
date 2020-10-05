%vars.DIC = ([prcdata.CO2SYS.DIC]);
DIC_Snorm_glider = (prcdata.timeseries.DIC./prcdata.timeseries.S) .* 38.3; 
DIC_Snorm_glider = DIC_Snorm_glider .* (prcdata.timeseries.sigma0 / 1000);
v.t = prcdata.timeseries.t;
time_grid = datenum(2016,03,01):1/24:datenum(2016,04,10);
for b = 1:numel(time_grid)
    bin(b).DIC_glider = nanmean(DIC_Snorm_glider(prcdata.timeseries.close_vector <= 5 & ...
        vars.t >= time_grid(b)-0.5 & vars.t < time_grid(b)+0.5 & prcdata.timeseries.depth > 9 ...
        & prcdata.timeseries.depth < 11));
    bin(b).DIC_buoy = nanmean(BOUSSOLE.DIC_CSYS_Snorm_mmolm3( ...
        BOUSSOLE.time_CSYS >= time_grid(b)-0.5 & BOUSSOLE.time_CSYS < time_grid(b)+0.5));
    bin(b).t = time_grid(b);
end
buoy = [bin.DIC_buoy]';
glider = [bin.DIC_glider]';
t = [bin.t];
checknans = isfinite(buoy) & isfinite(glider);

[fit1, gof1, out1] = fit(glider(checknans' & t <= datenum(2016,03,25)),buoy(checknans' & t <= datenum(2016,03,25)'),'Poly1');




