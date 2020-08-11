%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_inventory_LOWESS.m

% Script to get oxygen inventory term for NCP mass balance
% Using LOWESS method

% created by MPH in Norwich, 03/08/2020
% modified by MPH in Sydney, 1/08/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Oxygen inventory change | Calculating ');

%% testing local linear regression curve
% figure;
for n = 1:numel(options.dayrange)
    check = vars.t >= datenum(2016,03,options.dayrange(n) - 2) & vars.t <= datenum(2016,03,options.dayrange(n) + 2) & ...
        vars.depth <= 46 & isfinite(vars.O2);
    [lin_test(n).fit,lin_test(n).gof,lin_test(n).out] = fit(vars.t(check)',vars.O2(check)','poly1');
    lin_test(n).inv = lin_test(n).fit(datenum(2016,03,options.dayrange(n)));
    lin_test(n).inv_error = lin_test(n).gof.rmse;
    lin_test(n).inv_r2 = lin_test(n).gof.rsquare;    
%     plot(lin_test(n).fit,vars.t(check)',vars.O2(check)');
%     hold on
    % get standard error of fits
    mdl = fitlm(vars.t(check)',vars.O2(check)');
    today = datenum(2016,03,options.dayrange(n));
    error(n) = mdl.Coefficients.Estimate(2)-mdl.Coefficients.SE(2)*(datenum(2016,03,10))+mdl.Coefficients.SE(1);
    
end







%% Calculation

% grid oxygen in time (every minute)
time_grid = nanmin(vars.t):datenum(0,0,0,0,15,0):nanmax(vars.t);
check_nans = isfinite(vars.O2) & vars.depth <= 46;
O2_gridded = interp1(vars.t(check_nans),vars.O2(check_nans),time_grid,'linear');
% get LOES fit
[smooth_struct] = LOESS_window(time_grid,O2_gridded,4);
% interp to get inventory change on days
check_nans = isfinite(O2_gridded);
smooth_interp = interp1(smooth_struct.time,smooth_struct.smoothed,datenum(2016,03,09):1:datenum(2016,04,03),'Linear');

% NOTES:
% Use python, play around with loess fits, or bootstrap / perform multiple,
% could also try local linear regressions, DEPTH RESOLVED LIKE ADVECTION!
% loess fits and get percentiles as error
% https://stackoverflow.com/questions/31104565/confidence-interval-for-lowess-in-python

% get mean O2 per dive
for n_dive = 1:147
    check = isfinite(vars.O2) & vars.depth <=46 & vars.dive == n_dive;
    O2_gridded(n_dive) = nanmean(vars.O2(check));
    time_gridded(n_dive) = nanmedian(vars.t(check));
end
% grid O2 daily
check_nans = isfinite(O2_gridded);
O2_interp = interp1(time_gridded(check_nans),O2_gridded(check_nans),datenum(2016,03,09):1:datenum(2016,04,03),'Linear');

% get LOES fit
[smooth_struct] = LOESS_window(datenum(2016,03,09):1:datenum(2016,04,03),O2_interp,2);
% interp to get inventory change on days
check_nans = isfinite(O2_gridded);
smooth_interp = interp1(smooth_struct.time,smooth_struct.smoothed,datenum(2016,03,09):1:datenum(2016,04,03),'Linear');


% get LOWESS fit
check_nans = isfinite(O2_gridded);
f = fit([time_grid(check_nans)',O2_gridded(check_nans)'],ones(size(O2_gridded(check_nans)')),'lowess')

% get inventory differentials






O2_inv.inv_integral = [means_struct.O2_h].*options.h;

O2_inv.differentials = diff([means_struct.O2_h].*options.h); % mean concentration x depth of layer
% interpolate from between days to on the day (i.e. 0.5 -> 1, 1.5 -> 2 ..)
O2_inv.diffrange = (options.dayrange(1)+0.5*options.interval : ...
    options.interval : ...
    options.dayrange(end)-0.5*options.interval)';
O2_inv.wantedrange = (options.dayrange(1) + 0.5 * options.interval : ...
    options.interval : ...
    options.dayrange(end) + options.interval - 0.5 * options.interval)'- 0.5 * options.interval;
O2_inv.inv = interp1(O2_inv.diffrange, O2_inv.differentials, O2_inv.wantedrange,'linear','extrap');
% get standard deviation for estimate of errors, interpolate like means
O2_inv.inv_std = diff([means_struct.O2_std_h].*options.h); % mean concentration x depth of layer
O2_inv.inv_std = interp1(O2_inv.diffrange,O2_inv.inv_std,O2_inv.wantedrange,'linear','extrap');

disp('Oxygen inventory change | Finished. ');