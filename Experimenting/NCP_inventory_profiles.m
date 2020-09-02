%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_inventory.m

% Script to get oxygen inventory term for NCP mass balance

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Oxygen inventory change | Calculating ');

%% Calculation
% get inventory differentials
O2_inv.inv_integral = [means_struct.O2_h_profs].*options.h;
O2_inv.differentials = diff([means_struct.O2_h_profs].*options.h); % mean concentration x depth of layer
% interpolate from between days to on the day (i.e. 0.5 -> 1, 1.5 -> 2 ..)
O2_inv.diffrange = (options.dayrange(1)+0.5*options.interval : ...
    options.interval : ...
    options.dayrange(end)-0.5*options.interval)';
O2_inv.wantedrange = (options.dayrange(1) + 0.5 * options.interval : ...
    options.interval : ...
    options.dayrange(end) + options.interval - 0.5 * options.interval)'- 0.5 * options.interval;
O2_inv.inv = interp1(O2_inv.diffrange, O2_inv.differentials, O2_inv.wantedrange,'linear','extrap');
% get standard deviation for estimate of errors, interpolate like means
O2_inv.inv_std = diff([means_struct.O2_std_h_profs].*options.h); % mean concentration x depth of layer
O2_inv.inv_std = interp1(O2_inv.diffrange,O2_inv.inv_std,O2_inv.wantedrange,'linear','extrap');

disp('Oxygen inventory change | Finished. ');