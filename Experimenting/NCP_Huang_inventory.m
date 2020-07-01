%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_dens_inventory.m

% Script to get oxygen inventory term for NCP mass balance

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Oxygen inventory change | Calculating ');

%% Calculation
% get inventory differentials
O2_inv.inv_integral_h1 = [means_struct.O2_h1].*[means_struct.MLD];
O2_inv.differentials_h1 = diff([means_struct.O2_h1].*[means_struct.MLD]); % mean concentration x depth of layer
O2_inv.inv_integral_h2 = [means_struct.O2_h2].*options.h;
O2_inv.differentials_h2 = diff([means_struct.O2_h2].*options.h); % mean concentration x depth of layer
% interpolate from between days to on the day (i.e. 0.5 -> 1, 1.5 -> 2 ..)
O2_inv.diffrange = (options.dayrange(1)+0.5*options.interval : ...
    options.interval : ...
    options.dayrange(end)-0.5*options.interval)';
O2_inv.wantedrange = (options.dayrange(1) + 0.5 * options.interval : ...
    options.interval : ...
    options.dayrange(end) + options.interval - 0.5 * options.interval)'- 0.5 * options.interval;
O2_inv.inv_h1 = interp1(O2_inv.diffrange, O2_inv.differentials_h1, O2_inv.wantedrange,'linear','extrap');
O2_inv.inv_h2 = interp1(O2_inv.diffrange, O2_inv.differentials_h2, O2_inv.wantedrange,'linear','extrap');
% get standard deviation for estimate of errors, interpolate like means
O2_inv.inv_std_h1 = diff([means_struct.O2_std_h1].*[means_struct.MLD]); % mean concentration x depth of layer
O2_inv.inv_std_h1 = interp1(O2_inv.diffrange,O2_inv.inv_std_h1,O2_inv.wantedrange,'linear','extrap');
O2_inv.inv_std_h2 = diff([means_struct.O2_std_h2].*options.h); % mean concentration x depth of layer
O2_inv.inv_std_h2 = interp1(O2_inv.diffrange,O2_inv.inv_std_h2,O2_inv.wantedrange,'linear','extrap');

disp('Oxygen inventory change | Finished. ');