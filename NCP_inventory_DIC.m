%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_inventory.m

% Script to get oxygen inventory term for NCP mass balance

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('DIC inventory change | Calculating ');

%% Calculation
% get inventory differentials
DIC_inv.inv_integral = [means_struct.DIC_h].*options.h;
DIC_inv.differentials = diff([means_struct.DIC_h].*options.h); % mean concentration x depth of layer
% interpolate from between days to on the day (i.e. 0.5 -> 1, 1.5 -> 2 ..)
DIC_inv.diffrange = (options.dayrange(1)+0.5*options.interval : ...
    options.interval : ...
    options.dayrange(end)-0.5*options.interval)';
DIC_inv.wantedrange = (options.dayrange(1) + 0.5 * options.interval : ...
    options.interval : ...
    options.dayrange(end) + options.interval - 0.5 * options.interval)'- 0.5 * options.interval;
DIC_inv.inv = interp1(DIC_inv.diffrange, DIC_inv.differentials, DIC_inv.wantedrange,'linear','extrap');
% get standard deviation for estimate of errors, interpolate like means
DIC_inv.inv_std = diff([means_struct.DIC_std_h].*options.h); % mean concentration x depth of layer
DIC_inv.inv_std = interp1(DIC_inv.diffrange,DIC_inv.inv_std,DIC_inv.wantedrange,'linear','extrap');

disp('DIC inventory change | Finished. ');