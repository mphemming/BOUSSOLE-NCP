function [binned_profile] = bin_variable_profile(variable,vertical,vertical_grid,interp_option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% bin_variable_profile.m
%
% This script:
%
% 1. Bin a variable vertically using pressure, depth, or density
%
% Created 16/10/2018 by MPH, NSW-IMOS Sydney
% Last updated 16/10/2018
% Email: m.hemming@unsw.edu.au
%
%%%%%%%%%%%%%%%%
% Making a modification?
%%%%%%%%%%%%%%%%
% Add any update information here
% Name        Date        Information
%-----------------------------------------------------
%
% Input:
% ---------------------------------------------------------------------------------------------------------------------------------
% variable                               | any variable (T,S,O2)
% vertical                                | can be pressure (dbar), depth (m), or density (kgm3)
% vertical_grid                        | set vertical grid spacing (e.g. [0:1:100], for surf to 100m every metre)
% interp_option                      | set = 1 for simple linear interpolation, set 0 for just binning
%
% Output:
% ---------------------------------------------------------------------------------------------------------------------------------
%
% binned_profile   | structure containing vertically-binned mean, median,
%                             and stdev, and interpolated profiles if required
%
% ---------------------------------------------------------------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create binned profile

% what is the regular interval in the grid
grid_spacing = nanmean(diff(vertical_grid));
% produce binned profiles
for chosen_bin = 1:numel(vertical_grid)
    % find variable data points within the bin
    check_bin = vertical >= vertical_grid(chosen_bin)-grid_spacing/2 & vertical < vertical_grid(chosen_bin)+grid_spacing/2; 
    % save vertical parameter
    binned_profile.vertical(chosen_bin) = vertical_grid(chosen_bin)+(grid_spacing/2);
    % save mean, median, and stdev
    binned_profile.mean_var(chosen_bin) = nanmean(variable(check_bin));
    binned_profile.median_var(chosen_bin) = nanmedian(variable(check_bin));    
    binned_profile.stdev_var(chosen_bin) = nanstd(variable(check_bin)); 
end

%% interpolate if user requires

if interp_option == 1
    % find non NaNs in profiles
    not_nans = ~isnan(binned_profile.vertical) & ~isnan(binned_profile.mean_var) & ...
        ~isnan(binned_profile.median_var) & ~isnan(binned_profile.stdev_var);
    % interpolate between binned points using finite values only
    binned_profile.mean_var_int = interp1(binned_profile.vertical(not_nans),binned_profile.mean_var(not_nans),...
                                                        vertical_grid,'Linear');
    binned_profile.median_var_int = interp1(binned_profile.vertical(not_nans),binned_profile.median_var(not_nans),...
                                                        vertical_grid,'Linear');
    binned_profile.stdev_var_int = interp1(binned_profile.vertical(not_nans),binned_profile.stdev_var(not_nans),...
                                                        vertical_grid,'Linear');                                                    
end

%disp('Voila!');

end