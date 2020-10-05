function [handle] = add_bathymetry_to_plot(x,y)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add_bathymetry_to_plot.m
%
% This script:
%
% 1. adds patches to plot when given x and y coordinates. Most useful for
% adding bathymetry to plot
%
% Created 30/01/2019 by MPH, NSW-IMOS Sydney
% Last updated 05/08/2019
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
% 
% x              |   x values (e.g. longitude)
% y              |   y values (e.g. bathymetry)
%
% Output:
% ---------------------------------------------------------------------------------------------------------------------------------
%
% handle         |  handle associated with patch
% 
% patch has been added to current axes
%
% ---------------------------------------------------------------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use absolute y values (for case that bathymetry is negative)
if nansum(y) < 0
    y = y*-1;
end
% link y and x values for patch
x(end+1) = x(1);
y(end+1) = y(end);
hold on
%produce patch
handle = patch(x,y,'k');    

end