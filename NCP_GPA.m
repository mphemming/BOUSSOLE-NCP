%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_GPA.m

% Script to calculate binned geopotential anomaly profiles and other variables..
% Must be run within NCP.m

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 21/06/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get binned profiles per dive
disp('Geopotential anomalies | binning data profiles');
vertical_grid = 1:10:1000;
for dive = 1:147
    check_dive = vars.dive == dive;
    % temp
    [vars.bin(dive).T] = bin_variable_profile(vars.T(check_dive),...
        vars.depth(check_dive),vertical_grid,0);
    % psal
    [vars.bin(dive).S] = bin_variable_profile(vars.S(check_dive),...
        vars.depth(check_dive),vertical_grid,0);
    % time
    [vars.bin(dive).t] = bin_variable_profile(vars.t(check_dive),...
        vars.depth(check_dive),vertical_grid,0);  
    % press
    [vars.bin(dive).P] = bin_variable_profile(vars.P(check_dive),...
        vars.depth(check_dive),vertical_grid,0);
    % lon
    [vars.bin(dive).lon] = bin_variable_profile(vars.lon(check_dive),...
        vars.depth(check_dive),vertical_grid,0);    
    % lat
    [vars.bin(dive).lat] = bin_variable_profile(vars.lat(check_dive),...
        vars.depth(check_dive),vertical_grid,0);    
    % absS
    [vars.bin(dive).absS] = bin_variable_profile(vars.absS(check_dive),...
        vars.depth(check_dive),vertical_grid,0);      
    % consT
    [vars.bin(dive).consT] = bin_variable_profile(vars.consT(check_dive),...
        vars.depth(check_dive),vertical_grid,0);      
    % define bin depth
    vars.bin(dive).depth = vars.bin(dive).consT.vertical;
end
disp('Geopotential anomalies | data profiles binned');

%%  Calculate geopotential anomalies
disp('Geopotential anomalies | calculating GP anomalies');
for dive = 1:147
        vars.bin(dive).GPA = gsw_geo_strf_dyn_height(...
            vars.bin(dive).absS.median_var,...
            vars.bin(dive).consT.median_var,...
            vars.bin(dive).P.median_var,0); % use abs S and consT, relative to surface
end
disp('Geopotential anomalies | GP anomalies calculated');

clear ans check_dive dive vertical_grid
