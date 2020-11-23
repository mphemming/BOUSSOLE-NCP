%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_GPA.m

% Script to calculate binned geopotential anomaly profiles and other variables..
% Must be run within NCP.m

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 21/06/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get binned profiles per profile
disp('Geopotential anomalies | binning data profiles');
for n_prof = 1:294
    check_prof = vars.profile_n == n_prof;
    % temp
    [vars.bin(n_prof).T] = bin_variable_profile_NCP(vars.T(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);
    % psal
    [vars.bin(n_prof).S] = bin_variable_profile_NCP(vars.S(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);
    % time
    [vars.bin(n_prof).t] = bin_variable_profile_NCP(vars.t(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);  
    % press
    [vars.bin(n_prof).P] = bin_variable_profile_NCP(vars.P(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);
    % lon
    [vars.bin(n_prof).lon] = bin_variable_profile_NCP(vars.lon(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);    
    % lat
    [vars.bin(n_prof).lat] = bin_variable_profile_NCP(vars.lat(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);    
    % absS
    [vars.bin(n_prof).absS] = bin_variable_profile_NCP(vars.absS(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);      
    % consT
    [vars.bin(n_prof).consT] = bin_variable_profile_NCP(vars.consT(check_prof),...
        vars.depth(check_prof),options.vertical_grid,0);      
    % define bin depth
    vars.bin(n_prof).depth = vars.bin(n_prof).consT.vertical;
end
disp('Geopotential anomalies | data profiles binned');

%%  Calculate geopotential anomalies
disp('Geopotential anomalies | calculating GP anomalies');
for n_prof = 1:294
        vars.bin(n_prof).GPA = gsw_geo_strf_dyn_height(...
            vars.bin(n_prof).absS.median_var,...
            vars.bin(n_prof).consT.median_var,...
            vars.bin(n_prof).P.median_var,0); % use abs S and consT, relative to surface
end
disp('Geopotential anomalies | GP anomalies calculated');

clearvars -except options prcdata vars vars_profile_means Oscar AVISO
