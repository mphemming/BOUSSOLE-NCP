%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_plane_fits.m

% Script to calculate planes in oxygen and geopotential anomalies 
% needed to calculate the advection term

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 05/10/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if saved file already, don't need to do all of below
testing =1;
if exist([options.data_dir,'planes_',num2str(options.h),'.mat'],'file') | testing == 1 

%% get concatenated geopotential anomaly profile arrays for planes

disp(['O2 & DIC Plane-fits | organising data']);

GPA = [vars.bin.GPA];    
GPA_D = [vars.bin.depth];    
initial = [vars.bin.t]; GPA_t = [initial.median_var];
initial = [vars.bin.lon]; lon = [initial.median_var];
initial = [vars.bin.lat]; lat = [initial.median_var];

clear initial

%% extract profile data

fields = fieldnames(vars_profile_means);
for n = 1:numel(fields)
    eval(['profs.',cell2mat(fields(n)),' = [vars_profile_means.',cell2mat(fields(n)),'];']);
end


%% calculate cartesian coordinates for O2 Plane-fits
% Compute x (longitude) & y (latitude) distances relative to the origin:

disp(['O2 Plane-fits | Calculating cartesian coordinates']);

% chosen lat/lon outside of domain
lat0 = 43.2;
lon0 = 7.5;

% For O2
for i = 1:numel(profs.P)
    [profs.xa(i),~] = sw_dist([lat0; lat0],[lon0; profs.lon(i)],'km');
    [profs.ya(i),~] = sw_dist([lat0; profs.lat(i)],[lon0; lon0],'km');
end
for i = 1:numel(vars.P)
    [vars.xa(i),~] = sw_dist([lat0; lat0],[lon0; vars.lon(i)],'km');
    [vars.ya(i),~] = sw_dist([lat0; vars.lat(i)],[lon0; lon0],'km');
end

% For GPA
for i = 1:numel(lon)
    [xa_bin(i),~] = sw_dist([lat0; lat0],[lon0; lon(i)],'km');
    [ya_bin(i),~] = sw_dist([lat0; lat(i)],[lon0; lon0],'km');
end    

clear lat0 lon0 i

%% Loop | for each day +- n days moving window, 
% obtain O2 Plane-fits of oxygen and geopotential anomalies

disp(['O2 Plane-fits | Obtaining gradients for each time step window']);

for day = options.dayrange
    
    % get datenum for particular day
    if day <= 31
    planes_loop(day).date_num = datenum(2016,03,day,00,00,00);
    else
    planes_loop(day).date_num = datenum(2016,04,day-31,00,00,00);
    end    

    %% select data within moving time window
    
    planes_loop(day).day_selection_h = vars.t > (planes_loop(day).date_num-options.window) & ...
        vars.t < (planes_loop(day).date_num+options.window) & ...
        vars.P <= options.h;  % only top h metres for most variables 
   planes_loop(day).day_selection_h_profs = profs.t > (planes_loop(day).date_num-options.window) & ...
        profs.t < (planes_loop(day).date_num+options.window) & ...
        profs.P <= options.h;  % only top h metres for most variables    
    planes_loop(day).day_selection_h_adv = vars.t > (planes_loop(day).date_num-(options.window*options.window_adv)) & ...
        vars.t < (planes_loop(day).date_num+(options.window*options.window_adv)) & ...
        vars.P <= options.h;  % Using the longer time window for ADV calculation!        
    planes_loop(day).day_selection_GPA = ...
        GPA_t > (planes_loop(day).date_num-(options.window*options.window_adv)) ...
        & GPA_t < (planes_loop(day).date_num+(options.window*options.window_adv)); % for binned GPA profiles
    planes_loop(day).day_selection_DAC = ...
        vars.DACs.t > (planes_loop(day).date_num-(options.window*options.window_adv)) & ...
        vars.DACs.t < (planes_loop(day).date_num+(options.window*options.window_adv)); % for DACs
    % get number of days in window
    [~,m,d,~,~,~] = datevec(profs.t(planes_loop(day).day_selection_h_profs));
    planes_loop(day).days_in_window_prof = numel(unique(datenum(0,m,d)));

    %% Select data for O2 Plane-fits with NaNs removed
    % all measurements
    check_nan = isfinite(vars.O2) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).oxygen_planes.dive = vars.dive(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.Chl = vars.Chl(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.O2 = vars.O2(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.time = vars.t(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.xa = vars.xa(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.ya = vars.ya(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.P = vars.P(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.lon = vars.lon(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.lat = vars.lat(planes_loop(day).day_selection_h & check_nan);    
    % profiles
    check_nan = isfinite(profs.O2) & isfinite(profs.xa) & isfinite(profs.ya);
    planes_loop(day).oxygen_planes.dive_profs = profs.dive(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.Chl_profs = profs.Chl(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.O2_profs = profs.O2(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.time_profs = profs.t(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.xa_profs = profs.xa(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.ya_profs = profs.ya(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.P_profs = profs.P(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.lon_profs = profs.lon(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).oxygen_planes.lat_profs = profs.lat(planes_loop(day).day_selection_h_profs & check_nan);
    % advection
    check_nan = isfinite(vars.O2) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).oxygen_planes.O2_adv = vars.O2(planes_loop(day).day_selection_h_adv & check_nan);    
    planes_loop(day).oxygen_planes.time_adv = vars.t(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).oxygen_planes.P_adv = vars.P(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).oxygen_planes.lon_adv = vars.lon(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).oxygen_planes.lat_adv = vars.lat(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).oxygen_planes.xa_adv = vars.xa(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).oxygen_planes.ya_adv = vars.ya(planes_loop(day).day_selection_h_adv & check_nan);    

    %% Select data for DIC Plane-fits with NaNs removed
    % all measurements
    check_nan = isfinite(vars.DIC) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).DIC_planes.dive = vars.dive(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.Chl = vars.Chl(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.DIC = vars.DIC(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.time = vars.t(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.xa = vars.xa(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.ya = vars.ya(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.P = vars.P(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.lon = vars.lon(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.lat = vars.lat(planes_loop(day).day_selection_h & check_nan);    
    % profiles
    check_nan = isfinite(profs.DIC) & isfinite(profs.xa) & isfinite(profs.ya);
    planes_loop(day).DIC_planes.dive_profs = profs.dive(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.Chl_profs = profs.Chl(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.DIC_profs = profs.DIC(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.time_profs = profs.t(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.xa_profs = profs.xa(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.ya_profs = profs.ya(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.P_profs = profs.P(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.lon_profs = profs.lon(planes_loop(day).day_selection_h_profs & check_nan);
    planes_loop(day).DIC_planes.lat_profs = profs.lat(planes_loop(day).day_selection_h_profs & check_nan);
    % advection
    check_nan = isfinite(vars.DIC) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).DIC_planes.DIC_adv = vars.DIC(planes_loop(day).day_selection_h_adv & check_nan);    
    planes_loop(day).DIC_planes.time_adv = vars.t(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).DIC_planes.P_adv = vars.P(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).DIC_planes.lon_adv = vars.lon(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).DIC_planes.lat_adv = vars.lat(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).DIC_planes.xa_adv = vars.xa(planes_loop(day).day_selection_h_adv & check_nan);
    planes_loop(day).DIC_planes.ya_adv = vars.ya(planes_loop(day).day_selection_h_adv & check_nan);        
    
    %% Select data for geopotential anomaly Plane-fits with NaNs removed
    check_nan = isfinite(GPA) & isfinite(xa_bin) & isfinite(ya_bin);
    planes_loop(day).GPA_planes.GPA = GPA(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.depth = GPA_D(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.time = GPA_t(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.xa = xa_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.ya = ya_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.D_levels = unique(planes_loop(day).GPA_planes.depth);
   
%% Oxygen depth resolved gridded plane fits
    
% grid oxygen with similar bins to geopotential anomalies
xgrid = round(nanmin(planes_loop(day).oxygen_planes.xa)):0.5:round(nanmax(planes_loop(day).oxygen_planes.xa));
ygrid = round(nanmin(planes_loop(day).oxygen_planes.ya)):0.5:round(nanmax(planes_loop(day).oxygen_planes.ya));
[X,Y] = meshgrid(xgrid,ygrid);

for n_bin = 1:numel(options.vertical_grid)
    for n_x = 1:numel(xgrid)
        for n_y = 1:numel(ygrid)
            check_bin = planes_loop(day).oxygen_planes.P_adv >= options.vertical_grid(n_bin) - 6 & planes_loop(day).oxygen_planes.P_adv < options.vertical_grid(n_bin) + 4 ...
                & planes_loop(day).oxygen_planes.xa_adv >= xgrid(n_x) & planes_loop(day).oxygen_planes.xa_adv < xgrid(n_x)+0.49 & ...
                planes_loop(day).oxygen_planes.ya_adv >= ygrid(n_y) & planes_loop(day).oxygen_planes.ya_adv < ygrid(n_y)+0.49;
            planes_loop(day).oxygen_planes.gridded_O2_var(n_bin,n_x,n_y) = nanmedian(planes_loop(day).oxygen_planes.O2_adv(check_bin));
            planes_loop(day).oxygen_planes.gridded_O2_var_std(n_bin,n_x,n_y) = nanstd(planes_loop(day).oxygen_planes.O2_adv(check_bin));
            planes_loop(day).oxygen_planes.gridded_xa(n_bin,n_x,n_y) = xgrid(n_x);
            planes_loop(day).oxygen_planes.gridded_ya(n_bin,n_x,n_y) = ygrid(n_y);    
            planes_loop(day).oxygen_planes.gridded_z(n_bin) =  options.vertical_grid(n_bin);
        end
    end
end

% get fits
for n_bin = 1:numel(options.vertical_grid)
    xa = squeeze(planes_loop(day).oxygen_planes.gridded_xa(n_bin,:,:)); 
    ya = squeeze(planes_loop(day).oxygen_planes.gridded_ya(n_bin,:,:)); 
    O2 = squeeze(planes_loop(day).oxygen_planes.gridded_O2_var(n_bin,:,:));
    check_nans = isfinite(O2);
    if numel(O2(check_nans)) > 1
      % prepare for fit
      planes_loop(day).oxygen_planes.gridded_const = ones(size(O2(check_nans)));
      planes_loop(day).oxygen_planes.gridded_plane_matrix = [xa(check_nans), ya(check_nans), planes_loop(day).oxygen_planes.gridded_const];
      [planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).fit,planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).gof,...
          planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).out] = fit(planes_loop(day).oxygen_planes.gridded_plane_matrix(:,1:2),O2(check_nans),'poly11');    
      planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).O2_for_fit = O2(check_nans);
      planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).xa_for_fit = planes_loop(day).oxygen_planes.gridded_plane_matrix(:,1);
      planes_loop(day).oxygen_planes.gridded_O2fit(n_bin).ya_for_fit = planes_loop(day).oxygen_planes.gridded_plane_matrix(:,2);
      % robust fit using bisquare method (less weighting to outliers)
      [planes_loop(day).oxygen_planes.gridded_O2fit_robust(n_bin).vals,stats] = robustfit(planes_loop(day).oxygen_planes.gridded_plane_matrix(:,1:2),O2(check_nans),'bisquare');
      planes_loop(day).oxygen_planes.gridded_m_lon_error(n_bin) = stats.se(2);
      planes_loop(day).oxygen_planes.gridded_m_lat_error(n_bin) = stats.se(3);
      planes_loop(day).oxygen_planes.robust_stats = stats;
      planes_loop(day).oxygen_planes.gridded_m_lon(n_bin) = planes_loop(day).oxygen_planes.gridded_O2fit_robust(n_bin).vals(2);
      planes_loop(day).oxygen_planes.gridded_m_lat(n_bin)  = planes_loop(day).oxygen_planes.gridded_O2fit_robust(n_bin).vals(3);
      planes_loop(day).oxygen_planes.gridded_m_0(n_bin)  = planes_loop(day).oxygen_planes.gridded_O2fit_robust(n_bin).vals(1);           
    end
end


%% DIC depth resolved gridded plane fits
    
% grid oxygen with similar bins to geopotential anomalies
xgrid = round(nanmin(planes_loop(day).DIC_planes.xa)):0.5:round(nanmax(planes_loop(day).DIC_planes.xa));
ygrid = round(nanmin(planes_loop(day).DIC_planes.ya)):0.5:round(nanmax(planes_loop(day).DIC_planes.ya));
[X,Y] = meshgrid(xgrid,ygrid);

for n_bin = 1:numel(options.vertical_grid)
    for n_x = 1:numel(xgrid)
        for n_y = 1:numel(ygrid)
            check_bin = planes_loop(day).DIC_planes.P_adv >= options.vertical_grid(n_bin) - 6 & planes_loop(day).DIC_planes.P_adv < options.vertical_grid(n_bin) + 4 ...
                & planes_loop(day).DIC_planes.xa_adv >= xgrid(n_x) & planes_loop(day).DIC_planes.xa_adv < xgrid(n_x)+0.49 & ...
                planes_loop(day).DIC_planes.ya_adv >= ygrid(n_y) & planes_loop(day).DIC_planes.ya_adv < ygrid(n_y)+0.49;
            planes_loop(day).DIC_planes.gridded_DIC_var(n_bin,n_x,n_y) = nanmedian(planes_loop(day).DIC_planes.DIC_adv(check_bin));
            planes_loop(day).DIC_planes.gridded_DIC_var_std(n_bin,n_x,n_y) = nanstd(planes_loop(day).DIC_planes.DIC_adv(check_bin));
            planes_loop(day).DIC_planes.gridded_xa(n_bin,n_x,n_y) = xgrid(n_x);
            planes_loop(day).DIC_planes.gridded_ya(n_bin,n_x,n_y) = ygrid(n_y);    
            planes_loop(day).DIC_planes.gridded_z(n_bin) =  options.vertical_grid(n_bin);
        end
    end
end

% get fits
for n_bin = 1:numel(options.vertical_grid)
    xa = squeeze(planes_loop(day).DIC_planes.gridded_xa(n_bin,:,:)); 
    ya = squeeze(planes_loop(day).DIC_planes.gridded_ya(n_bin,:,:)); 
    DIC = squeeze(planes_loop(day).DIC_planes.gridded_DIC_var(n_bin,:,:));
    check_nans = isfinite(DIC);
    if numel(DIC(check_nans)) > 1
      % prepare for fit
      planes_loop(day).DIC_planes.gridded_const = ones(size(DIC(check_nans)));
      planes_loop(day).DIC_planes.gridded_plane_matrix = [xa(check_nans), ya(check_nans), planes_loop(day).DIC_planes.gridded_const];
      % robust fit using bisquare method (less weighting to outliers)
      [planes_loop(day).DIC_planes.gridded_DICfit_robust(n_bin).vals,stats] = robustfit(planes_loop(day).DIC_planes.gridded_plane_matrix(:,1:2),DIC(check_nans),'bisquare');
      planes_loop(day).DIC_planes.gridded_m_lon_error(n_bin) = stats.se(2);
      planes_loop(day).DIC_planes.gridded_m_lat_error(n_bin) = stats.se(3);
      planes_loop(day).DIC_planes.robust_stats = stats;
      planes_loop(day).DIC_planes.gridded_m_lon(n_bin) = planes_loop(day).DIC_planes.gridded_DICfit_robust(n_bin).vals(2);
      planes_loop(day).DIC_planes.gridded_m_lat(n_bin)  = planes_loop(day).DIC_planes.gridded_DICfit_robust(n_bin).vals(3);
      planes_loop(day).DIC_planes.gridded_m_0(n_bin)  = planes_loop(day).DIC_planes.gridded_DICfit_robust(n_bin).vals(1);           
    end
end

    
    %% Geopotential anomalies
     ii = 1;   
      for bn = planes_loop(day).GPA_planes.D_levels
        % get values for particular bin depth only
        bin_depth_check = planes_loop(day).GPA_planes.depth == bn;

        % get coefficients  
        planes_loop(day).GPA_planes.gpafit(ii).GPA = planes_loop(day).GPA_planes.GPA(bin_depth_check);
        planes_loop(day).GPA_planes.gpafit(ii).const = ...
            ones(size(planes_loop(day).GPA_planes.GPA(bin_depth_check)));
        planes_loop(day).GPA_planes.gpafit(ii).plane_matrix =...
            [planes_loop(day).GPA_planes.xa(bin_depth_check)', ...
            planes_loop(day).GPA_planes.ya(bin_depth_check)', ...
            planes_loop(day).GPA_planes.gpafit(ii).const'];
        % standard fit
        [planes_loop(day).GPA_planes.standard_fit(ii).fit,planes_loop(day).GPA_planes.standard_fit(ii).gof,...
             planes_loop(day).GPA_planes.standard_fit(ii).out] = fit(planes_loop(day).GPA_planes.gpafit(ii).plane_matrix(:,1:2),planes_loop(day).GPA_planes.GPA(bin_depth_check)','poly11');     
        % robust fit using bisquare method (less weighting to outliers)
        [vals,stats] = robustfit(planes_loop(day).GPA_planes.gpafit(ii).plane_matrix(:,1:2),planes_loop(day).GPA_planes.GPA(bin_depth_check)','bisquare');
        planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix = vals;
        % normal fit function
        [planes_loop(day).GPA_planes.gpafit(ii).fit, planes_loop(day).GPA_planes.gpafit(ii).gof,~]  = ...
            fit(planes_loop(day).GPA_planes.gpafit(ii).plane_matrix(:,1:2),planes_loop(day).GPA_planes.GPA(bin_depth_check)','poly11');
    
        planes_loop(day).GPA_planes.gpafit(ii).error_x = stats.se(2);
        planes_loop(day).GPA_planes.gpafit(ii).error_y = stats.se(3);            
        planes_loop(day).GPA_planes.gpafit(ii).gpa = planes_loop(day).GPA_planes.GPA;
        planes_loop(day).GPA_planes.gpa_xA(ii) = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(2);
        planes_loop(day).GPA_planes.gpa_yA(ii)  = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(3);
        planes_loop(day).GPA_planes.gpa_0A(ii)  = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(1);
        planes_loop(day).GPA_planes.gpafit(ii).gpa_predictions = ...
            (planes_loop(day).GPA_planes.xa(bin_depth_check)*planes_loop(day).GPA_planes.gpa_xA(ii)) ...
            + (planes_loop(day).GPA_planes.ya(bin_depth_check)*planes_loop(day).GPA_planes.gpa_yA(ii)) ...
            + planes_loop(day).GPA_planes.gpa_0A(ii);
        % get statistical information
        planes_loop(day).GPA_planes.gpafit(ii).stats = stats;
        
        ii = ii+1;
                    
      end
    
    %% Across transects method 
    
    %% Oxygen
    
    planes_loop(day).across_O2.S = planes_loop(day).oxygen_planes.O2_adv(...
        planes_loop(day).oxygen_planes.lat > 43.2 & planes_loop(day).oxygen_planes.lat  <=43.3);
    planes_loop(day).across_O2.N = planes_loop(day).oxygen_planes.O2_adv(...
        planes_loop(day).oxygen_planes.lat  > 43.43 & planes_loop(day).oxygen_planes.lat  < 43.5);
    planes_loop(day).across_O2.W = planes_loop(day).oxygen_planes.O2_adv(...
        planes_loop(day).oxygen_planes.lon  > 7.64 & planes_loop(day).oxygen_planes.lon  < 7.725);
    planes_loop(day).across_O2.E = planes_loop(day).oxygen_planes.O2_adv(...
        planes_loop(day).oxygen_planes.lon  > 7.925 & planes_loop(day).oxygen_planes.lon  < 8);
    
    planes_loop(day).across_O2.Sx = planes_loop(day).oxygen_planes.xa(...
        planes_loop(day).oxygen_planes.lat > 43.2 & planes_loop(day).oxygen_planes.lat  <=43.3);
    planes_loop(day).across_O2.Nx = planes_loop(day).oxygen_planes.xa(...
        planes_loop(day).oxygen_planes.lat  > 43.43 & planes_loop(day).oxygen_planes.lat  < 43.5);
    planes_loop(day).across_O2.Wx = planes_loop(day).oxygen_planes.xa(...
        planes_loop(day).oxygen_planes.lon  > 7.64 & planes_loop(day).oxygen_planes.lon  < 7.725);
    planes_loop(day).across_O2.Ex = planes_loop(day).oxygen_planes.xa(...
        planes_loop(day).oxygen_planes.lon  > 7.925 & planes_loop(day).oxygen_planes.lon  < 8);    
    
    planes_loop(day).across_O2.Sy = planes_loop(day).oxygen_planes.ya(...
        planes_loop(day).oxygen_planes.lat > 43.2 & planes_loop(day).oxygen_planes.lat  <=43.3);
    planes_loop(day).across_O2.Ny = planes_loop(day).oxygen_planes.ya(...
        planes_loop(day).oxygen_planes.lat  > 43.43 & planes_loop(day).oxygen_planes.lat  < 43.5);
    planes_loop(day).across_O2.Wy = planes_loop(day).oxygen_planes.ya(...
        planes_loop(day).oxygen_planes.lon  > 7.64 & planes_loop(day).oxygen_planes.lon  < 7.725);
    planes_loop(day).across_O2.Ey = planes_loop(day).oxygen_planes.ya(...
        planes_loop(day).oxygen_planes.lon  > 7.925 & planes_loop(day).oxygen_planes.lon  < 8);       
    
    planes_loop(day).across_O2.gradient_x_across = (nanmedian(planes_loop(day).across_O2.E) ...
        - nanmedian(planes_loop(day).across_O2.W)) / ...
        (nanmedian(planes_loop(day).across_O2.Ex) - nanmedian(planes_loop(day).across_O2.Wx));            
    
     planes_loop(day).across_O2.gradient_y_across = (nanmedian(planes_loop(day).across_O2.N) ...
         - nanmedian(planes_loop(day).across_O2.S)) / ...
         (nanmedian(planes_loop(day).across_O2.Ny) - nanmedian(planes_loop(day).across_O2.Sy));    
     
    %% DIC
    
    planes_loop(day).across_DIC.S = planes_loop(day).DIC_planes.DIC_adv(...
        planes_loop(day).DIC_planes.lat > 43.2 & planes_loop(day).DIC_planes.lat  <=43.3);
    planes_loop(day).across_DIC.N = planes_loop(day).DIC_planes.DIC_adv(...
        planes_loop(day).DIC_planes.lat  > 43.43 & planes_loop(day).DIC_planes.lat  < 43.5);
    planes_loop(day).across_DIC.W = planes_loop(day).DIC_planes.DIC_adv(...
        planes_loop(day).DIC_planes.lon  > 7.64 & planes_loop(day).DIC_planes.lon  < 7.725);
    planes_loop(day).across_DIC.E = planes_loop(day).DIC_planes.DIC_adv(...
        planes_loop(day).DIC_planes.lon  > 7.925 & planes_loop(day).DIC_planes.lon  < 8);
    
    planes_loop(day).across_DIC.Sx = planes_loop(day).DIC_planes.xa(...
        planes_loop(day).DIC_planes.lat > 43.2 & planes_loop(day).DIC_planes.lat  <=43.3);
    planes_loop(day).across_DIC.Nx = planes_loop(day).DIC_planes.xa(...
        planes_loop(day).DIC_planes.lat  > 43.43 & planes_loop(day).DIC_planes.lat  < 43.5);
    planes_loop(day).across_DIC.Wx = planes_loop(day).DIC_planes.xa(...
        planes_loop(day).DIC_planes.lon  > 7.64 & planes_loop(day).DIC_planes.lon  < 7.725);
    planes_loop(day).across_DIC.Ex = planes_loop(day).DIC_planes.xa(...
        planes_loop(day).DIC_planes.lon  > 7.925 & planes_loop(day).DIC_planes.lon  < 8);    
    
    planes_loop(day).across_DIC.Sy = planes_loop(day).DIC_planes.ya(...
        planes_loop(day).DIC_planes.lat > 43.2 & planes_loop(day).DIC_planes.lat  <=43.3);
    planes_loop(day).across_DIC.Ny = planes_loop(day).DIC_planes.ya(...
        planes_loop(day).DIC_planes.lat  > 43.43 & planes_loop(day).DIC_planes.lat  < 43.5);
    planes_loop(day).across_DIC.Wy = planes_loop(day).DIC_planes.ya(...
        planes_loop(day).DIC_planes.lon  > 7.64 & planes_loop(day).DIC_planes.lon  < 7.725);
    planes_loop(day).across_DIC.Ey = planes_loop(day).DIC_planes.ya(...
        planes_loop(day).DIC_planes.lon  > 7.925 & planes_loop(day).DIC_planes.lon  < 8);       
    
    planes_loop(day).across_DIC.gradient_x_across = (nanmedian(planes_loop(day).across_DIC.E) ...
        - nanmedian(planes_loop(day).across_DIC.W)) / ...
        (nanmedian(planes_loop(day).across_DIC.Ex) - nanmedian(planes_loop(day).across_DIC.Wx));            
    
     planes_loop(day).across_DIC.gradient_y_across = (nanmedian(planes_loop(day).across_DIC.N) ...
         - nanmedian(planes_loop(day).across_DIC.S)) / ...
         (nanmedian(planes_loop(day).across_DIC.Ny) - nanmedian(planes_loop(day).across_DIC.Sy));        

    %% get time window means
    
    % all measurements
    planes_loop(day).means.O2_h = nanmean(planes_loop(day).oxygen_planes.O2);  % mmol m^-3
    planes_loop(day).means.O2_h_median = nanmedian(planes_loop(day).oxygen_planes.O2);  % mmol m^-3
    planes_loop(day).means.O2_std_h = nanstd(planes_loop(day).oxygen_planes.O2);  % mmol m^-3 
    planes_loop(day).means.O2_standard_error = planes_loop(day).means.O2_std_h/sqrt(numel(planes_loop(day).oxygen_planes.O2));
    planes_loop(day).means.O2_MAD_h = mad(planes_loop(day).oxygen_planes.O2,1);  % mmol m^-3     
    planes_loop(day).means.O2_surf = nanmean(vars.O2(planes_loop(day).day_selection_h & vars.depth < 10));  % mmol m^-3    
    planes_loop(day).means.O2_surf_std = nanstd(planes_loop(day).day_selection_h & vars.depth < 10);  % mmol m^-3     
    planes_loop(day).means.O2_inv_h = planes_loop(day).means.O2_h*options.h;  % mmol m^-2
    planes_loop(day).means.DIC_h = nanmean(planes_loop(day).DIC_planes.DIC);  % mmol m^-3
    planes_loop(day).means.DIC_h_median = nanmedian(planes_loop(day).DIC_planes.DIC);  % mmol m^-3
    planes_loop(day).means.DIC_std_h = nanstd(planes_loop(day).DIC_planes.DIC);  % mmol m^-3 
    planes_loop(day).means.DIC_standard_error = planes_loop(day).means.DIC_std_h/sqrt(numel(planes_loop(day).DIC_planes.DIC));
    planes_loop(day).means.DIC_MAD_h = mad(planes_loop(day).DIC_planes.DIC,1);  % mmol m^-3     
    planes_loop(day).means.DIC_surf = nanmean(vars.DIC(planes_loop(day).day_selection_h & vars.depth < 10));  % mmol m^-3    
    planes_loop(day).means.DIC_surf_std = nanstd(planes_loop(day).day_selection_h & vars.depth < 10);  % mmol m^-3     
    planes_loop(day).means.DIC_inv_h = planes_loop(day).means.DIC_h*options.h;  % mmol m^-2    
    planes_loop(day).means.DACu_h = nanmean(vars.DACs.DACu(planes_loop(day).day_selection_DAC));   
    planes_loop(day).means.DACv_h = nanmean(vars.DACs.DACv(planes_loop(day).day_selection_DAC));  
    planes_loop(day).means.DACu_std_h = nanstd(vars.DACs.DACu(planes_loop(day).day_selection_DAC));   
    planes_loop(day).means.DACv_std_h = nanstd(vars.DACs.DACv(planes_loop(day).day_selection_DAC));    
    planes_loop(day).means.lon_h = nanmedian(vars.lon(planes_loop(day).day_selection_h));
    planes_loop(day).means.lat_h = nanmedian(vars.lat(planes_loop(day).day_selection_h));
    planes_loop(day).means.fCO2_h =nanmean(vars.fCO2(planes_loop(day).day_selection_h));    
    planes_loop(day).means.fCO213_h =nanmean(vars.fCO213(planes_loop(day).day_selection_h));
    planes_loop(day).means.fCO2_surf =nanmean(vars.fCO2(planes_loop(day).day_selection_h & vars.depth < 10));    
    planes_loop(day).means.fCO213_surf =nanmean(vars.fCO213(planes_loop(day).day_selection_h & vars.depth < 10));    
    planes_loop(day).means.P_h = nanmean(vars.P(planes_loop(day).day_selection_h));
    planes_loop(day).means.MLD_h = nanmean(vars.MLD_combined(planes_loop(day).day_selection_h));
    planes_loop(day).means.MLD_std_h = nanstd(vars.MLD_combined(planes_loop(day).day_selection_h)); 
    planes_loop(day).means.S_h = nanmean(vars.S(planes_loop(day).day_selection_h));
    planes_loop(day).means.S_std_h = nanstd(vars.S(planes_loop(day).day_selection_h));
    planes_loop(day).means.S_surf = nanmean(vars.S(planes_loop(day).day_selection_h & vars.depth < 10));
    planes_loop(day).means.S_surf_std = nanstd(vars.S(planes_loop(day).day_selection_h & vars.depth < 10));   
    planes_loop(day).means.T_h = nanmean(vars.T(planes_loop(day).day_selection_h));
    planes_loop(day).means.T_std_h = nanstd(vars.T(planes_loop(day).day_selection_h));
    planes_loop(day).means.T_surf = nanmean(vars.T(planes_loop(day).day_selection_h & vars.depth < 10));
    planes_loop(day).means.T_surf_std = nanstd(vars.T(planes_loop(day).day_selection_h & vars.depth < 10));   
    planes_loop(day).means.sig0_h = nanmean(vars.sigma0(planes_loop(day).day_selection_h));
    planes_loop(day).means.sig0_surf = nanmean(vars.sigma0(planes_loop(day).day_selection_h & vars.depth < 10));    
    planes_loop(day).means.wind = nanmean(vars.Wind(planes_loop(day).day_selection_h));
    planes_loop(day).means.wind_squared = nanmean(planes_loop(day).means.wind.^2);       
    % Satellite currents
    % AVISO
    check_sat_time =  AVISO.time >= datenum(2016,03,day-options.window*options.window_adv,00,00,00) ...
        &  AVISO.time <= datenum(2016,03,day+options.window*options.window_adv,00,00,00);
    planes_loop(day).means.GVsat_U = nanmean(AVISO.absU_middle(check_sat_time));
    planes_loop(day).means.GVsat_V = nanmean(AVISO.absV_middle(check_sat_time));    
    % Oscar
    check_sat_time =  Oscar.time >= datenum(2016,03,day-options.window*options.window_adv,00,00,00) ...
        &  Oscar.time <= datenum(2016,03,day+options.window*options.window_adv,00,00,00);    
    planes_loop(day).means.GVsat_Oscar_U = nanmean(Oscar.U_middle(check_sat_time));
    planes_loop(day).means.GVsat_Oscar_V = nanmean(Oscar.V_middle(check_sat_time));  
    % sea level pressure
    planes_loop(day).means.sea_level_pressure = nanmean(vars.sea_level_press(planes_loop(day).day_selection_h))/100; % in mbar
    planes_loop(day).means.sea_level_pressure_atm = planes_loop(day).means.sea_level_pressure/1013.25; % in atm
    % profs
    planes_loop(day).means.O2_h_profs = nanmean(planes_loop(day).oxygen_planes.O2_profs);  % mmol m^-3
    planes_loop(day).means.O2_h_median_profs = nanmedian(planes_loop(day).oxygen_planes.O2_profs);  % mmol m^-3
    planes_loop(day).means.O2_std_h_profs = nanstd(planes_loop(day).oxygen_planes.O2_profs);  % mmol m^-3 
    planes_loop(day).means.O2_standard_error_profs = planes_loop(day).means.O2_std_h_profs/sqrt(numel(planes_loop(day).oxygen_planes.O2_profs));
    planes_loop(day).means.O2_MAD_h_profs = mad(planes_loop(day).oxygen_planes.O2_profs,1);  % mmol m^-3     
    planes_loop(day).means.O2_surf_profs = nanmean(profs.O2_surf(planes_loop(day).day_selection_h_profs));  % mmol m^-3    
    planes_loop(day).means.O2_surf_std_profs = nanstd(profs.O2_surf(planes_loop(day).day_selection_h_profs));   % mmol m^-3     
    planes_loop(day).means.O2_inv_h_profs = planes_loop(day).means.O2_h_profs*options.h;  % mmol m^-2
    planes_loop(day).means.DIC_h_profs = nanmean(planes_loop(day).DIC_planes.DIC_profs);  % mmol m^-3
    planes_loop(day).means.DIC_h_median_profs = nanmedian(planes_loop(day).DIC_planes.DIC_profs);  % mmol m^-3
    planes_loop(day).means.DIC_std_h_profs = nanstd(planes_loop(day).DIC_planes.DIC_profs);  % mmol m^-3 
    planes_loop(day).means.DIC_standard_error_profs = planes_loop(day).means.DIC_std_h_profs/sqrt(numel(planes_loop(day).DIC_planes.DIC_profs));
    planes_loop(day).means.DIC_MAD_h_profs = mad(planes_loop(day).DIC_planes.DIC_profs,1);  % mmol m^-3     
    planes_loop(day).means.DIC_surf_profs = nanmean(profs.DIC_surf(planes_loop(day).day_selection_h_profs));  % mmol m^-3    
    planes_loop(day).means.DIC_surf_std_profs = nanstd(profs.DIC_surf(planes_loop(day).day_selection_h_profs));   % mmol m^-3     
    planes_loop(day).means.DIC_inv_h_profs = planes_loop(day).means.DIC_h_profs*options.h;  % mmol m^-2    
    planes_loop(day).means.lon_h_profs = nanmedian(profs.lon(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.lat_h_profs = nanmedian(profs.lat(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.fCO2_h_profs =nanmean(profs.fCO2(planes_loop(day).day_selection_h_profs));    
    planes_loop(day).means.fCO213_h_profs =nanmean(profs.fCO213(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.fCO2_surf_profs =nanmean(profs.fCO2_surf(planes_loop(day).day_selection_h_profs));    
    planes_loop(day).means.fCO213_surf_profs =nanmean(profs.fCO213_surf(planes_loop(day).day_selection_h_profs));    
    planes_loop(day).means.P_h_profs = nanmean(profs.P(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.MLD_h_profs = nanmean(profs.MLD_combined(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.MLD_std_h_profs = nanstd(profs.MLD_combined(planes_loop(day).day_selection_h_profs)); 
    planes_loop(day).means.S_h_profs = nanmean(profs.S(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.S_std_h_profs = nanstd(profs.S(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.S_surf_profs = nanmean(profs.S_surf(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.S_surf_std_profs = nanstd(profs.S_surf(planes_loop(day).day_selection_h_profs));   
    planes_loop(day).means.T_h_profs = nanmean(profs.T(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.T_std_h_profs = nanstd(profs.T(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.T_surf_profs = nanmean(profs.T_surf(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.T_surf_std_profs = nanstd(profs.T_surf(planes_loop(day).day_selection_h_profs));   
    planes_loop(day).means.sig0_h_profs = nanmean(profs.sigma0(planes_loop(day).day_selection_h_profs));
    planes_loop(day).means.sig0_surf_profs = nanmean(profs.sigma0_surf(planes_loop(day).day_selection_h_profs));    
    
end

%% Loop for entrainment variables


for day = options.dayrange(2):options.interval:options.dayrange(end)-options.interval

    
    %% data selection logical vector to include data 2 days either side of day
    % different vectors because of different dimensions
    
    if day <= 31
        date_num = datenum(2016,03,day,00,00,00);
    else
        date_num = datenum(2016,04,day-31,00,00,00);
    end

    tdayPUmb = vars.t > planes_loop(day).date_num - options.window ...
        & vars.t < planes_loop(day).date_num + options.window ...
        & vars.depth < planes_loop(day+1).means.MLD_h;  % only top h metres
    tdayPUmb_profs = profs.t > planes_loop(day).date_num - options.window ...
        & profs.t < planes_loop(day).date_num + options.window ...
        & profs.depth < planes_loop(day+1).means.MLD_h;  % only top h metres
    
    O2_dayUmb = vars.O2(tdayPUmb);    
    O2_dayUmb_profs = profs.O2(tdayPUmb_profs); 
    DIC_dayUmb = vars.DIC(tdayPUmb);    
    DIC_dayUmb_profs = profs.DIC(tdayPUmb_profs);       
    
    planes_loop(day).means.O2invt1MLDt2 = nanmean(O2_dayUmb)*planes_loop(day+1).means.MLD_h;
    planes_loop(day).means.O2invt1MLDt2std = nanstd(O2_dayUmb)*planes_loop(day+1).means.MLD_h;
    planes_loop(day).means.O2invt1MLDt2_profs = nanmean(O2_dayUmb_profs)*planes_loop(day+1).means.MLD_h_profs;
    planes_loop(day).means.O2invt1MLDt2std_profs = nanstd(O2_dayUmb_profs)*planes_loop(day+1).means.MLD_h_profs;    

    planes_loop(day).means.DICinvt1MLDt2 = nanmean(DIC_dayUmb)*planes_loop(day+1).means.MLD_h;
    planes_loop(day).means.DICinvt1MLDt2std = nanstd(DIC_dayUmb)*planes_loop(day+1).means.MLD_h;
    planes_loop(day).means.DICinvt1MLDt2_profs = nanmean(DIC_dayUmb_profs)*planes_loop(day+1).means.MLD_h_profs;
    planes_loop(day).means.DICinvt1MLDt2std_profs = nanstd(DIC_dayUmb_profs)*planes_loop(day+1).means.MLD_h_profs;    
    
end
% need to add NaN dummies so that horzcat works with structure 
for day = [options.dayrange(1) options.dayrange(end)]
    planes_loop(day).means.O2invt1MLDt2 = NaN;
    planes_loop(day).means.O2invt1MLDt2std = NaN;
    planes_loop(day).means.O2invt1MLDt2_profs = NaN;
    planes_loop(day).means.O2invt1MLDt2std_profs = NaN;    

    planes_loop(day).means.DICinvt1MLDt2 = NaN;
    planes_loop(day).means.DICinvt1MLDt2std = NaN;
    planes_loop(day).means.DICinvt1MLDt2_profs = NaN;
    planes_loop(day).means.DICinvt1MLDt2std_profs = NaN;        
    
end

clear O2_dayUmb *dayPUmb ii bin_depth_check bn check* date_num day GPA GPA_* lat lon press_selection *_bin
%% save file for later use
save([options.data_dir,'planes_',num2str(options.h),'.mat'],'planes_loop');

%% if file already created
else
    disp(['O2 Plane-fits | Loading previously created file containing O2 Plane-fits']);
    load([options.data_dir,'planes_',num2str(options.h),'.mat']);
end
    means_struct = [planes_loop.means];
    disp(['O2 Plane-fits | Finished']);