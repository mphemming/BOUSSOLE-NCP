%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_planes.m

% Script to calculate planes in oxygen and geopotential anomalies 
% needed to calculate the advection term

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 21/06/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if saved file already, don't need to do all of below
testing =1;
if sum(~exist('planes.mat','file') + testing == 1) ~= 0

%% get concatenated geopotential anomaly profile arrays for planes

disp(['Plane-fits | organising data']);

GPA = [vars.bin.GPA];    
GPA_D = [vars.bin.depth];    
initial = [vars.bin.t]; GPA_t = [initial.median_var];
initial = [vars.bin.lon]; lon = [initial.median_var];
initial = [vars.bin.lat]; lat = [initial.median_var];

clear initial
%% calculate cartesian coordinates for plane-fits
% Compute x (longitude) & y (latitude) distances relative to the origin:

disp(['Plane-fits | Calculating cartesian coordinates']);

% chosen lat/lon outside of domain
lat0 = 43.2;
lon0 = 7.5;

% For O2
for i = 1:numel(vars.P)
    [vars.xa(i),~] = sw_dist([lat0; lat0],[lon0; vars.lon(i)],'km');
    [vars.ya(i),~] = sw_dist([lat0; vars.lat(i)],[lon0; lon0],'km');
end

% For GPA
for i = 1:numel(lon);
    [xa_bin(i),~] = sw_dist([lat0; lat0],[lon0; lon(i)],'km');
    [ya_bin(i),~] = sw_dist([lat0; lat(i)],[lon0; lon0],'km');
end    

clear lat0 lon0 i
%% Loop | for each day +- n days moving window, 
% obtain plane-fits of oxygen and geopotential anomalies

disp(['Plane-fits | Obtaining gradients for each time step window']);

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
    planes_loop(day).day_selection_GPA = ...
        GPA_t > (planes_loop(day).date_num-options.window) ...
        & GPA_t < (planes_loop(day).date_num+options.window); % for binned GPA profiles
    planes_loop(day).day_selection_DAC = ...
        vars.DACs.t > (planes_loop(day).date_num-options.window) & ...
        vars.DACs.t < (planes_loop(day).date_num+options.window); % for DACs
    % Select data for oxygen plane-fits with NaNs removed
    check_nan = isfinite(vars.O2) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).oxygen_planes.O2 = vars.O2(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.xa = vars.xa(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.ya = vars.ya(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.P = vars.P(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.lon = vars.lon(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).oxygen_planes.lat = vars.lat(planes_loop(day).day_selection_h & check_nan);
    % Select data for geopotential anomaly plane-fits with NaNs removed
    check_nan = isfinite(GPA) & isfinite(xa_bin) & isfinite(ya_bin);
    planes_loop(day).GPA_planes.GPA = GPA(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.depth = GPA_D(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.time = GPA_t(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.xa = xa_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.ya = ya_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.D_levels = unique(planes_loop(day).GPA_planes.depth);
    % get smoothed oxygen
    xgrid = 10:0.05:45; ygrid = 0:0.05:40;
    planes_loop(day).oxygen_planes.O2_bin = bindata2(planes_loop(day).oxygen_planes.O2,...
        planes_loop(day).oxygen_planes.xa,planes_loop(day).oxygen_planes.ya,xgrid,ygrid);
    NaN_check = reshape(isnan(planes_loop(day).oxygen_planes.O2_bin),size(planes_loop(day).oxygen_planes.O2_bin));
    planes_loop(day).oxygen_planes.O2_bin_smooth = smoothn(planes_loop(day).oxygen_planes.O2_bin,20/0.05);   
    planes_loop(day).oxygen_planes.O2_bin_smooth(NaN_check) = NaN;
    planes_loop(day).oxygen_planes.xa_smooth = repmat(xgrid(1:end-1),[800 1]);
    planes_loop(day).oxygen_planes.ya_smooth = repmat(ygrid(1:end-1),[700 1])';
    planes_loop(day).oxygen_planes.O2_bin_smooth = planes_loop(day).oxygen_planes.O2_bin_smooth(:);
    planes_loop(day).oxygen_planes.xa_smooth = planes_loop(day).oxygen_planes.xa_smooth(:);
    planes_loop(day).oxygen_planes.ya_smooth = planes_loop(day).oxygen_planes.ya_smooth(:);
    check = isnan(planes_loop(day).oxygen_planes.O2_bin_smooth);
    planes_loop(day).oxygen_planes.O2_bin_smooth(check) = [];
    planes_loop(day).oxygen_planes.xa_smooth(check) = [];
    planes_loop(day).oxygen_planes.ya_smooth(check) = [];
    % get lat and lon info
    planes_loop(day).oxygen_planes.lat_bin = bindata2(planes_loop(day).oxygen_planes.lat,...
        planes_loop(day).oxygen_planes.xa,planes_loop(day).oxygen_planes.ya,xgrid,ygrid);
    planes_loop(day).oxygen_planes.lon_bin = bindata2(planes_loop(day).oxygen_planes.lon,...
        planes_loop(day).oxygen_planes.xa,planes_loop(day).oxygen_planes.ya,xgrid,ygrid);    
    planes_loop(day).oxygen_planes.lat_bin = planes_loop(day).oxygen_planes.lat_bin(:);
    planes_loop(day).oxygen_planes.lon_bin = planes_loop(day).oxygen_planes.lon_bin(:);
    planes_loop(day).oxygen_planes.lat_bin(check) = [];
    planes_loop(day).oxygen_planes.lon_bin(check) = [];
    
    clear check NaN_check xgrid ygrid
    
    %% Plane-fits
    
    %% Oxygen
    
    planes_loop(day).oxygen_planes.const = ones(size(planes_loop(day).oxygen_planes.O2_bin_smooth));
    planes_loop(day).oxygen_planes.plane_matrix = [planes_loop(day).oxygen_planes.xa_smooth, ...
        planes_loop(day).oxygen_planes.ya_smooth, planes_loop(day).oxygen_planes.const];
    planes_loop(day).oxygen_planes.coefficients_from_matrix = planes_loop(day).oxygen_planes.plane_matrix\planes_loop(day).oxygen_planes.O2_bin_smooth;
    % fit function below gets same result as line above. 
    [planes_loop(day).oxygen_planes.O2fit, planes_loop(day).oxygen_planes.gof,~]  = fit(planes_loop(day).oxygen_planes.plane_matrix(:,1:2),...
        planes_loop(day).oxygen_planes.O2_bin_smooth,'poly11');
    % get important variables
    planes_loop(day).oxygen_planes.day = day;
    planes_loop(day).oxygen_planes.datestr = datestr(planes_loop(day).date_num); 
    % get confidence boundaries (95%) for fit
    planes_loop(day).oxygen_planes.confidence_boundaries = confint(planes_loop(day).oxygen_planes.O2fit);
    planes_loop(day).oxygen_planes.error_x = abs(planes_loop(day).oxygen_planes.O2fit.p10 - planes_loop(day).oxygen_planes.confidence_boundaries(1,2));      
    planes_loop(day).oxygen_planes.error_y = abs(planes_loop(day).oxygen_planes.O2fit.p01 - planes_loop(day).oxygen_planes.confidence_boundaries(1,3));         
    % get statistical information
    planes_loop(day).oxygen_planes.m_lon = planes_loop(day).oxygen_planes.O2fit.p10;
    planes_loop(day).oxygen_planes.m_lat  = planes_loop(day).oxygen_planes.O2fit.p01;
    planes_loop(day).oxygen_planes.m_0  = planes_loop(day).oxygen_planes.O2fit.p00;
    planes_loop(day).oxygen_planes.O2_predictions = (planes_loop(day).oxygen_planes.xa*planes_loop(day).oxygen_planes.m_lon) + ...
        (planes_loop(day).oxygen_planes.ya*planes_loop(day).oxygen_planes.m_lat) + planes_loop(day).oxygen_planes.m_0;
    planes_loop(day).oxygen_planes.O2_resid = planes_loop(day).oxygen_planes.O2_bin_smooth - planes_loop(day).oxygen_planes.O2_predictions;
    planes_loop(day).oxygen_planes.O2_SSR = sum(planes_loop(day).oxygen_planes.O2_resid.^2);
    planes_loop(day).oxygen_planes.O2_SSE =  planes_loop(day).oxygen_planes.O2_predictions - mean(planes_loop(day).oxygen_planes.O2_bin_smooth); 
    planes_loop(day).oxygen_planes.O2_SSE =  planes_loop(day).oxygen_planes.O2_SSE.^2; 
    planes_loop(day).oxygen_planes.O2_SSE = sum(planes_loop(day).oxygen_planes.O2_SSE);
    planes_loop(day).oxygen_planes.O2_SST = planes_loop(day).oxygen_planes.O2_SSE + planes_loop(day).oxygen_planes.O2_SSR;
    %planes_loop(day).oxygen_planes.O2_R = planes_loop(day).oxygen_planes.O2_SSE/planes_loop(day).oxygen_planes.O2_SST;    
    planes_loop(day).oxygen_planes.O2_RMSE = planes_loop(day).oxygen_planes.gof.rmse;
    % get standard errors (may be incorrect)
    planes_loop(day).oxygen_planes.m_lon_standarderror = sqrt( nansum( ...
        (planes_loop(day).oxygen_planes.O2_bin_smooth' - (planes_loop(day).oxygen_planes.m_0 + ...
        planes_loop(day).oxygen_planes.m_lon.*planes_loop(day).oxygen_planes.plane_matrix(:,1) )).^2) ) ...
        ./ length(planes_loop(day).oxygen_planes.O2_bin_smooth);    
    planes_loop(day).oxygen_planes.m_lat_standarderror = sqrt( nansum( ...
        (planes_loop(day).oxygen_planes.O2_bin_smooth' - (planes_loop(day).oxygen_planes.m_0 + ...
        planes_loop(day).oxygen_planes.m_lat.*planes_loop(day).oxygen_planes.plane_matrix(:,2) )).^2) ) ...
        ./ length(planes_loop(day).oxygen_planes.O2_bin_smooth) ;    
    planes_loop(day).oxygen_planes.plane_standarderror = sqrt( nansum( ...
        (planes_loop(day).oxygen_planes.O2_bin_smooth' - planes_loop(day).oxygen_planes.O2fit(planes_loop(day).oxygen_planes.plane_matrix(:,1:2))).^2 ) ...
        ./ length(planes_loop(day).oxygen_planes.O2_bin_smooth) ) ; 
    
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
        planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix = ...
            planes_loop(day).GPA_planes.gpafit(ii).plane_matrix\planes_loop(day).GPA_planes.GPA(bin_depth_check)';
        % fit function below gets same result as line above. 
        [planes_loop(day).GPA_planes.gpafit(ii).fit, planes_loop(day).GPA_planes.gpafit(ii).gof,~]  = ...
            fit(planes_loop(day).GPA_planes.gpafit(ii).plane_matrix(:,1:2),planes_loop(day).GPA_planes.GPA(bin_depth_check)','poly11');
      
        planes_loop(day).GPA_planes.gpafit(ii).confidence_boundaries = ...
            confint(planes_loop(day).GPA_planes.gpafit(ii).fit);
        planes_loop(day).GPA_planes.gpafit(ii).error_x = ...
            abs(planes_loop(day).GPA_planes.gpafit(ii).fit.p10 - planes_loop(day).GPA_planes.gpafit(ii).confidence_boundaries(1,2));      
        planes_loop(day).GPA_planes.gpafit(ii).error_y = ...
            abs(planes_loop(day).GPA_planes.gpafit(ii).fit.p01 - planes_loop(day).GPA_planes.gpafit(ii).confidence_boundaries(1,3));                
        planes_loop(day).GPA_planes.gpafit(ii).gpa = planes_loop(day).GPA_planes.GPA;
        planes_loop(day).GPA_planes.gpa_xA(ii) = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(1);
        planes_loop(day).GPA_planes.gpa_yA(ii)  = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(2);
        planes_loop(day).GPA_planes.gpa_0A(ii)  = planes_loop(day).GPA_planes.gpafit(ii).coefficients_from_matrix(3);
        planes_loop(day).GPA_planes.gpafit(ii).gpa_predictions = ...
            (planes_loop(day).GPA_planes.xa(bin_depth_check)*planes_loop(day).GPA_planes.gpa_xA(ii)) ...
            + (planes_loop(day).GPA_planes.ya(bin_depth_check)*planes_loop(day).GPA_planes.gpa_yA(ii)) ...
            + planes_loop(day).GPA_planes.gpa_0A(ii);
        % get statistical information
           planes_loop(day).GPA_planes.gpafit(ii).GPA_resid = ...
                planes_loop(day).GPA_planes.gpafit(ii).GPA -  ...
                planes_loop(day).GPA_planes.gpafit(ii).gpa_predictions;
            planes_loop(day).GPA_planes.gpafit(ii).GPA_SSR = ...
                sum(planes_loop(day).GPA_planes.gpafit(ii).GPA_resid.^2);
            planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE = ...
                planes_loop(day).GPA_planes.gpafit(ii).gpa_predictions - ...
                mean(planes_loop(day).GPA_planes.gpafit(ii).GPA); 
            planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE =  planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE.^2; 
            planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE = sum( planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE);
            planes_loop(day).GPA_planes.gpafit(ii).GPA_SST = planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE + planes_loop(day).GPA_planes.gpafit(ii).GPA_SSR;
            planes_loop(day).GPA_planes.gpafit(ii).GPA_R = planes_loop(day).GPA_planes.gpafit(ii).GPA_SSE/planes_loop(day).GPA_planes.gpafit(ii).GPA_SST;
            planes_loop(day).GPA_planes.gpafit(ii).GPA_RMSE = planes_loop(day).GPA_planes.gpafit(ii).gof.rmse;        
            planes_loop(day).GPA_planes.gpafit(ii).GPA_D = bn;
            
                    ii = ii+1;
                    
      end
    
    %% Across transects method 
    
    %% Oxygen
    
    planes_loop(day).across_O2.S = planes_loop(day).oxygen_planes.O2_bin_smooth(...
        planes_loop(day).oxygen_planes.lat_bin > 43.2 & planes_loop(day).oxygen_planes.lat_bin  <=43.3);
    planes_loop(day).across_O2.N = planes_loop(day).oxygen_planes.O2_bin_smooth(...
        planes_loop(day).oxygen_planes.lat_bin  > 43.43 & planes_loop(day).oxygen_planes.lat_bin  < 43.5);
    planes_loop(day).across_O2.W = planes_loop(day).oxygen_planes.O2_bin_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.64 & planes_loop(day).oxygen_planes.lon_bin  < 7.725);
    planes_loop(day).across_O2.E = planes_loop(day).oxygen_planes.O2_bin_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.925 & planes_loop(day).oxygen_planes.lon_bin  < 8);
    
    planes_loop(day).across_O2.Sx = planes_loop(day).oxygen_planes.xa_smooth(...
        planes_loop(day).oxygen_planes.lat_bin > 43.2 & planes_loop(day).oxygen_planes.lat_bin  <=43.3);
    planes_loop(day).across_O2.Nx = planes_loop(day).oxygen_planes.xa_smooth(...
        planes_loop(day).oxygen_planes.lat_bin  > 43.43 & planes_loop(day).oxygen_planes.lat_bin  < 43.5);
    planes_loop(day).across_O2.Wx = planes_loop(day).oxygen_planes.xa_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.64 & planes_loop(day).oxygen_planes.lon_bin  < 7.725);
    planes_loop(day).across_O2.Ex = planes_loop(day).oxygen_planes.xa_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.925 & planes_loop(day).oxygen_planes.lon_bin  < 8);    
    
    planes_loop(day).across_O2.Sy = planes_loop(day).oxygen_planes.ya_smooth(...
        planes_loop(day).oxygen_planes.lat_bin > 43.2 & planes_loop(day).oxygen_planes.lat_bin  <=43.3);
    planes_loop(day).across_O2.Ny = planes_loop(day).oxygen_planes.ya_smooth(...
        planes_loop(day).oxygen_planes.lat_bin  > 43.43 & planes_loop(day).oxygen_planes.lat_bin  < 43.5);
    planes_loop(day).across_O2.Wy = planes_loop(day).oxygen_planes.ya_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.64 & planes_loop(day).oxygen_planes.lon_bin  < 7.725);
    planes_loop(day).across_O2.Ey = planes_loop(day).oxygen_planes.ya_smooth(...
        planes_loop(day).oxygen_planes.lon_bin  > 7.925 & planes_loop(day).oxygen_planes.lon_bin  < 8);       
    
    planes_loop(day).across_O2.gradient_x_across = (nanmedian(planes_loop(day).across_O2.E) ...
        - nanmedian(planes_loop(day).across_O2.W)) / ...
        (nanmedian(planes_loop(day).across_O2.Ex) - nanmedian(planes_loop(day).across_O2.Wx));            
    
     planes_loop(day).across_O2.gradient_y_across = (nanmedian(planes_loop(day).across_O2.N) ...
         - nanmedian(planes_loop(day).across_O2.S)) / ...
         (nanmedian(planes_loop(day).across_O2.Ny) - nanmedian(planes_loop(day).across_O2.Sy));    

    %% get time window means
    planes_loop(day).means.O2_h = nanmean(planes_loop(day).oxygen_planes.O2);  % mmol m^-3    
    planes_loop(day).means.O2_std_h = nanstd(planes_loop(day).oxygen_planes.O2);  % mmol m^-3 
    planes_loop(day).means.O2_surf = nanmean(vars.O2(planes_loop(day).day_selection_h & vars.depth < 10));  % mmol m^-3    
    planes_loop(day).means.O2_surf_std = nanstd(planes_loop(day).day_selection_h & vars.depth < 10);  % mmol m^-3     
    planes_loop(day).means.O2_inv_h = planes_loop(day).means.O2_h*options.h;  % mmol m^-2
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
    planes_loop(day).means.MLD_h = nanmean(vars.MLD(planes_loop(day).day_selection_h));
    planes_loop(day).means.MLD_std_h = nanstd(vars.MLD(planes_loop(day).day_selection_h)); 
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
    check_sat_time =  vars.GVsat_time >= datenum(2016,03,day-options.window,00,00,00) ...
        &  vars.GVsat_time <= datenum(2016,03,day+options.window,00,00,00);
    planes_loop(day).means.GVsat_U = nanmean(vars.GVsat_U(check_sat_time));
    planes_loop(day).means.GVsat_V = nanmean(vars.GVsat_V(check_sat_time));    
    planes_loop(day).means.sea_level_pressure = nanmean(vars.sea_level_press(planes_loop(day).day_selection_h))/100; % in mbar
    planes_loop(day).means.sea_level_pressure_atm = planes_loop(day).means.sea_level_pressure/1013.25; % in atm
    
end

%% Loop for entrainment variables


for day = options.dayrange(2):options.interval:options.dayrange(end)-options.interval

    
    %% data selection logical vector to include data 2 days either side of day
    % different vectors because of different dimensions
    
    if day <= 31;
    date_num = datenum(2016,03,day,00,00,00);
    else
    date_num = datenum(2016,04,day-31,00,00,00);
    end

tdayPUmb = vars.t > planes_loop(day).date_num - options.window ...
    & vars.t < planes_loop(day).date_num + options.window ...
    & vars.depth < planes_loop(day+1).means.MLD_h;  % only top h metres
O2_dayUmb = vars.O2(tdayPUmb);    
planes_loop(day).means.O2invt1MLDt2 = nanmean(O2_dayUmb)*planes_loop(day+1).means.MLD_h;
planes_loop(day).means.O2invt1MLDt2std = nanstd(O2_dayUmb)*planes_loop(day+1).means.MLD_h;

end
% need to add NaN dummies so that horzcat works with structure 
for day = [options.dayrange(1) options.dayrange(end)]
    planes_loop(day).means.O2invt1MLDt2 = NaN;
    planes_loop(day).means.O2invt1MLDt2std = NaN;
end

clear O2_dayUmb *dayPUmb ii bin_depth_check bn check* date_num day GPA GPA_* lat lon press_selection *_bin
%% save file for later use
save('planes','planes_loop');

%% if file already created
else
    disp(['Plane-fits | Loading previously created file containing plane-fits']);
    load('planes.mat');
end
    means_struct = [planes_loop.means];
    disp(['Plane-fits | Finished']);