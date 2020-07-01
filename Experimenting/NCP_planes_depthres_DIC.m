%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_planes.m

% Script to calculate planes in DIC and geopotential anomalies 
% needed to calculate the advection term

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 21/06/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% if saved file already, don't need to do all of below
% testing =1;
% if sum(~exist('planes.mat','file') + testing == 1) ~= 0

%% get concatenated geopotential anomaly profile arrays for planes

disp(['DIC Plane-fits | organising data']);

GPA = [vars.bin.GPA];    
GPA_D = [vars.bin.depth];    
initial = [vars.bin.t]; GPA_t = [initial.median_var];
initial = [vars.bin.lon]; lon = [initial.median_var];
initial = [vars.bin.lat]; lat = [initial.median_var];

clear initial
%% calculate cartesian coordinates for DIC Plane-fits
% Compute x (longitude) & y (latitude) distances relative to the origin:

disp(['DIC Plane-fits | Calculating cartesian coordinates']);

% chosen lat/lon outside of domain
lat0 = 43.2;
lon0 = 7.5;

% For DIC
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
% obtain DIC Plane-fits of DIC and geopotential anomalies

disp(['DIC Plane-fits | Obtaining gradients for each time step window']);

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
    % Select data for DIC DIC Plane-fits with NaNs removed
    check_nan = isfinite(vars.DIC) & isfinite(vars.xa) & isfinite(vars.ya);
    planes_loop(day).DIC_planes.DIC = vars.DIC(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.xa = vars.xa(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.ya = vars.ya(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.P = vars.P(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.lon = vars.lon(planes_loop(day).day_selection_h & check_nan);
    planes_loop(day).DIC_planes.lat = vars.lat(planes_loop(day).day_selection_h & check_nan);
    % Select data for geopotential anomaly DIC Plane-fits with NaNs removed
    check_nan = isfinite(GPA) & isfinite(xa_bin) & isfinite(ya_bin);
    planes_loop(day).GPA_planes.GPA = GPA(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.depth = GPA_D(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.time = GPA_t(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.xa = xa_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.ya = ya_bin(planes_loop(day).day_selection_GPA & check_nan);  
    planes_loop(day).GPA_planes.D_levels = unique(planes_loop(day).GPA_planes.depth);
    
    
    %% DIC Plane-fits
    
    %% DIC
    
    planes_loop(day).DIC_planes.const = ones(size(planes_loop(day).DIC_planes.DIC));
    planes_loop(day).DIC_planes.plane_matrix = [planes_loop(day).DIC_planes.xa', planes_loop(day).DIC_planes.ya', planes_loop(day).DIC_planes.const'];
    planes_loop(day).DIC_planes.coefficients_from_matrix = planes_loop(day).DIC_planes.plane_matrix\planes_loop(day).DIC_planes.DIC';
    % fit function below gets same result as line above. 
    [planes_loop(day).DIC_planes.DICfit, planes_loop(day).DIC_planes.gof,~]  = fit(planes_loop(day).DIC_planes.plane_matrix(:,1:2),planes_loop(day).DIC_planes.DIC','poly11');
    % get important variables
    planes_loop(day).DIC_planes.day = day;
    planes_loop(day).DIC_planes.datestr = datestr(planes_loop(day).date_num); 
    % get confidence boundaries (95%) for fit
    planes_loop(day).DIC_planes.confidence_boundaries = confint(planes_loop(day).DIC_planes.DICfit);
    planes_loop(day).DIC_planes.error_x = abs(planes_loop(day).DIC_planes.DICfit.p10 - planes_loop(day).DIC_planes.confidence_boundaries(1,2));      
    planes_loop(day).DIC_planes.error_y = abs(planes_loop(day).DIC_planes.DICfit.p01 - planes_loop(day).DIC_planes.confidence_boundaries(1,3));         
    % get statistical information
    planes_loop(day).DIC_planes.m_lon = planes_loop(day).DIC_planes.DICfit.p10;
    planes_loop(day).DIC_planes.m_lat  = planes_loop(day).DIC_planes.DICfit.p01;
    planes_loop(day).DIC_planes.m_0  = planes_loop(day).DIC_planes.DICfit.p00;
    planes_loop(day).DIC_planes.DIC_predictions = (planes_loop(day).DIC_planes.xa*planes_loop(day).DIC_planes.m_lon) + ...
        (planes_loop(day).DIC_planes.ya*planes_loop(day).DIC_planes.m_lat) + planes_loop(day).DIC_planes.m_0;
    planes_loop(day).DIC_planes.DIC_resid = planes_loop(day).DIC_planes.DIC - planes_loop(day).DIC_planes.DIC_predictions;
    planes_loop(day).DIC_planes.DIC_SSR = sum(planes_loop(day).DIC_planes.DIC_resid.^2);
    planes_loop(day).DIC_planes.DIC_SSE =  planes_loop(day).DIC_planes.DIC_predictions - mean(planes_loop(day).DIC_planes.DIC); 
    planes_loop(day).DIC_planes.DIC_SSE =  planes_loop(day).DIC_planes.DIC_SSE.^2; 
    planes_loop(day).DIC_planes.DIC_SSE = sum(planes_loop(day).DIC_planes.DIC_SSE);
    planes_loop(day).DIC_planes.DIC_SST = planes_loop(day).DIC_planes.DIC_SSE + planes_loop(day).DIC_planes.DIC_SSR;
    planes_loop(day).DIC_planes.DIC_R = planes_loop(day).DIC_planes.DIC_SSE/planes_loop(day).DIC_planes.DIC_SST;    
    planes_loop(day).DIC_planes.DIC_RMSE = planes_loop(day).DIC_planes.gof.rmse;
    % get standard errors (may be incorrect)
    planes_loop(day).DIC_planes.m_lon_standarderror = sqrt( nansum( ...
        (planes_loop(day).DIC_planes.DIC' - (planes_loop(day).DIC_planes.m_0 + ...
        planes_loop(day).DIC_planes.m_lon.*planes_loop(day).DIC_planes.plane_matrix(:,1) )).^2) ) ...
        ./ length(planes_loop(day).DIC_planes.DIC);    
    planes_loop(day).DIC_planes.m_lat_standarderror = sqrt( nansum( ...
        (planes_loop(day).DIC_planes.DIC' - (planes_loop(day).DIC_planes.m_0 + ...
        planes_loop(day).DIC_planes.m_lat.*planes_loop(day).DIC_planes.plane_matrix(:,2) )).^2) ) ...
        ./ length(planes_loop(day).DIC_planes.DIC) ;    
    planes_loop(day).DIC_planes.plane_standarderror = sqrt( nansum( ...
        (planes_loop(day).DIC_planes.DIC' - planes_loop(day).DIC_planes.DICfit(planes_loop(day).DIC_planes.plane_matrix(:,1:2))).^2 ) ...
        ./ length(planes_loop(day).DIC_planes.DIC) ) ; 
    
%% DIC depth resolved gridded  
    
% grid oxygen with similar bins to geopotential anomalies
xgrid = round(nanmin(planes_loop(day).DIC_planes.xa)):0.5:round(nanmax(planes_loop(day).DIC_planes.xa));
ygrid = round(nanmin(planes_loop(day).DIC_planes.ya)):0.5:round(nanmax(planes_loop(day).DIC_planes.ya));
[X,Y] = meshgrid(xgrid,ygrid);
zgrid = 1:5:60; % same as GPA

for n_bin = 1:numel(zgrid)
    for n_x = 1:numel(xgrid)
        for n_y = 1:numel(ygrid)
            check_bin = planes_loop(day).DIC_planes.P >= zgrid(n_bin) - 6 & planes_loop(day).DIC_planes.P < zgrid(n_bin) + 4 ...
                & planes_loop(day).DIC_planes.xa >= xgrid(n_x) & planes_loop(day).DIC_planes.xa < xgrid(n_x)+0.49 & ...
                planes_loop(day).DIC_planes.ya >= ygrid(n_y) & planes_loop(day).DIC_planes.ya < ygrid(n_y)+0.49;
            planes_loop(day).DIC_planes.gridded_DIC_var(n_bin,n_x,n_y) = nanmedian(planes_loop(day).DIC_planes.DIC(check_bin));
            planes_loop(day).DIC_planes.gridded_xa(n_bin,n_x,n_y) = xgrid(n_x);
            planes_loop(day).DIC_planes.gridded_ya(n_bin,n_x,n_y) = ygrid(n_y);    
            planes_loop(day).DIC_planes.gridded_z(n_bin) =  zgrid(n_bin);
        end
    end
end

% get fits
for n_bin = 1:numel(zgrid)
    xa = squeeze(planes_loop(day).DIC_planes.gridded_xa(n_bin,:,:)); 
    ya = squeeze(planes_loop(day).DIC_planes.gridded_ya(n_bin,:,:)); 
    DIC = squeeze(planes_loop(day).DIC_planes.gridded_DIC_var(n_bin,:,:));
    check_nans = isfinite(DIC);
    if numel(DIC(check_nans)) > 1
        % prepare for fit
        planes_loop(day).DIC_planes.gridded_const = ones(size(DIC(check_nans)));
        planes_loop(day).DIC_planes.gridded_plane_matrix = [xa(check_nans), ya(check_nans), planes_loop(day).DIC_planes.gridded_const];
        % fit function below gets same result as line above. 
        [planes_loop(day).DIC_planes.gridded_DICfit(n_bin).vals, planes_loop(day).DIC_planes.gridded_gof(n_bin).vals,~]  = ...
            fit(planes_loop(day).DIC_planes.gridded_plane_matrix(:,1:2),DIC(check_nans),'poly11');
        % get important variables   
        planes_loop(day).DIC_planes.gridded_m_lon(n_bin) = planes_loop(day).DIC_planes.gridded_DICfit(n_bin).vals.p10;
        planes_loop(day).DIC_planes.gridded_m_lat(n_bin)  = planes_loop(day).DIC_planes.gridded_DICfit(n_bin).vals.p01;
        planes_loop(day).DIC_planes.gridded_m_0(n_bin)  = planes_loop(day).DIC_planes.gridded_DICfit(n_bin).vals.p00;        
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
    
    %% DIC
    
    planes_loop(day).across_DIC.S = planes_loop(day).DIC_planes.DIC(...
        planes_loop(day).DIC_planes.lat > 43.2 & planes_loop(day).DIC_planes.lat  <=43.3);
    planes_loop(day).across_DIC.N = planes_loop(day).DIC_planes.DIC(...
        planes_loop(day).DIC_planes.lat  > 43.43 & planes_loop(day).DIC_planes.lat  < 43.5);
    planes_loop(day).across_DIC.W = planes_loop(day).DIC_planes.DIC(...
        planes_loop(day).DIC_planes.lon  > 7.64 & planes_loop(day).DIC_planes.lon  < 7.725);
    planes_loop(day).across_DIC.E = planes_loop(day).DIC_planes.DIC(...
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

    %% get DIC-specific time window means
    planes_loop(day).means.DIC_h = nanmean(planes_loop(day).DIC_planes.DIC);  % mmol m^-3    
    planes_loop(day).means.DIC_std_h = nanstd(planes_loop(day).DIC_planes.DIC);  % mmol m^-3 
    planes_loop(day).means.DIC_surf = nanmean(vars.DIC(planes_loop(day).day_selection_h & vars.depth < 10));  % mmol m^-3    
    planes_loop(day).means.DIC_surf_std = nanstd(planes_loop(day).day_selection_h & vars.depth < 10);  % mmol m^-3     
    planes_loop(day).means.DIC_inv_h = planes_loop(day).means.DIC_h*options.h;  % mmol m^-2
    planes_loop(day).means.DACu_h = nanmean(vars.DACs.DACu(planes_loop(day).day_selection_DAC));   
    planes_loop(day).means.DACv_h = nanmean(vars.DACs.DACv(planes_loop(day).day_selection_DAC));  
    planes_loop(day).means.DACu_std_h = nanstd(vars.DACs.DACu(planes_loop(day).day_selection_DAC));   
    planes_loop(day).means.DACv_std_h = nanstd(vars.DACs.DACv(planes_loop(day).day_selection_DAC));    
    
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
DIC_dayUmb = vars.DIC(tdayPUmb);    
planes_loop(day).means.DICinvt1MLDt2 = nanmean(DIC_dayUmb)*planes_loop(day+1).means.MLD_h;
planes_loop(day).means.DICinvt1MLDt2std = nanstd(DIC_dayUmb)*planes_loop(day+1).means.MLD_h;

end
% need to add NaN dummies so that horzcat works with structure 
for day = [options.dayrange(1) options.dayrange(end)]
    planes_loop(day).means.DICinvt1MLDt2 = NaN;
    planes_loop(day).means.DICinvt1MLDt2std = NaN;
end

clear DIC_dayUmb *dayPUmb ii bin_depth_check bn check* date_num day GPA GPA_* lat lon press_selection *_bin
%% save file for later use
save('planes','planes_loop');

%% if file already created
% else
%     disp(['DIC Plane-fits | Loading previously created file containing DIC Plane-fits']);
%     load('planes.mat');
% end
%     means_struct = [planes_loop.means];
%     disp(['DIC Plane-fits | Finished']);