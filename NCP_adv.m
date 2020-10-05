%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_adv.m

% Script to calculate O2 Advection, following method by Alkire et al 2014

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 05/10/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('O2 & DIC Advection | Calculating absolute velocities');

%% calculate geostrophic absolute velocitie

for day = options.dayrange-8
    % Use only dives that go as far down as 1000m (ish)
     O2_adv(day).time_selection =  vars.DACs.t > planes_loop(day+8).date_num - (options.window*options.window_adv) ...
         &  vars.DACs.t < planes_loop(day+8).date_num + (options.window*options.window_adv);
    O2_adv(day).lat = vars.DACs.lat(O2_adv(day).time_selection);
    O2_adv(day).f = gsw_f(nanmedian(O2_adv(day).lat)); % coriolis parameter
    O2_adv(day).GV_D = [planes_loop(day+8).GPA_planes.D_levels]; % obtain bin depths within time window
    % concatenate GPA gradient profiles within time window
    O2_adv(day).D_U = [planes_loop(day+8).GPA_planes.gpa_xA]; % remember, 1day+-nday mean already done, so use just 'day'
    O2_adv(day).D_V = [planes_loop(day+8).GPA_planes.gpa_yA]; % remember, 1day+-nday mean already done, so use just 'day'
    O2_adv(day).D_U_errors = [planes_loop(day+8).GPA_planes.gpafit.error_x]; % remember, 1day+-nday mean already done, so use just 'day'
    O2_adv(day).D_V_errors = [planes_loop(day+8).GPA_planes.gpafit.error_y]; % remember, 1day+-nday mean already done, so use just 'day'    
    % Calculate velocity relative to the sea-surface:
    % 1000 is conversion between km and m, don't multiply by distance because gradients are per km
    O2_adv(day).DU = -O2_adv(day).D_V / (1000*O2_adv(day).f); % average over n day window, 
    O2_adv(day).DV = O2_adv(day).D_U / (1000*O2_adv(day).f); % average over n day window
    O2_adv(day).DU_errors = -O2_adv(day).D_V_errors / (1000*O2_adv(day).f); % average over n day window, 
    O2_adv(day).DV_errors = O2_adv(day).D_U_errors / (1000*O2_adv(day).f); % average over n day window    
    % Remove the mean from the rel. velocity to obtain mean anomalies    
    O2_adv(day).DU_anom = O2_adv(day).DU - mean(O2_adv(day).DU);
    O2_adv(day).DV_anom = O2_adv(day).DV - mean(O2_adv(day).DV);
    O2_adv(day).DU_anom_errors = O2_adv(day).DU_errors - mean(O2_adv(day).DU_errors);
    O2_adv(day).DV_anom_errors = O2_adv(day).DV_errors - mean(O2_adv(day).DV_errors);    
    % Add rel. vel. to time interval mean DAC to obtain approximation of sea-surface velocity:    
    O2_adv(day).DU_abs = O2_adv(day).DU_anom + means_struct(day).DACu_h; % using n-day averaged DAC
    O2_adv(day).DV_abs = O2_adv(day).DV_anom + means_struct(day).DACv_h; % using n-day averaged DAC
    % calc mean for h 
    O2_adv(day).DU_abs_h = nanmean(O2_adv(day).DU_abs(O2_adv(day).GV_D <= options.h));
    O2_adv(day).DV_abs_h = nanmean(O2_adv(day).DV_abs(O2_adv(day).GV_D <= options.h));
    % calc mean for surface 
    O2_adv(day).DU_abs_surf = nanmean(O2_adv(day).DU_abs(O2_adv(day).GV_D <= 10));
    O2_adv(day).DV_abs_surf = nanmean(O2_adv(day).DV_abs(O2_adv(day).GV_D <= 10));    
end

%% get oxygen gradients and U and V, calculate O2 Advection
% oxygen_planes = [planes_loop.oxygen_planes];
for day = options.dayrange
oxygen_planes = planes_loop(day).oxygen_planes;
U = [O2_adv(day-8).DU_abs];
V = [O2_adv(day-8).DV_abs];
% oxygen gradients
O2_adv(day-8).oxy_x = oxygen_planes.gridded_m_lon./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
O2_adv(day-8).oxy_y = oxygen_planes.gridded_m_lat./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
O2_adv(day-8).oxy_x_errors = oxygen_planes.gridded_m_lon_error./1000; % mmol m^-3 m^-1 d^-1, divide by 1000 t get m instead of km
O2_adv(day-8).oxy_y_errors = oxygen_planes.gridded_m_lat_error./1000; % mmol m^-3 m^-1 d^-1, divide by 1000 t get m instead of km
% U and V
O2_adv(day-8).U = U(1:numel(O2_adv(day-8).oxy_x)); % m s^-1
O2_adv(day-8).V = V(1:numel(O2_adv(day-8).oxy_x)); % m s^-1
% calculate O2 Advection using oxygen gradients and U and V velocities
O2_adv(day-8).adv = ((O2_adv(day-8).oxy_x) .* O2_adv(day-8).U) + ...
    ((O2_adv(day-8).oxy_y) .* O2_adv(day-8).V); % mmol m^-3 s^-1 
O2_adv(day-8).adv = nanmean(O2_adv(day-8).adv) * 86400 * options.h; % mmol m^-2 d^-1
O2_adv(day-8).adv_std = nanstd(((O2_adv(day-8).oxy_x) .* O2_adv(day-8).U) + ...
    ((O2_adv(day-8).oxy_y) .* O2_adv(day-8).V)) * 86400 * options.h;
end


%% get DIC gradients and U and V, calculate DIC Advection
% oxygen_planes = [planes_loop.oxygen_planes];
for day = options.dayrange
DIC_planes = planes_loop(day).DIC_planes;
U = [O2_adv(day-8).DU_abs];
V = [O2_adv(day-8).DV_abs];
% oxygen gradients
DIC_adv(day-8).DIC_x = DIC_planes.gridded_m_lon./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
DIC_adv(day-8).DIC_y = DIC_planes.gridded_m_lat./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
DIC_adv(day-8).DIC_x_errors = DIC_planes.gridded_m_lon_error./1000; % mmol m^-3 m^-1 d^-1, divide by 1000 t get m instead of km
DIC_adv(day-8).DIC_y_errors = DIC_planes.gridded_m_lat_error./1000; % mmol m^-3 m^-1 d^-1, divide by 1000 t get m instead of km
% U and V
DIC_adv(day-8).U = U(1:numel(DIC_adv(day-8).DIC_x)); % m s^-1
DIC_adv(day-8).V = V(1:numel(DIC_adv(day-8).DIC_x)); % m s^-1
% calculate O2 Advection using oxygen gradients and U and V velocities
DIC_adv(day-8).adv = ((DIC_adv(day-8).DIC_x) .* DIC_adv(day-8).U) + ...
    ((DIC_adv(day-8).DIC_y) .* DIC_adv(day-8).V); % mmol m^-3 s^-1 
DIC_adv(day-8).adv = nanmean(DIC_adv(day-8).adv) * 86400 * options.h; % mmol m^-2 d^-1
DIC_adv(day-8).adv_std = nanstd(((DIC_adv(day-8).DIC_x) .* DIC_adv(day-8).U) + ...
    ((DIC_adv(day-8).DIC_y) .* DIC_adv(day-8).V)) * 86400 * options.h;
end

disp('O2 & DIC Advection | Term calculated');

clear day press_selection