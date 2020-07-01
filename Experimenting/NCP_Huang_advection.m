%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_Huang_advection.m

% Script to calculate Advection, following method by Alkire et al 2014

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Advection | Calculating absolute velocities');

%% calculate geostrophic absolute velocitie

for day = options.dayrange-8
    % Use only dives that go as far down as 1000m (ish)
     O2_adv(day).time_selection =  vars.DACs.t > planes_loop(day+8).date_num - options.window ...
         &  vars.DACs.t < planes_loop(day+8).date_num + options.window;
    O2_adv(day).lat = vars.DACs.lat(O2_adv(day).time_selection);
    O2_adv(day).f = gsw_f(nanmedian(O2_adv(day).lat)); % coriolis parameter
    O2_adv(day).GV_D = [planes_loop(day+8).GPA_planes.D_levels]; % obtain bin depths within time window
    % concatenate GPA gradient profiles within time window
    O2_adv(day).D_U = [planes_loop(day+8).GPA_planes.gpa_xA]; % remember, 1day+-2day mean already done, so use just 'day'
    O2_adv(day).D_V = [planes_loop(day+8).GPA_planes.gpa_yA]; % remember, 1day+-2day mean already done, so use just 'day'
    % Calculate velocity relative to the sea-surface:
    % 1000 is conversion between km and m, don't multiply by distance because gradients are per km
    O2_adv(day).DU = -O2_adv(day).D_V / (1000*O2_adv(day).f); % average over 3 day window, 
    O2_adv(day).DV = O2_adv(day).D_U / (1000*O2_adv(day).f); % average over 3 day window
    % Remove the mean from the rel. velocity to obtain mean anomalies    
    O2_adv(day).DU_anom = O2_adv(day).DU - mean(O2_adv(day).DU);
    O2_adv(day).DV_anom = O2_adv(day).DV - mean(O2_adv(day).DV);
    % Add rel. vel. to time interval mean DAC to obtain approximation of sea-surface velocity:    
    O2_adv(day).DU_abs = O2_adv(day).DU_anom + means_struct(day).DACu_h; % using 4-day averaged DAC
    O2_adv(day).DV_abs = O2_adv(day).DV_anom + means_struct(day).DACv_h; % using 4-day averaged DAC
    % calc mean for h 
    O2_adv(day).DU_abs_h = nanmean(O2_adv(day).DU_abs(O2_adv(day).GV_D <= options.h));
    O2_adv(day).DV_abs_h = nanmean(O2_adv(day).DV_abs(O2_adv(day).GV_D <= options.h));
    % calc mean for surface 
    O2_adv(day).DU_abs_surf = nanmean(O2_adv(day).DU_abs(O2_adv(day).GV_D <= 10));
    O2_adv(day).DV_abs_surf = nanmean(O2_adv(day).DV_abs(O2_adv(day).GV_D <= 10));    
end

%% get oxygen gradients and U and V, calculate advection
oxygen_planes = [planes_loop.oxygen_planes];
for day = options.dayrange-8
% oxygen gradients
O2_adv(day).oxy_x = oxygen_planes(day).m_lon./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
O2_adv(day).oxy_y = oxygen_planes(day).m_lat./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
% U and V
O2_adv(day).U = [O2_adv(day).DU_abs_h]; % m s^-1
O2_adv(day).V = [O2_adv(day).DV_abs_h]; % m s^-1 
% error in gradients
O2_adv(day).oxy_x_error = oxygen_planes(day).m_lon_standarderror./1000;
O2_adv(day).oxy_y_error = oxygen_planes(day).m_lat_standarderror./1000;
% calculate advection using oxygen gradients and U and V velocities
O2_adv(day).adv = ((O2_adv(day).oxy_x) .* O2_adv(day).U) + ...
    ((O2_adv(day).oxy_y) .* O2_adv(day).V);
O2_adv(day).adv = O2_adv(day).adv * 86400 * options.h; 
end

disp('Advection | Advection term calculated');

clear day press_selection