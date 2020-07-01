%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_DIC Advection.m

% Script to calculate DIC Advection, following method by Alkire et al 2014

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get DIC gradients and U and V, calculate DIC Advection
% DIC_planes = [planes_loop.DIC_planes];
% for day = options.dayrange-8
% % DIC gradients
% DIC_adv(day).DIC_x = DIC_planes(day).m_lon./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
% DIC_adv(day).DIC_y = DIC_planes(day).m_lat./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
% % U and V
% DIC_adv(day).U = [DIC_adv(day).DU_abs_h]; % m s^-1
% DIC_adv(day).V = [DIC_adv(day).DV_abs_h]; % m s^-1 
% % error in gradients
% DIC_adv(day).DIC_x_error = DIC_planes(day).m_lon_standarderror./1000;
% DIC_adv(day).DIC_y_error = DIC_planes(day).m_lat_standarderror./1000;
% % calculate DIC Advection using DIC gradients and U and V velocities
% DIC_adv(day).adv = ((DIC_adv(day).DIC_x) .* DIC_adv(day).U) + ...
%     ((DIC_adv(day).DIC_y) .* DIC_adv(day).V);
% DIC_adv(day).adv = DIC_adv(day).adv * 86400 * options.h; 
% end

for day = options.dayrange
DIC_planes = planes_loop(day).DIC_planes;
% oxygen gradients
DIC_adv(day-8).oxy_x = DIC_planes.gridded_m_lon./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
DIC_adv(day-8).oxy_y = DIC_planes.gridded_m_lat./1000; % mmol m^-3 m^-1 d^-1, multiply by 1000 t get m instead of km
% U and V
DIC_adv(day-8).U = U(1:numel(DIC_adv(day-8).oxy_x)); % m s^-1
DIC_adv(day-8).V = V(1:numel(DIC_adv(day-8).oxy_x)); % m s^-1 
% error in gradients
% DIC_adv(day).oxy_x_error = DIC_planes(day).m_lon_standarderror./1000;
% DIC_adv(day).oxy_y_error = DIC_planes(day).m_lat_standarderror./1000;
% calculate O2 Advection using oxygen gradients and U and V velocities
DIC_adv(day-8).adv = ((DIC_adv(day-8).oxy_x) .* DIC_adv(day-8).U) + ...
    ((DIC_adv(day-8).oxy_y) .* DIC_adv(day-8).V); % mmol m^-3 s^-1 d^-1
DIC_adv(day-8).adv = DIC_adv(day-8).adv; % mmol m^-2 d^-1
% DIC_adv(day-8).adv = cumtrapz(DIC_planes.gridded_z(1:numel(DIC_adv(day-8).oxy_x)), DIC_adv(day-8).adv); % mmol m^-2 d^-1
% DIC_adv(day-8).adv = DIC_adv(day-8).adv(end); 
% DIC_adv(day).adv = DIC_adv(day).adv * 86400 * options.h;
DIC_adv(day-8).adv = nanmean(DIC_adv(day-8).adv) * 86400 * options.h; 
% DIC_adv(day-8).adv = cumtrapz(DIC_adv(day-8).adv);DIC_inv
% DIC_adv(day-8).adv = DIC_adv(day-8).adv(end)/options.h * 86400 * options.h;
end


disp('DIC Advection | Term calculated');

clear day press_selection