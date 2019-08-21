% test script to see what effect kz might have on flux

%% get binned profiles per timestep
for day = options.dayrange
    planes_loop(day).day_selection= vars.t > (planes_loop(day).date_num-options.window) & ...
    vars.t < (planes_loop(day).date_num+options.window);  % only top h metres for most variables
    kz(day-8).O2 = bin_variable_profile(vars.O2(planes_loop(day).day_selection),...
        vars.P(planes_loop(day).day_selection),0:2:300,0); 
    kz(day-8).O2grad = interp1(kz(day-8).O2.vertical, kz(day-8).O2.mean_var,options.h+5,'Linear') ...
        - interp1(kz(day-8).O2.vertical, kz(day-8).O2.mean_var,options.h-5,'Linear');
    kz(day-8).time =ones(size(kz(day-8).O2.vertical)).*planes_loop(day).date_num;
    kz(day-8).kz = (3E-4*86400)* kz(day-8).O2grad / 10;
    kz(day-8).kz_alkire = (0.0001*86400)* kz(day-8).O2grad / 10;
end

% kz units: m2s-1 x mmol m-3 = mmol m-1 s-1 (need /10 for correct units) * m
% = mmol m2 s-1

%% get Kz error estimate

for day = options.dayrange
        kz(day-8).stdev1 = interp1(kz(day-8).O2.vertical, kz(day-8).O2.stdev_var,options.h+5,'Linear');
        kz(day-8).stdev2 = interp1(kz(day-8).O2.vertical, kz(day-8).O2.stdev_var,options.h-5,'Linear');
end
