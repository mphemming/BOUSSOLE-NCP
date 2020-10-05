% test script to see what effect kz might have on flux

%% get binned profiles per timestep

% O2
for day = options.dayrange
    planes_loop(day).day_selection= vars.t > (planes_loop(day).date_num-options.window) & ...
    vars.t < (planes_loop(day).date_num+options.window);  % only top h metres for most variables
    Unique_profs = unique(vars.profile_n(planes_loop(day).day_selection));
    for n = 1:numel(Unique_profs)
        kz(day-8).O2(n).bin = bin_variable_profile(vars.O2(planes_loop(day).day_selection & vars.profile_n == Unique_profs(n)),...
            vars.P(planes_loop(day).day_selection & vars.profile_n == Unique_profs(n)),0:2:300,0);         

        kz(day-8).O2grad(n) = interp1(kz(day-8).O2(n).bin.vertical, kz(day-8).O2(n).bin.mean_var,options.h+5,'Linear') ...
            - interp1(kz(day-8).O2(n).bin.vertical, kz(day-8).O2(n).bin.mean_var,options.h-5,'Linear');
        kz(day-8).time(n,:) =ones(size(kz(day-8).O2(n).bin.vertical)).*nanmedian(vars.t(planes_loop(day).day_selection & vars.profile_n == Unique_profs(n)));        
        kz(day-8).kz(n) = (0.0003*86400)* kz(day-8).O2grad(n) / 10;    
        kz(day-8).kz_alkire(n) = (0.0001*86400)* kz(day-8).O2grad(n) / 10;
        kz(day-8).kz_liliane_slides(n) = (0.001*86400)* kz(day-8).O2grad(n) / 10;   
        kz(day-8).kz_Huang(n) = (0.00001*86400)* kz(day-8).O2grad(n) / 10;   
        kz(day-8).kz_Jan(n) = (0.00003*86400)* kz(day-8).O2grad(n) / 10;  
    end
    kz(day-8).means.Kz = nanmean(kz(day-8).kz);
    kz(day-8).means.Kz_std = nanstd(kz(day-8).kz);
    kz(day-8).means.Kz_alkire = nanmean(kz(day-8).kz_alkire);
    kz(day-8).means.Kz_alkire_std = nanstd(kz(day-8).kz_alkire);    
    kz(day-8).means.Kz_alkire_std_error = kz(day-8).means.Kz_alkire_std/sqrt(numel(kz(day-8).kz_alkire(isfinite(kz(day-8).kz_alkire))));
    kz(day-8).means.Kz_liliane_slides = nanmean(kz(day-8).kz_liliane_slides);
    kz(day-8).means.Kz_liliane_slides_std = nanstd(kz(day-8).kz_liliane_slides);       
    kz(day-8).means.Kz_Huang = nanmean(kz(day-8).kz_Huang);
    kz(day-8).means.Kz_Huang_std = nanstd(kz(day-8).kz_Huang);          
    kz(day-8).means.Kz_Jan = nanmean(kz(day-8).kz_Jan);
    kz(day-8).means.Kz_Jan_std = nanstd(kz(day-8).kz_Jan); 
    kz(day-8).means.Kz_Jan_std_err = nanstd(kz(day-8).kz_Jan) / sqrt(numel(kz(day-8).kz_Jan)); 
end

a = [kz.means];
kz_ = [a.Kz_Jan];
kz_std = [a.Kz_Jan_std_err];
% references:
% [Jurado et al., 2012] for Alkire, North East Atlantic Ocean
% I think the 0.0003 kz comes from Copin-Montegut 1999



% DIC
for day = options.dayrange
    planes_loop(day).day_selection= vars.t > (planes_loop(day).date_num-options.window) & ...
    vars.t < (planes_loop(day).date_num+options.window);  % only top h metres for most variables
    kz(day-8).DIC = bin_variable_profile(vars.DIC(planes_loop(day).day_selection),...
        vars.P(planes_loop(day).day_selection),0:2:300,0); 
    kz(day-8).DICgrad = interp1(kz(day-8).DIC.vertical, kz(day-8).DIC.mean_var,options.h+5,'Linear') ...
        - interp1(kz(day-8).DIC.vertical, kz(day-8).DIC.mean_var,options.h-5,'Linear');
    kz(day-8).time =ones(size(kz(day-8).DIC.vertical)).*planes_loop(day).date_num;
    kz(day-8).kz_DIC = (3E-4*86400)* kz(day-8).DICgrad / 10;
    kz(day-8).kz_alkire_DIC = (0.0001*86400)* kz(day-8).DICgrad / 10;
end



% kz units: m2s-1 x mmol m-3 = mmol m-1 s-1 (need /10 for correct units) * m
% = mmol m-2 s-1

%% get Kz error estimate

for day = options.dayrange
        kz(day-8).stdev1_O2 = interp1(kz(day-8).O2.vertical, kz(day-8).O2.stdev_var,options.h+5,'Linear');
        kz(day-8).stdev2_O2 = interp1(kz(day-8).O2.vertical, kz(day-8).O2.stdev_var,options.h-5,'Linear');
        kz(day-8).stdev1_DIC = interp1(kz(day-8).DIC.vertical, kz(day-8).DIC.stdev_var,options.h+5,'Linear');
        kz(day-8).stdev2_DIC = interp1(kz(day-8).DIC.vertical, kz(day-8).DIC.stdev_var,options.h-5,'Linear');        
end
