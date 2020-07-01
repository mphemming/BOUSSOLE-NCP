%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to calculate toptions.he propogated errors in the
% NCP calculation

% created by MPH in Norwich, 08.03.2018
% modified by MPH in Sydney, 03/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('NCP and errors | Calculating ');

%% resources

% options.http://lectureonline.cl.msu.edu/~mmp/lvars.DACss/error/e2.options.htm 
% options.https://poptions.hysics.appstate.edu/undergraduate-programs/lvars.DACsoratory/resources/error-propagation

%% obtain standard deviations in space for time window
    
for day = options.dayrange

    %% data selection logical vector to include data 2 days eitoptions.her side of day
    % different vectors because of different dimensions
    errors.selection(day-8).tday = vars.t > planes_loop(day).date_num - options.window ...
        & vars.t < planes_loop(day).date_num + options.window;    % full profile deptoptions.h
    errors.selection(day-8).tdayP = vars.t > planes_loop(day).date_num - options.window ...
        & vars.t < planes_loop(day).date_num + options.window & vars.P < options.h;  % only top options.h metres
    errors.selection(day-8).tdayPsurf = vars.t > planes_loop(day).date_num - options.window ... 
        & vars.t < planes_loop(day).date_num + options.window & vars.P < 10;  % only tsurface 
    errors.selection(day-8).tdayDAC = vars.DACs.t > planes_loop(day).date_num - options.window ...
        & vars.DACs.t < planes_loop(day).date_num + options.window; % for DACs

    %% DAC errors

    errors.ADV(day-8).meanlon = nanmean(vars.DACs.lon(errors.selection(day-8).tdayDAC));
    errors.ADV(day-8).meanlat = nanmean(vars.DACs.lat(errors.selection(day-8).tdayDAC));    
    errors.ADV(day-8).meanDACu = nanmean(vars.DACs.DACu(errors.selection(day-8).tdayDAC));
    errors.ADV(day-8).meanDACv = nanmean(vars.DACs.DACv(errors.selection(day-8).tdayDAC));
    errors.ADV(day-8).stdDACu = nanstd(vars.DACs.DACu(errors.selection(day-8).tdayDAC));    
    errors.ADV(day-8).stdDACv = nanstd(vars.DACs.DACv(errors.selection(day-8).tdayDAC));     
    
    %% oxygen spatial standard deviation
    % mean for each profile within layer options.h
    
    errors.dives(day-8).divesO2 = vars.dive(errors.selection(day-8).tdayP); 
    errors.dives(day-8).divesO2unique = unique(errors.dives(day-8).divesO2);
    errors.dives(day-8).divesO2surf = vars.dive(errors.selection(day-8).tdayPsurf);     
    
    for divenumber = errors.dives(day-8).divesO2unique
        errors.spatial(day-8).O2sp = vars.O2((errors.selection(day-8).tdayP));     
        errors.spatial(day-8).O2splon = vars.lon((errors.selection(day-8).tdayP));
        errors.spatial(day-8).O2splat = vars.lat((errors.selection(day-8).tdayP));
        errors.spatial(day-8).O2sptemp = vars.T((errors.selection(day-8).tdayPsurf));        
        errors.spatial(day-8).O2spsal = vars.S((errors.selection(day-8).tdayPsurf)); 
        errors.spatial(day-8).O2spt = datevec(vars.t((errors.selection(day-8).tdayP)));
        errors.spatial(day-8).O2spt = errors.spatial(day-8).O2spt(:,3);
        errors.spatial(day-8).O2sptsurf = datevec(vars.t((errors.selection(day-8).tdayPsurf)));
        errors.spatial(day-8).O2sptsurf = errors.spatial(day-8).O2sptsurf(:,3);       
        errors.spatial(day-8).O2space(divenumber).mean = ...
            nanmean(errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == divenumber));  
        
        if day >= 9 & day <= 34  
        errors.spatial(day-8).O2space(divenumber).inv_zmix2 = nanmean(...
            errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == divenumber)) .* ...
            means_struct(day-7).MLD_h; 
        errors.spatial(day-8).O2space(divenumber).inv_zmix2_diff =  (nanmean(...
            vars.O2(vars.dive == divenumber & vars.P > options.h & ...
            vars.P <  means_struct(day-7).MLD_h)));
        end
        
        errors.spatial(day-8).O2space(divenumber).lonmean = ...
            nanmean(errors.spatial(day-8).O2splon(errors.dives(day-8).divesO2 == divenumber)); 
        errors.spatial(day-8).O2space(divenumber).latmean = ...
            nanmean(errors.spatial(day-8).O2splat(errors.dives(day-8).divesO2 == divenumber)); 
        errors.spatial(day-8).O2space(divenumber).Tmean = ...
            nanmean(errors.spatial(day-8).O2sptemp(errors.dives(day-8).divesO2surf == divenumber));         
        errors.spatial(day-8).O2space(divenumber).Smean = ...
            nanmean(errors.spatial(day-8).O2spsal(errors.dives(day-8).divesO2surf == divenumber));  

        % only use values at edges of time window ( Toptions.his is wrong!?)
        
        errors.spatial(day-8).O2space(divenumber).mean_end = ...
            nanmean(errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == ...
            divenumber & errors.spatial(day-8).O2spt' > nanmax(errors.spatial(day-8).O2spt)-1));        
        errors.spatial(day-8).O2space(divenumber).mean_begin = ...
            nanmean(errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == ...
            divenumber & errors.spatial(day-8).O2spt' < nanmin(errors.spatial(day-8).O2spt)+1));            
        errors.spatial(day-8).O2space(divenumber).std_end = ...
            nanstd(errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == ...
            divenumber & errors.spatial(day-8).O2spt' > nanmax(errors.spatial(day-8).O2spt)-1));        
        errors.spatial(day-8).O2space(divenumber).std_begin = ...
            nanstd(errors.spatial(day-8).O2sp(errors.dives(day-8).divesO2 == ...
            divenumber & errors.spatial(day-8).O2spt' < nanmin(errors.spatial(day-8).O2spt)+1));  

    end
    %% inventory change errors 
    errors.inv(day-8).stdO2space = nanstd([errors.spatial(day-8).O2space.mean]);
    errors.inv(day-8).meanO2space = nanmean([errors.spatial(day-8).O2space.mean]);        
    errors.inv(day-8).stdTspace = nanstd([errors.spatial(day-8).O2space.Tmean]);
    errors.inv(day-8).meanTspace = nanmean([errors.spatial(day-8).O2space.Tmean]);       
    errors.inv(day-8).stdSspace = nanstd([errors.spatial(day-8).O2space.Smean]);
    errors.inv(day-8).meanSspace = nanmean([errors.spatial(day-8).O2space.Smean]); 
    % Only using values at beginning and end of time window
    errors.inv(day-8).stdO2space_begin = nanstd([errors.spatial(day-8).O2space.mean_begin]);
    errors.inv(day-8).stdO2space_end = nanstd([errors.spatial(day-8).O2space.mean_end]);       

    if day >=1 & day <=34;    
    errors.inv(day-8).mean_inv_t1zmix2 = nanmean([errors.spatial(day-8).O2space.inv_zmix2]);     
    errors.inv(day-8).std_inv_t1zmix2 = nanstd([errors.spatial(day-8).O2space.inv_zmix2]);    
    errors.inv(day-8).mean_inv_t1zmix2_diff = nanmean([errors.spatial(day-8).O2space.inv_zmix2_diff]);     
    errors.inv(day-8).std_inv_t1zmix2_diff = nanstd([errors.spatial(day-8).O2space.inv_zmix2_diff]);    
    end
    
%     clear O2space
end



%% Inventory change

% error of selecting Z_eu deptoptions.h

errors.O2cal.medianO2 = 182; % µmol kg
errors.O2cal.bottlermse_march = 2.39; % µmol kg
errors.O2cal.bottlermse_April = 0.85; % µmol kg
errors.O2cal.Gliderrmse_march = 2.60; % µmol kg
errors.O2cal.Gliderrmse_April = 4.32; % µmol kg
errors.O2Zeu.std_dev = O2_inv.inv_std;

for day = 10:34
errors.inv(day-9).errors = sqrt( (errors.inv(2).stdO2space_begin)^2 + (errors.inv(25).stdO2space_end)^2)*options.h;
end

%% Advection
for day = 1:numel(errors.ADV)
    errors.ADV(day).U_error = errors.ADV(day).stdDACu;
    errors.ADV(day).V_error = errors.ADV(day).stdDACv;
    errors.ADV(day).E_adv = sqrt( (O2_adv(day).oxy_x .* errors.ADV(day).U_error).^2 + ...
        (O2_adv(day).U .* O2_adv(day).oxy_x_error).^2 + ...
        (O2_adv(day).V .* O2_adv(day).oxy_y_error).^2 + ...
        (O2_adv(day).oxy_y .* errors.ADV(day).V_error).^2  );
    errors.ADV(day).errors_adv = errors.ADV(day).E_adv .* 86400 .* options.h;
end

%% Air-sea exchange
% load bottle and CTD-glider O2 relationship data
load('Bottle')
load('CTDgli','CTDgli')

errors.ASE.KO2 = O2_ase.ASE_uncertainty.ASEAlkireKO2/100*24; % determined from 20% uncertainty of KO2 and using ASEflux function, mmol m^-2 d^-1
errors.ASE.o2sat = o2satSTP(O2_ase.Temp, O2_ase.Salt, O2_ase.press*1013.25) - ...
    o2satSTP(O2_ase.Temp+[errors.inv.stdTspace], O2_ase.Salt+[errors.inv.stdSspace], O2_ase.press*1013.25)
errors.ASE.o2sat = ([means_struct.sig0_surf]/1000) .* errors.ASE.o2sat; 
errors.ASE.bottlermse_march = 2.39; % µmol kg
errors.ASE.bottlermse_April = 0.85; % µmol kg
errors.ASE.Gliderrmse_march = 2.60; % µmol kg
errors.ASE.Gliderrmse_April = 4.32; % µmol kg
errors.ASE.bottle_march = abs(Bottle.O2_Mar_fit.p1 - 0.7858);
errors.ASE.bottle_april = abs(Bottle.O2_Apr_fit.p1 - 0.7559);
errors.ASE.Mar_bottSE = Bottle.O2_Mar_standarderror;
errors.ASE.Apr_bottSE = Bottle.O2_Apr_standarderror;
errors.ASE.Mar_GliSE = CTDgli.Mar_standarderror;
errors.ASE.Apr_GliSE = CTDgli.Apr_standarderror;
errors.ASE.KO2vals = O2_ase.ASE_uncertainty.ASEAlkireKO2val/100*24;
errors.ASE.Schmidt = O2_ase.ASE_uncertainty.ASEAlkireSch;

errors.ASE.error_c = sqrt( (errors.ASE.Mar_bottSE)^2 + (errors.ASE.Mar_GliSE)^2);
errors.ASE.errors_ASE = sqrt(  (([means_struct.O2_surf] - O2_ase.O2_saturation) ...
    .* errors.ASE.KO2).^2 + (errors.ASE.KO2vals .* errors.ASE.error_c).^2 + ...
    (errors.ASE.KO2vals .* errors.ASE.o2sat).^2);

%% Entrainment

errors.ent.Inv_zlim = options.h .* [errors.inv.stdO2space];
errors.ent.Inv_zmix2 =[errors.inv.std_inv_t1zmix2];

errors.ent.errors_ent = options.h * [errors.inv(2:end-1).std_inv_t1zmix2_diff];
errors.ent.errors_ent([O2_ent.ent] == 0) = 0;

%% Overall error

E_inv = [errors.inv.errors]; E_ADV = [errors.ADV.errors_adv]; 
E_ase = [errors.ASE.errors_ASE]; E_ent = errors.ent.errors_ent;

errors.error_NCP =  sqrt( (E_inv).^2 + (E_ase(2:end-1)).^2 + ...
    (E_ADV(2:end-1)).^2 +  (E_ent).^2 );

clear E_inv E_ADV E_ase Bottle CTDgli day divenumber ADV E_ent 

disp('NCP and errors | finished ');
disp('Script | finished :-) ');