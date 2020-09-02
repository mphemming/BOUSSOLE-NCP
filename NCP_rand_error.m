%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to calculate the NCP erros using Pearson 
% random number generation

% created by MPH in Norwich, 04/07/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('NCP and errors (O2 & DIC) | Calculating ');

%% get errors for some variables

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
    
   %% DIC spatial standard deviation
    % mean for each profile within layer options.h
    
    errors.dives(day-8).divesDIC = vars.dive(errors.selection(day-8).tdayP); 
    errors.dives(day-8).divesDICunique = unique(errors.dives(day-8).divesDIC);
    errors.dives(day-8).divesDICsurf = vars.dive(errors.selection(day-8).tdayPsurf);     
    
    for divenumber = errors.dives(day-8).divesDICunique
        errors.spatial(day-8).DICsp = vars.DIC((errors.selection(day-8).tdayP));     
        errors.spatial(day-8).DICsplon = vars.lon((errors.selection(day-8).tdayP));
        errors.spatial(day-8).DICsplat = vars.lat((errors.selection(day-8).tdayP));
        errors.spatial(day-8).DICsptemp = vars.T((errors.selection(day-8).tdayPsurf));        
        errors.spatial(day-8).DICspsal = vars.S((errors.selection(day-8).tdayPsurf)); 
        errors.spatial(day-8).DICspt = datevec(vars.t((errors.selection(day-8).tdayP)));
        errors.spatial(day-8).DICspt = errors.spatial(day-8).DICspt(:,3);
        errors.spatial(day-8).DICsptsurf = datevec(vars.t((errors.selection(day-8).tdayPsurf)));
        errors.spatial(day-8).DICsptsurf = errors.spatial(day-8).DICsptsurf(:,3);       
        errors.spatial(day-8).DICspace(divenumber).mean = ...
            nanmean(errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == divenumber));  
        
        if day >= 9 & day <= 34  
        errors.spatial(day-8).DICspace(divenumber).inv_zmix2 = nanmean(...
            errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == divenumber)) .* ...
            means_struct(day-7).MLD_h; 
        errors.spatial(day-8).DICspace(divenumber).inv_zmix2_diff =  (nanmean(...
            vars.DIC(vars.dive == divenumber & vars.P > options.h & ...
            vars.P <  means_struct(day-7).MLD_h)));
        end
        
        errors.spatial(day-8).DICspace(divenumber).lonmean = ...
            nanmean(errors.spatial(day-8).DICsplon(errors.dives(day-8).divesDIC == divenumber)); 
        errors.spatial(day-8).DICspace(divenumber).latmean = ...
            nanmean(errors.spatial(day-8).DICsplat(errors.dives(day-8).divesDIC == divenumber)); 
        errors.spatial(day-8).DICspace(divenumber).Tmean = ...
            nanmean(errors.spatial(day-8).DICsptemp(errors.dives(day-8).divesDICsurf == divenumber));         
        errors.spatial(day-8).DICspace(divenumber).Smean = ...
            nanmean(errors.spatial(day-8).DICspsal(errors.dives(day-8).divesDICsurf == divenumber));  

        % only use values at edges of time window ( Toptions.his is wrong!?)
        
        errors.spatial(day-8).DICspace(divenumber).mean_end = ...
            nanmean(errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == ...
            divenumber & errors.spatial(day-8).DICspt' > nanmax(errors.spatial(day-8).DICspt)-1));        
        errors.spatial(day-8).DICspace(divenumber).mean_begin = ...
            nanmean(errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == ...
            divenumber & errors.spatial(day-8).DICspt' < nanmin(errors.spatial(day-8).DICspt)+1));            
        errors.spatial(day-8).DICspace(divenumber).std_end = ...
            nanstd(errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == ...
            divenumber & errors.spatial(day-8).DICspt' > nanmax(errors.spatial(day-8).DICspt)-1));        
        errors.spatial(day-8).DICspace(divenumber).std_begin = ...
            nanstd(errors.spatial(day-8).DICsp(errors.dives(day-8).divesDIC == ...
            divenumber & errors.spatial(day-8).DICspt' < nanmin(errors.spatial(day-8).O2spt)+1));  

    end
    
    %% inventory change errors (Oxygen)
    errors.invO2(day-8).stdO2space = nanstd([errors.spatial(day-8).O2space.mean]);
    errors.invO2(day-8).meanO2space = nanmean([errors.spatial(day-8).O2space.mean]);        
    errors.invO2(day-8).stdTspace = nanstd([errors.spatial(day-8).O2space.Tmean]);
    errors.invO2(day-8).meanTspace = nanmean([errors.spatial(day-8).O2space.Tmean]);       
    errors.invO2(day-8).stdSspace = nanstd([errors.spatial(day-8).O2space.Smean]);
    errors.invO2(day-8).meanSspace = nanmean([errors.spatial(day-8).O2space.Smean]); 
    % Only using values at beginning and end of time window
    errors.invO2(day-8).stdO2space_begin = nanstd([errors.spatial(day-8).O2space.mean_begin]);
    errors.invO2(day-8).stdO2space_end = nanstd([errors.spatial(day-8).O2space.mean_end]);       

    if day >=1 & day <=34;    
    errors.invO2(day-8).mean_inv_t1zmix2 = nanmean([errors.spatial(day-8).O2space.inv_zmix2]);     
    errors.invO2(day-8).std_inv_t1zmix2 = nanstd([errors.spatial(day-8).O2space.inv_zmix2]);    
    errors.invO2(day-8).mean_inv_t1zmix2_diff = nanmean([errors.spatial(day-8).O2space.inv_zmix2_diff]);     
    errors.invO2(day-8).std_inv_t1zmix2_diff = nanstd([errors.spatial(day-8).O2space.inv_zmix2_diff]);    
    end
    
    %% inventory change errors (DIC)
    errors.invDIC(day-8).stdDICspace = nanstd([errors.spatial(day-8).DICspace.mean]);
    errors.invDIC(day-8).meanDICspace = nanmean([errors.spatial(day-8).DICspace.mean]);        
    errors.invDIC(day-8).stdTspace = nanstd([errors.spatial(day-8).DICspace.Tmean]);
    errors.invDIC(day-8).meanTspace = nanmean([errors.spatial(day-8).DICspace.Tmean]);       
    errors.invDIC(day-8).stdSspace = nanstd([errors.spatial(day-8).DICspace.Smean]);
    errors.invDIC(day-8).meanSspace = nanmean([errors.spatial(day-8).DICspace.Smean]); 
    % Only using values at beginning and end of time window
    errors.invDIC(day-8).stdDICspace_begin = nanstd([errors.spatial(day-8).DICspace.mean_begin]);
    errors.invDIC(day-8).stdDICspace_end = nanstd([errors.spatial(day-8).DICspace.mean_end]);       

    if day >=1 & day <=34  
    errors.invDIC(day-8).mean_inv_t1zmix2 = nanmean([errors.spatial(day-8).DICspace.inv_zmix2]);     
    errors.invDIC(day-8).std_inv_t1zmix2 = nanstd([errors.spatial(day-8).DICspace.inv_zmix2]);    
    errors.invDIC(day-8).mean_inv_t1zmix2_diff = nanmean([errors.spatial(day-8).DICspace.inv_zmix2_diff]);     
    errors.invDIC(day-8).std_inv_t1zmix2_diff = nanstd([errors.spatial(day-8).DICspace.inv_zmix2_diff]);    
    end    
    
end

%% ASE errors

%--------------------------------
% Oxygen
%--------------------------------

errors.O2_ASE.calibration = 3.3; % mmol m-3
errors.O2_ASE.gas = 20; % 20%
simulated_errors.O2_ASE.calibration = pearsrnd(0,3.3,0,3,100000,1); % mmol m-3
simulated_errors.O2_ASE.gas = (pearsrnd(0,20,0,3,100000,1) / 100) + 1; % percent
simulated_errors.O2_ASE.O2_surf = [means_struct.O2_surf];
% simulated_errors.O2_ASE.calibration = random_range(0, 3.3,3000);% mmol m-3
% simulated_errors.O2_ASE.gas = random_range(0, 0.2,3000)+1;% percent
errors.O2_ASE.O2_sat = o2satSTP(O2_ase.Temp, O2_ase.Salt, O2_ase.press*1013.25) - ...
    o2satSTP(O2_ase.Temp+[errors.invO2.stdTspace], O2_ase.Salt+[errors.invO2.stdSspace], O2_ase.press*1013.25);
errors.O2_ASE.O2_sat = ([means_struct.sig0_surf]/1000) .* errors.O2_ASE.O2_sat; 
for ts = 1:numel(O2_ase.ASE)
    errors.O2_ASE.O2_sat_sim(ts,:) = pearsrnd(0,errors.O2_ASE.O2_sat(ts),0,3,100000,1); % mmol m-3
end
% ASE calculation
for ts = 1:numel(O2_ase.ASE)
    for n = 1:100000
        simulated_errors.O2_ASE.ASE_sim(ts,n) = ((O2_ase.KO2(ts))/100*24).*simulated_errors.O2_ASE.gas(n) ...
            .*((simulated_errors.O2_ASE.O2_surf(ts)+simulated_errors.O2_ASE.calibration(n) - ...
            (O2_ase.O2_saturation(ts)+simulated_errors.O2_ASE.calibration(n)+errors.O2_ASE.O2_sat_sim(ts,n) .* O2_ase.bub(ts)))) .*O2_ase.correction(ts) ;
    end
end

%% Advection errors

for day = options.dayrange-8
    
    % create DACs simulation
    simulated_errors.ADV(day).DACu = pearsrnd(0,means_struct(day).DACu_std_h,0,3,100000,1);
    simulated_errors.ADV(day).DACv = pearsrnd(0,means_struct(day).DACv_std_h,0,3,100000,1);   
    % create DU and DV anom, plus O2 plane errors
    for n_layer = 1:numel(O2_adv(day).oxy_x_errors)
        %
        DU_anom(day).vals(n_layer).vals = pearsrnd(0,abs(O2_adv(day).DU_anom_errors(n_layer))*0.341,0,3,100000,1);
        DV_anom(day).vals(n_layer).vals = pearsrnd(0,abs(O2_adv(day).DV_anom_errors(n_layer))*0.341,0,3,100000,1);
        %        
        simulated_errors.ADV(day).DU_abs_errors(n_layer,:) = DU_anom(day).vals(n_layer).vals + simulated_errors.ADV(day).DACu;
        simulated_errors.ADV(day).DV_abs_errors(n_layer,:) = squeeze([DV_anom(day).vals(n_layer).vals]) + simulated_errors.ADV(day).DACv;
        %
        simulated_errors.ADV(day).oxy_x(day,n_layer,:) = pearsrnd(0,abs(O2_adv(day).oxy_x_errors(n_layer))*0.341,0,3,100000,1);
        simulated_errors.ADV(day).oxy_y(day,n_layer,:) = pearsrnd(0,abs(O2_adv(day).oxy_y_errors(n_layer))*0.341,0,3,100000,1);        
    end
    
    % calculate ADV errors
    for n_layer = 1:numel(O2_adv(day).oxy_x_errors)    
        simulated_errors.ADV(day).ADV(n_layer,:) = ( squeeze(O2_adv(day).oxy_x(n_layer) + simulated_errors.ADV(day).oxy_x(day,n_layer,:)) .* ...
                                                    (O2_adv(day).U(n_layer) + simulated_errors.ADV(day).DU_abs_errors(n_layer,:))' ) + ...
                                                    ( squeeze(O2_adv(day).oxy_y(n_layer) + simulated_errors.ADV(day).oxy_y(day,n_layer,:)) .* ...
                                                        (O2_adv(day).V(n_layer) + simulated_errors.ADV(day).DV_abs_errors(n_layer,:))' ); % mmol m^-3 s^-1 
        simulated_errors.ADV(day).ADV_estimate = nanmean(simulated_errors.ADV(day).ADV,1)' * 86400 * options.h;
    end
    
%     simulated_errors.ADV(day).ADV_estimate_final = ...
%         nanmean(simulated_errors.ADV(day).ADV_estimate(simulated_errors.ADV(day).ADV_estimate ~= 0))* 86400 * options.h;
end

% old method
% for ts = 1:numel([O2_adv.adv])
%     errors.O2_ADV.ADV_std(ts,:) = pearsrnd(0,O2_adv(ts).adv_std,0,3,100000,1); % mmol m-3
%     simulated_errors.O2_ADV.ADV_std(ts,:) = errors.O2_ADV.ADV_std(ts,:) + O2_adv(ts).adv;
% end

%% Inventory change

% for day = 10:34
% errors.invO2(day-9).errors = sqrt( (errors.invO2(2).stdO2space_begin)^2 + (errors.invO2(25).stdO2space_end)^2)*options.h;
% errors.invDIC(day-9).errors = sqrt( (errors.invDIC(2).stdDICspace_begin)^2 + (errors.invDIC(25).stdDICspace_end)^2)*options.h;
% end
% 
% for ts = 1:numel([O2_adv.adv])
%     errors.O2_inv.inv(ts,:) = pearsrnd(0,errors.invO2(1).errors,0,3,100000,1); % mmol m-3
%     simulated_errors.O2_inv.inv(ts,:) = errors.O2_inv.inv(ts,:) + O2_inv.inv(ts);
% end

for day = options.dayrange-8
    errors.invO2(day).O2_std = means_struct(day).dives_O2_h_standard_error;
    simulated_errors.O2_inv.std(day,:) = pearsrnd(0,errors.invO2(day).O2_std,0,3,100000,1);
    % calculate errors
    for ts = 1:100000
        simulated_errors.O2_inv.inv_integral(ts,day) = simulated_errors.O2_inv.std(day,ts).*options.h;
    end
end
for ts = 1:100000
    simulated_errors.O2_inv.differentials(ts,:) = diff(simulated_errors.O2_inv.inv_integral(ts,:));
    simulated_errors.O2_inv.inv(ts,:) = interp1(O2_inv.diffrange, simulated_errors.O2_inv.differentials(ts,:), O2_inv.wantedrange,'linear','extrap');
end


%% entrainment

% errors.ENT.MLD_std = [means_struct.MLD_std_h]
errors.ENT.MLD_std = [means_struct.MLD_standard_error_h];

for ts = 2:27
    simulated_errors.ENT.MLD_error(ts,:) = pearsrnd(0,0.5,0,3,100000,1);
end

for ts = 2:26
    simulated_errors.ENT.MLDt1(ts,:) = means_struct(ts).MLD_h + simulated_errors.ENT.MLD_error(ts,:);
    simulated_errors.ENT.MLDt2(ts,:) = means_struct(ts+options.interval).MLD_h + simulated_errors.ENT.MLD_error(ts+1,:);
    simulated_errors.ENT.O2invt1MLDt2(ts,:) = simulated_errors.ENT.MLDt2(ts,:).*simulated_errors.O2_inv.inv(:,ts)'+([means_struct(ts-1).O2_h]*46);
    simulated_errors.ENT.O2invht1(ts,:) = simulated_errors.O2_inv.inv(:,ts);
    
    
  if O2_ent(ts).MLDt2 > options.h
       if O2_ent(ts).MLDt1 < O2_ent(ts).MLDt2
           
%        s = simulated_errors.ENT.MLDt2(ts,:);
%        s(s < 46) = NaN;
%         change = (options.h ./ s);
        change = (options.h / O2_ent(day).MLDt2);

        simulated_errors.ENT.ent(ts,:) = ((O2_ent(day).O2invt1MLDt2 .* change) - simulated_errors.ENT.O2invht1(ts,:)) / 1; % mmol m^-2   
        simulated_errors.ENT.ent(ts,:) = simulated_errors.ENT.ent(ts,:) - nanmean(simulated_errors.ENT.ent(ts,:));
        
       else
       simulated_errors.ENT.ent(ts,:) = zeros(1,100000);           
       end
   else
   simulated_errors.ENT.ent(ts,:) = zeros(1,100000);     
   end    
end

simulated_errors.ENT.ent(1,:) = zeros(1,100000); 
simulated_errors.ENT.ent(27,:) = zeros(1,100000); 

%% Estimate NCP error

ENT = [0 , [O2_ent.ent], 0];

for ts = 1:27
    % with inventory change
%     simulated_errors.O2_NCP(ts,:) = simulated_errors.O2_inv.inv(ts,:)  + simulated_errors.ADV(ts).ADV_estimate ...
%         + simulated_errors.O2_ASE.ASE_sim(ts,:) - ENT(ts);
    % without inventory change
    simulated_errors.O2_NCP(ts,:) = simulated_errors.O2_inv.inv(:,ts)' + simulated_errors.ADV(ts).ADV_estimate' ...
        + simulated_errors.O2_ASE.ASE_sim(ts,:) - ENT(ts);    
    % old method ADV: simulated_errors.O2_ADV.ADV_std(ts,:)
end








