%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to calculate the NCP erros using Pearson 
% random number generation

% created by MPH in Norwich, 04/07/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('NCP and errors (O2 & DIC) | Calculating ');

%% get errors for some variables


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

    if day >=1 & day <=34;    
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
simulated_errors.O2_ASE.calibration = pearsrnd(0,3.3*0.33,0,3,100000,1); % mmol m-3
simulated_errors.O2_ASE.gas = (pearsrnd(0,20*0.33,0,3,100000,1) / 100) + 1; % percent
simulated_errors.O2_ASE.O2_surf = [means_struct.O2_surf];
% simulated_errors.O2_ASE.calibration = random_range(0, 3.3,3000);% mmol m-3
% simulated_errors.O2_ASE.gas = random_range(0, 0.2,3000)+1;% percent
errors.O2_ASE.O2_sat = o2satSTP(O2_ase.Temp, O2_ase.Salt, O2_ase.press*1013.25) - ...
    o2satSTP(O2_ase.Temp+[errors.invO2.stdTspace], O2_ase.Salt+[errors.invO2.stdSspace], O2_ase.press*1013.25);
errors.O2_ASE.O2_sat = ([means_struct.sig0_surf]/1000) .* errors.O2_ASE.O2_sat; 
for ts = 1:numel(O2_ase.ASE)
    errors.O2_ASE.O2_sat_sim(ts,:) = pearsrnd(0,simulated_errors.O2_ASE.O2_sat(ts)*0.33,0,3,100000,1); % mmol m-3
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

for ts = 1:numel([O2_adv.adv])
    errors.O2_ADV.ADV_std(ts,:) = pearsrnd(0,O2_adv(ts).adv_std,0,3,100000,1); % mmol m-3
    simulated_errors.O2_ADV.ADV_std(ts,:) = errors.O2_ADV.ADV_std(ts,:) + O2_adv(ts).adv;
end

%% Inventory change

for day = 10:34
errors.invO2(day-9).errors = sqrt( (errors.invO2(2).stdO2space_begin)^2 + (errors.invO2(25).stdO2space_end)^2)*options.h;
errors.invDIC(day-9).errors = sqrt( (errors.invDIC(2).stdDICspace_begin)^2 + (errors.invDIC(25).stdDICspace_end)^2)*options.h;
end

for ts = 1:numel([O2_adv.adv])
    errors.O2_inv.inv(ts,:) = pearsrnd(0,errors.invO2(1).errors,0,3,100000,1); % mmol m-3
    simulated_errors.O2_inv.inv(ts,:) = errors.O2_inv.inv(ts,:) + O2_inv.inv(ts);
end


%% Estimate NCP error

ENT = [0 , [O2_ent.ent], 0];

for ts = 1:27
    % with inventory change
    simulated_errors.O2_NCP(ts,:) = simulated_errors.O2_inv.inv(ts,:)  + simulated_errors.O2_ADV.ADV_std(ts,:) ...
        + simulated_errors.O2_ASE.ASE_sim(ts,:) - ENT(ts);
    % without inventory change
%     simulated_errors.O2_NCP(ts,:) = O2_inv.inv(ts)  + simulated_errors.O2_ADV.ADV_std(ts,:) ...
%         + simulated_errors.O2_ASE.ASE_sim(ts,:) - ENT(ts);    
end








