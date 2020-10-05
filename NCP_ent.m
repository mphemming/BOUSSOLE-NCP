%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_ent.m

% Script to calculate entrainment 
% Based on equation used by Umberto in his thesis

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 05/102019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Entrainment | Calculating effects of a changing MLD ');

%% calculate entrainment for each day in a loop

for day = (options.dayrange(2):options.interval:options.dayrange(end)-options.interval)-8 % because of MLD diffs different range

    O2_ent(day).MLDt1 = means_struct(day).MLD_h_profs;
    O2_ent(day).MLDt2 = means_struct(day+options.interval).MLD_h_profs;
    O2_ent(day).O2invt1MLDt2 = means_struct(day).O2invt1MLDt2_profs;
    O2_ent(day).O2invht1 =  O2_inv.inv_integral(day);
    DIC_ent(day).DICinvt1MLDt2 = means_struct(day).DICinvt1MLDt2_profs;
    DIC_ent(day).DICinvht1 =  DIC_inv.inv_integral(day);    

    %% oxygen
   if O2_ent(day).MLDt2 > options.h
       if O2_ent(day).MLDt1 < O2_ent(day).MLDt2
          
           O2_ent(day).dayindex(day) = day;

        O2_ent(day).ent = ((O2_ent(day).O2invt1MLDt2 * ...
            (options.h / O2_ent(day).MLDt2)) ...
            - O2_ent(day).O2invht1) / 1; % mmol m^-2   
        
%         O2_ent(day).ent_2 = options.h * ( (O2_ent(day).O2invt1MLDt2 - O2_ent(day).O2invht1) )/ 1;
        % I don't know what Jan was talking about - look up entrainment
        % flux in Emmerson paper!
       else
       O2_ent(day).ent = 0;           
       end
   else
   O2_ent(day).ent = 0;     
   end
   %% DIC
     if O2_ent(day).MLDt2 > options.h
       if O2_ent(day).MLDt1 < O2_ent(day).MLDt2
          
           DIC_ent(day).dayindex(day) = day;

        DIC_ent(day).ent = ((DIC_ent(day).DICinvt1MLDt2 * ...
            (options.h / O2_ent(day).MLDt2)) ...
            - DIC_ent(day).DICinvht1) / 1; % mmol m^-2   
        
%         DIC_ent(day).ent_2 = options.h * ( (DIC_ent(day).DICinvt1MLDt2 - DIC_ent(day).DICinvht1) )/ 1;
        % I don't know what Jan was talking about - look up entrainment
        % flux in Emmerson paper!
       else
       DIC_ent(day).ent = 0;           
       end
   else
   DIC_ent(day).ent = 0;     
   end

end

disp('Entrainment | Effects calculated');