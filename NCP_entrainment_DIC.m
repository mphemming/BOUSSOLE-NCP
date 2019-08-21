%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_DIC_Entrainment.m

% Script to calculate DIC_Entrainment 
% Based on equation used by Umberto in his thesis

% created by MPH in Norwich, 14/11/2017
% modified by MPH in Sydney, 01/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('DIC_Entrainment | Calculating effects of a changing MLD ');

%% calculate DIC_Entrainment for each day in a loop

for day = (options.dayrange(2):options.interval:options.dayrange(end)-options.interval)-8 % because of MLD diffs different range

    DIC_ent(day).MLDt1 = means_struct(day).MLD_h;
    DIC_ent(day).MLDt2 = means_struct(day+options.interval).MLD_h;
    DIC_ent(day).DICinvt1MLDt2 = means_struct(day).DICinvt1MLDt2;
    DIC_ent(day).DICinvht1 =  DIC_inv.inv_integral(day);

   if DIC_ent(day).MLDt2 > options.h
       if DIC_ent(day).MLDt1 < DIC_ent(day).MLDt2
          
           DIC_ent(day).dayindex(day) = day;

        DIC_ent(day).ent = ((DIC_ent(day).DICinvt1MLDt2 * ...
            (options.h / DIC_ent(day).MLDt2)) ...
            - DIC_ent(day).DICinvht1) / 1; % mmol m^-2   

       else
       DIC_ent(day).ent = 0;           
       end
   else
   DIC_ent(day).ent = 0;     
   end

end

disp('DIC_Entrainment | Effects calculated');