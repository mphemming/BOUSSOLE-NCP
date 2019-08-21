

function [F_CO2, dpCO2,Sc,K,a]=FCO2(pCO2_agua, pCO2_atm,T,S,u,umeansquared)
   
    %%Function for the calculation or air-sea CO2 flux
        %   Cecilia Chapa Balcorta 
        %   Ensenada, Baja California
        %   UNIVERSIDAD AUTÓNOMA DE BAJA CALIFORNIA/ UNIVERSIDAD DEL MAR
    
        %   filename: FCO2.m
        %   created: 07/04/2014
        %   modified: 03/23/2015
        %   copyright @ 2015 Cecilia Chapa-Balcorta
        
       
        %INPUT
    
   %pCO2_agua= seawater pCO2 (uatm)
   %pCO2_atm=  atmospheric pCO2 (uatm)
   %T=  Temperature (Celsius)
   %S=  Salinity 
   %u = Wind speed (m/s) %%% THIS IS INCORRECT, SHOULD be mean(U^2) not
   %(mean(U))^2, INPUT MUST BE WIND SQUARED, NOT NORMAL WIND
   
    
   %Air-sea CO2 is calculated as follows:
   
           % FCO2 =K*a(dpCO2) 
    
    %Where
    %K=is the transfer velocity according to Wanninkhof (2014).
    %a = CO2 solibility constant according to Weiss (1974)
    %dpCO2 is the difference of air and seawater pCO2 
    
    
    %%%%% CO2 Transfer velocity calculation %%%%%%%%%
       
   Sc=Schmidt(T);

%    K=0.251*(u.^2).*((Sc./660).^-0.5); % Wanninkhof 2014 %% THIS BIT WAS
%    INCORRECT

   K=0.251*umeansquared.*((Sc./660).^-0.5); % Wanninkhof 2014 
   
   dpCO2=pCO2_agua-pCO2_atm;     %%calculation of delta pCO2
   a=Ko_weiss(T,S);    %Solibility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1
    
F_CO2 =0.24*K.*a.*dpCO2; %CO2 flux (mmol m^-2 d^-1)
end

%%%%%Subrutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%****Solibuility constant (Weiss, 1974) ***********************************

function [Ko]=Ko_weiss(T,S)
%
A=[-60.2409, 93.4517, 23.3585];  %mol/Kg.atm
B=[0.023517, -0.023656, 0.0047036]; %mol/Kg.atm
T=T+273.15; %Conversion from Celsius degrees to Kelvins
Ln_Ko=A(1)+(A(2).*(100./T))+(A(3).*log(T./100))+S.*(B(1)+(B(2).*(T./100))+(B(3).*(T./100).^2));
Ko=exp(Ln_Ko);

end

%******** Schmidt Number*********

    %For water of salinity=35 and temperature range 0-30°C    %%%%%%%%%%%%%
    
    function [Sc]=Schmidt(T)
            A = 2116.8;     B = 136.25;     C = 4.7353;     D = 0.092307; E = 0.0007555;
            Sc= A - (B.*T)+(C.*T.^2)-(D.*T.^3)+(E.*T.^4);
            
    end
        
