%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to calculate the air-sea exchange flux using a range of
% different parameterisations, uncertainty calculated using Monte-Carlo
% method
%

% created by MPH 23/03/2018 University of East Anglia, Norwich
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% temp = in situ temperature at surface (deg C)
% wind10m = wind speed at 10m above surface (m/s) 
% O2 = oxygen concentration at surface (mmol m^3); 
% conversion from µmol/kg to mmol m^3 = O2(µmol/kg) *(sigmatheta/1000) where 
% O2 and sigmatheta are vectors
% Orerror = standard deviation of mean surface oxygen values needed for
% uncertainty calculation (mmol m^3)
% O2sat = oxygen saturation (mmol m^3) (calculated using Garcia and Gordon (1992))
% Patm = atmospheric sea level pressure (atm)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% method = 
% 1. Alkire et al., (2014) paper using Wanninkhof et al., (2014) parameterizations
% 2. Wanninkhof et al., (1992)
% 3. Wanninkhof and Mcgills (1999)
% 4. Nightingale et al., (2000)
% 5. Stanley et al., (2009)
%
% bubble = yes/no 1/0 , if bubble injection wanted. Methods 1-4 use Wolf
% and Thorpe (1992), method 5 uses Stanley et al., (2009)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
% ASE = O2 air-sea exchange flux in mmol m^-2 d^-1
% Uncertainty = Calculated uncertainty using Monte-Carlo method


function [ASE uncertainty KO2 Sch bub] = ASEflux(temp,wind10m,meanwind10msquared,O2,O2error,O2sat,Patm,method,bubble)
%% bubble relationships using Fig. 5 from http://www.ldeo.columbia.edu/~spk/Papers/Nicholsonetal_InverseAirSeaFlux_11.pdf 

% bubblewind = [0 2.5 5 7.1 7.5 9 11.3 11.5]; % m/s
% Stanley = [0 1 5 13.1 22.5 54 85.9 105]*1000/365;% mmol mO2^-2 d-1
% Nicholson = [0 1 4 9.5 18 41.5 62.7 79]*1000/365;% mmol m^-2 d-1
% 
% [fitStanley,gofStanley,outStanley] = fit(bubblewind',Stanley','poly3');
% [fitNicholson,gofNicholson,outNicholson] = fit(bubblewind',Nicholson','poly3');

%% Alkire / Wanninkhof 2014

% The Schmidt number is the kinematic viscosity of water divided by the
% molecular diffusion coefficient of the gas, should = 568 at 20C in 35 salinity water
% The kinematic viscosity for fresh water and seawater are from Sharqawy et al. (2010). 
% The diffusion coefficients of gases are from the following: Ar, O2, N2, N2O, and CCl4 fit
% using Wilke and Chang (1955) as adapted by Hayduk and Laudie (1974).

% Schmidt number
ScO2Alkire =1920.4 - (135.6*temp) + (5.2122*temp.^2) -(0.10939*temp.^3) + (0.0009377*temp.^4); % Schmidt no.
% gas diffusivity
KO2Alkire = 0.251 .* meanwind10msquared .* ((ScO2Alkire/660).^-0.5); 
% bubble injection
if bubble == 1;
bub = (1+(0.01.*((sqrt(wind10m)/9).^2))); % bubble injection % Woolf and Thorpe
end
if bubble == 0;
    bub = 1;
end

ASEAlkire = (KO2Alkire/100*24).*(O2 - (O2sat .* bub));
bub_alkire = bub;
% uncertainty.ASEAlkireKO2= ((KO2Alkire/100*24).*(O2 - (O2sat .* bub))) - ((KO2Alkire.*1.2/100*24).*(O2 - (O2sat .* bub)));
uncertainty.ASEAlkireKO2 = 0.2*KO2Alkire;
uncertainty.ASEAlkireKO2val = KO2Alkire;
uncertainty.ASEAlkireSch = ScO2Alkire;

%% Wanninkhof 1992

ScO2Wann1992 = 1953.4 - (128*temp) + (3.9918*temp.^2) -(0.050091*temp.^3); % Schmidt no. table A1
KO2Wann1992 = 0.39 .* wind10m.^2 .* ((ScO2Wann1992/660).^-0.5); 

ScArWann1992 =1909.1 - (125.09*temp) + (3.9012*temp.^2) - (0.048953*temp.^3);

if bubble == 1;
bub = (1+(0.01.*((wind10m/9).^2))); % bubble injection % Woolf and Thorpe
end
if bubble == 0;
    bub = 1;
end

ASEWann1992 = (KO2Wann1992/100*24).*(O2 - (O2sat .* bub));


%% Nightingale 2000

KO2Night = 0.222  .* wind10m.^2 + 0.333  .* wind10m;

if bubble == 1;
bub = (1+(0.01.*((wind10m/9).^2))); % bubble injection % Woolf and Thorpe
end
if bubble == 0;
    bub = 1;
end

ASENight = (KO2Night/100*24).*(O2 - (O2sat .* bub));

%% Wanninkhof and McGills 1999 

KO2Wann1999 = 0.0283  .* wind10m.^3  .* ((ScO2Wann1992/660).^-0.5); 

if bubble == 1;
bub = (1+(0.01.*((wind10m/9).^2))); % bubble injection % Woolf and Thorpe
end
if bubble == 0;
    bub = 1;
end

ASEWann1999 = (KO2Wann1999/100*24).*(O2 - (O2sat .* bub));

%% Stanley 2009

%% no bubbles

Yg = 0.97; %+- 0.14

% Fge diffusive gas exchange flux
ASEStanley = (Yg * 8.6E-7) .* (ScO2Wann1992/660).^-0.5 .* wind10m.^2 .* (O2sat - O2); % mmol m^-2 s^-1
ASEStanley = ASEStanley*(60*60*24); % mmol m^-2 d^-1
ASEStanley = ASEStanley.*-1;

%% with bubbles

% trapped bubbles Fc

Ac = 9.1E-11; % tunable model parameter s^2m^-2
R = 8.31; % gas constant m^3 Pa mol^-1 K^-1
T = temp+273.15; % temperature Kelvin
PO2 = 0.21.*(Patm*101325); % Pa

Fc = Ac.*(wind10m - 2.27).^3 .* (PO2./(R.*T)); % mol m^2 s^s
Fc = Fc *1000; % convert to mmol m^2 s^-1
Fc = Fc*(60*60*24); % mmol m^-2 d^-1

% partially dissolved bubbles Fp

Ap = 2.3E-3; % S^-2 m^-2
Do = 1; % m^-2 S^-1
Di = 2.10E-5; % https://en.wikipedia.org/wiki/Mass_diffusivity   Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems (2nd ed.). New York: Cambridge University Press. ISBN 0-521-45078-0.
                     % at 25C temp
Di = Di./100^2; % convert to m^-2 s^-1 from cm^-2 s^-1, double-checked on internet and is correct

% for Pib
p = 1026; % density of seawater kg m^-3
g = 9.81; % gravity m^-1 s^-1
zbub = (0.15 .* wind10m - 0.55); % bubble sink depth in metres
Xi = 0.21; % mole fraction of O2 in dry air
% bunsen = 30; % Broecker, Wallace S., and Tsung-Hung Peng. "Tracers in the Sea." (1982): 125-159. https://www.ocean.washington.edu/courses/oc400/Lecture_Notes/CHPT11.pdf 
bunsen = 0.02764; % Weiss 1970, T = 14C, S = 38 ; https://sci-hub.la/https://www.sciencedirect.com/science/article/pii/0011747170900379


Pib = Xi.*((Patm*101325) + (p*g*zbub)); % pressure of gas in bubble (Pascal)

% partial pressure of O2
% DOmg = O2 / (44.66/1.42903);
% k = 1.25E-3*exp(-1700.*(1./(temp+273.15) - 1/298.15)); %  https://www.atmos-chem-phys.net/15/4399/2015/acp-15-4399-2015.pdf  or researchgate below, used latter for now
% % k = 1.2E-5; % (between  20C  and  25C  and  1 atm)
% Piw = DOmg./k; % https://www.researchgate.net/post/How_to_convert_partial_pressure_of_dissolved_oxygen_into_concentration_of_dissolved_oxygen                     

% PV = nRT, P = nRT/V
% http://web.mnstate.edu/jasperse/Chem210/Handouts/Ch%2010%20Gas%20Math%20Summary%20(p8).pdf 

Piw = (O2 .* 0.0821 .* 293.15) * 101325; % Pa

Fp = Ap * (wind10m - 2.27).^3 .* (bunsen*(Di/Do).^0.66) .* ((Xi.*(Patm*101325))./(R.*(temp+273.15))) .* (1 + ((p*g*zbub)./(Patm*101325)) - ((O2./O2sat)/1000));
Fp = Fp*1000; % mmol m^-2 s^-1
Fp = Fp*(60*60*24); % mmol m^-2 d^-1

ASEStanleybub = ASEStanley.*-1 + Fp + Fc;
ASEStanleybub = ASEStanleybub.*-1;

%% calculate uncertainty using monte-Carlo approach

% calculate errors in calibration/derivation of O2 and DIC
% using Monte-Carlo approach

% 
% bub1 = (1+(0.01.*((wind10m/9).^2))); % bubble injection % Woolf and Thorpe
% bub2 = 1 + (O2sat.*46).\fitStanley(wind10m)';
% bub3 = 1 + (O2sat.*46).\fitNicholson(wind10m)';    
% 
% KO2error = [nanmean([KO2Alkire,KO2Wann1992,KO2Night,KO2Wann1999]),nanstd([KO2Alkire,KO2Wann1992,KO2Night,KO2Wann1999])];
% bubbleerror = [nanmean([bub1,bub2,bub3]),nanstd([bub1,bub2,bub3])];
% % O2error is defined at start input
% O2saterror = [nanmean(O2sat),nanmean(0.02*(O2sat))]; % using 0.2% from Quay 2010 paper - for now
% 
% funToProp = @(x) x(1).*(x(3) - (x(4).*x(2)))
% % funToProp = @(x) x(2).^2 ./ x(1);
% 
% e1 = generateMCparameters('gaussian',[KO2error]); 
% e2 = generateMCparameters('gaussian',[bubbleerror]);
% e3 = generateMCparameters('gaussian',[nanmean(O2),nanmean(O2error)]);
% e4 = generateMCparameters('gaussian',[O2saterror]);
% paramMatrix = [e1;e2;e3;e4];
% 
% [MonteCarlo,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);


%% create figure





%% Output

% method = 1/2/3/4 = Alkire paper using Wanninkhof 2014/Wanninkhof 1992/Wanninkhof and
% Mcgills 1999/Nightingale 2000/ 

if method == 1;
    ASE = ASEAlkire;
    KO2 = KO2Alkire;
    Sch = ScO2Alkire;
    bub = bub_alkire;
end
if method == 2;
    ASE = ASEAWann1992;
end
if method == 3;
    ASE = ASEWann1999;
end
if method == 4;
    ASE = ASENight;
end
if method == 5 & bubble == 0;
    ASE = ASEStanley;
end

if method == 5 & bubble == 1;
    ASE = ASEStanleybub;
end

for step = 1:length(ASE);

uncertainty.mean(step) = nanmean([ASEAlkire(step), ASEWann1992(step),ASEWann1999(step),ASENight(step)]);
uncertainty.std(step) =  nanstd([ASEAlkire(step), ASEWann1992(step),ASEWann1999(step),ASENight(step)]);
uncertainty.bubble = bub;
end

end