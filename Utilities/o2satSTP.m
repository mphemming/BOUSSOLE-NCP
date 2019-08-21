function c_sat = o2sat(temp, sal, pres)
% calculate O2 saturation concentration (in micromol/kg)
% as a function of temperature, salinity and pressure
%
% taken from Garcia & Gordon, Limnology & Oceanography 1992
% good for t_F <= t <= 40¼C and 0 < sal < 42 psu
% re-analysed data from Benson & Krause (1984)
% mixing ratio: 0.20946
%
% temp: temperature in degrees Celsius
% sal:  salinity in practical salinity units
% pres: pressure in hPa
% o2sat: O2 saturation concentration in µmol/kg

A = [3.80369 -0.0986643 5.10006 4.17887 3.20291 5.80871];
B = [-0.00951519 -0.0113864 -0.00770028 -0.00701577];
C = -2.75915e-07;

% Calculate scaled temperature following Hernan and Garcia (1992), p. 1310
Ts = log((298.15-temp) ./ (273.15+temp));

c_sat = (100*pres-vp_h2o(temp, sal)) ./ (101325-vp_h2o(temp, sal)).*exp(polyval(A,Ts) + sal.* polyval(B,Ts) + C.* sal.^2);
