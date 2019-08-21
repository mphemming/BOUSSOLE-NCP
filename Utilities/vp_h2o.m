function vp = vp_h2o(temp, sal)
% Vapor pressure of seawater in Pa (Green & Carritt, 1967)
T = temp + 273.15;
vp = 101325 * (1-0.000537.*sal) .* exp(18.1973*(1-373.16./T) + 3.1813e-07*(1- exp(26.1205*(1-T/373.16))) - 0.018726*(1-exp(8.03945*(1-373.16./T))) + 5.02802*log(373.16./T));

