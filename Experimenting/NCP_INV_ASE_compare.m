%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCP_INV_ASE_compare.m

% created by MPH in Sydney, 18/08/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Inventory change

% buoy
load([options.data_dir,'BOUSSOLE.mat']);  
BUOY.O2 = BOUSSOLE.O2_raw_10m_Liliane;
check = isfinite(BOUSSOLE.time_CSYS) & isfinite(BOUSSOLE.rho_CSYS);
BUOY.dens = interp1(BOUSSOLE.time_CSYS(check),BOUSSOLE.rho_CSYS(check),BOUSSOLE.O2_raw_10m_Liliane_date,'Linear');
BUOY.O2 = (BUOY.dens/1000) .* BUOY.O2;
BUOY.O2_date = BOUSSOLE.O2_raw_10m_Liliane_date;
% get 4-day means
for n = 1:27
    check_t = BUOY.O2_date >= datenum(2016,03,options.dayrange(n)) - 2 & BUOY.O2_date < datenum(2016,03,options.dayrange(n)) + 2;
    BUOY.O2_bin(n) = nanmean(BUOY.O2(check_t));
end
% diff
BUOY.O2_diff = interp1(1.5:26.5,diff(BUOY.O2_bin),1:27);

figure('units','normalized','position',[.1 .1 .7 .6])
p1 = plot(datenum(2016,03,options.dayrange),[O2_inv.inv],'LineWidth',2);
hold on;
p2 = plot(datenum(2016,03,options.dayrange),BUOY.O2_diff*46,'LineWidth',2);
datetick('x','dd/mm');
set(gca,'FontSize',16,'LineWidth',2);
ylabel('Inventory Change [mmol m^{-3}]');
add_zero
set(gcf,'Color','W');
leg = legend([p1 p2],'Glider','Buoy');
set(leg,'Location','NorthWest','Box','Off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: complete plot, and also look at distance between buoy and glider. can this explain differences?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASE

%buoy
load([options.directory,'/data/3hrbins'],'bins')
t = [bins.Bt];
O2 = [bins.O2];
dens =[bins.dens];
T = [bins.BT];
S = [bins.BS];
wind10 = [bins.wind10];
atmpress = [bins.press];
% get 4-day binned data
for n = 1:27
    check_t = t >= datenum(2016,03,options.dayrange(n)) - 2 & t < datenum(2016,03,options.dayrange(n)) + 2;
    O2_d(n) = nanmean(O2(check_t));
    dens_d(n) = nanmean(dens(check_t));    
    T_d(n) = nanmean(T(check_t));    
    S_d(n) = nanmean(S(check_t));    
    wind10_d(n) = nanmean(wind10(check_t));    
    atmpress_d(n) = nanmean(atmpress(check_t));    
end

O2_saturation = o2satSTP(T_d,S_d, atmpress_d/100);
% convert to mmol m^-3
O2_saturation = (dens_d/1000) .* O2_saturation;
% schmidt
ScO2 = 1920.4 - (135.6 * T_d) + (5.2122 * T_d.^2) ...
    - (0.10939 * T_d.^3) + (0.0009377 * T_d.^4); 
% Gas diffusivity
KO2 = 0.251 .* wind10_d.^2 .* ((ScO2/660).^-0.5); % gas diffusivity
% Calculate airsea exchange
% gas diffusivity, bubble, and schmidt also calculated in funciton 'ASEflux'
[ASE ASE_uncertainty] = ASEflux(T_d, wind10_d, wind10_d.^2, ...
    O2_d,ones(size(O2_d)),O2_saturation,atmpress_d,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: input data correct? would need to compare over same 4-day
% averaging, why differences in ASE at some times???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('units','normalized','position',[.1 .1 .7 .6])
p1 = plot(datenum(2016,03,options.dayrange),[O2_ase.ASE],'LineWidth',2);
hold on;
p2 = plot(datenum(2016,03,options.dayrange),ASE,'LineWidth',2);
datetick('x','dd/mm');
set(gca,'FontSize',16,'LineWidth',2);
ylabel('ASE [mmol m^{-3}]');
add_zero
set(gcf,'Color','W');
leg = legend([p1 p2],'Glider','Buoy');
set(leg,'Location','NorthWest','Box','Off');
ylim([-120 60])




