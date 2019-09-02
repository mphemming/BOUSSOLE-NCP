%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N_table.m

% Script to get NCP mean estimates and errors for paper table

% created by MPH in Norwich, 26/08/2019
% updated by MPH in Sydney, 31/08/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% glider NCP

%% O2

N_table.O2_glider.ADV = [O2_adv(2:end-1).adv];
N_table.O2_glider.ASE = [O2_ase.ASE(2:end-1)];
N_table.O2_glider.ADV_error = [errors.ADV(2:end-1).errors_adv];
N_table.O2_glider.ASE_error = [errors.ASE.errors_ASE(2:end-1)];
N_table.O2_glider.NCP = NCP_est;
N_table.O2_glider.NCP_no_adv = NCP_est_no_adv;
N_table.O2_glider.NCP_error = [errors.error_NCP];
N_table.O2_glider.NCP_no_adv_error = [errors.error_NCP_no_adv];
N_table.O2_glider.days = datenum(2016,03,10):1:datenum(2016,04,03);

% 10-19 March
check = N_table.O2_glider.days >= datenum(2016,03,10) & N_table.O2_glider.days <= datenum(2016,03,19);
N_table.O2_glider.M1019.days = N_table.O2_glider.days(check);
N_table.O2_glider.M1019.FADV = nanmean(N_table.O2_glider.ADV(check));
N_table.O2_glider.M1019.FADV_error = nanmean(N_table.O2_glider.ADV_error(check));
N_table.O2_glider.M1019.FASE = nanmean(N_table.O2_glider.ASE(check));
N_table.O2_glider.M1019.FASE_error = nanmean(N_table.O2_glider.ASE_error(check));
N_table.O2_glider.M1019.NCP = nanmean(N_table.O2_glider.NCP(check));
N_table.O2_glider.M1019.NCP_error = nanmean(N_table.O2_glider.NCP_error(check));
N_table.O2_glider.M1019.NCP_no_adv = nanmean(N_table.O2_glider.NCP_no_adv(check));
N_table.O2_glider.M1019.NCP_no_adv_error = nanmean(N_table.O2_glider.NCP_no_adv_error(check));
N_table.O2_glider.M1019.NCP_C = N_table.O2_glider.M1019.NCP/1.45;
N_table.O2_glider.M1019.NCP_error_C = N_table.O2_glider.M1019.NCP_error/1.45;
N_table.O2_glider.M1019.NCP_no_adv_C = N_table.O2_glider.M1019.NCP_no_adv/1.45;
N_table.O2_glider.M1019.NCP_no_adv_error_C = N_table.O2_glider.M1019.NCP_no_adv_error/1.45;

% 20-25 March
check = N_table.O2_glider.days > datenum(2016,03,20) & N_table.O2_glider.days < datenum(2016,03,25);
N_table.O2_glider.M2025.days = N_table.O2_glider.days(check);
N_table.O2_glider.M2025.FADV = nanmean(N_table.O2_glider.ADV(check));
N_table.O2_glider.M2025.FADV_error = nanmean(N_table.O2_glider.ADV_error(check));
N_table.O2_glider.M2025.FASE = nanmean(N_table.O2_glider.ASE(check));
N_table.O2_glider.M2025.FASE_error = nanmean(N_table.O2_glider.ASE_error(check));
N_table.O2_glider.M2025.NCP = nanmean(N_table.O2_glider.NCP(check));
N_table.O2_glider.M2025.NCP_error = nanmean(N_table.O2_glider.NCP_error(check));
N_table.O2_glider.M2025.NCP_no_adv = nanmean(N_table.O2_glider.NCP_no_adv(check));
N_table.O2_glider.M2025.NCP_no_adv_error = nanmean(N_table.O2_glider.NCP_no_adv_error(check));
N_table.O2_glider.M2025.NCP_C = N_table.O2_glider.M2025.NCP/1.45;
N_table.O2_glider.M2025.NCP_error_C = N_table.O2_glider.M2025.NCP_error/1.45;
N_table.O2_glider.M2025.NCP_no_adv_C = N_table.O2_glider.M2025.NCP_no_adv/1.45;
N_table.O2_glider.M2025.NCP_no_adv_error_C = N_table.O2_glider.M2025.NCP_no_adv_error/1.45;


% 29-01 April
check = N_table.O2_glider.days > datenum(2016,03,29) & N_table.O2_glider.days < datenum(2016,04,02);
N_table.O2_glider.M2901.days = N_table.O2_glider.days(check);
N_table.O2_glider.M2901.FADV = nanmean(N_table.O2_glider.ADV(check));
N_table.O2_glider.M2901.FADV_error = nanmean(N_table.O2_glider.ADV_error(check));
N_table.O2_glider.M2901.FASE = nanmean(N_table.O2_glider.ASE(check));
N_table.O2_glider.M2901.FASE_error = nanmean(N_table.O2_glider.ASE_error(check));
N_table.O2_glider.M2901.NCP = nanmean(N_table.O2_glider.NCP(check));
N_table.O2_glider.M2901.NCP_error = nanmean(N_table.O2_glider.NCP_error(check));
N_table.O2_glider.M2901.NCP_no_adv = nanmean(N_table.O2_glider.NCP_no_adv(check));
N_table.O2_glider.M2901.NCP_no_adv_error = nanmean(N_table.O2_glider.NCP_no_adv_error(check));
N_table.O2_glider.M2901.NCP_C = N_table.O2_glider.M2901.NCP/1.45;
N_table.O2_glider.M2901.NCP_error_C = N_table.O2_glider.M2901.NCP_error/1.45;
N_table.O2_glider.M2901.NCP_no_adv_C = N_table.O2_glider.M2901.NCP_no_adv/1.45;
N_table.O2_glider.M2901.NCP_no_adv_error_C = N_table.O2_glider.M2901.NCP_no_adv_error/1.45;

% 10-25 March
check = N_table.O2_glider.days >= datenum(2016,03,10) & N_table.O2_glider.days <= datenum(2016,03,25);
N_table.O2_glider.M1025.days = N_table.O2_glider.days(check);
N_table.O2_glider.M1025.FADV = nanmean(N_table.O2_glider.ADV(check));
N_table.O2_glider.M1025.FADV_error = nanmean(N_table.O2_glider.ADV_error(check));
N_table.O2_glider.M1025.FASE = nanmean(N_table.O2_glider.ASE(check));
N_table.O2_glider.M1025.FASE_error = nanmean(N_table.O2_glider.ASE_error(check));
N_table.O2_glider.M1025.NCP = nanmean(N_table.O2_glider.NCP(check));
N_table.O2_glider.M1025.NCP_error = nanmean(N_table.O2_glider.NCP_error(check));
N_table.O2_glider.M1025.NCP_no_adv = nanmean(N_table.O2_glider.NCP_no_adv(check));
N_table.O2_glider.M1025.NCP_no_adv_error = nanmean(N_table.O2_glider.NCP_no_adv_error(check));
N_table.O2_glider.M1025.NCP_C = N_table.O2_glider.M1025.NCP/1.45;
N_table.O2_glider.M1025.NCP_error_C = N_table.O2_glider.M1025.NCP_error/1.45;
N_table.O2_glider.M1025.NCP_no_adv_C = N_table.O2_glider.M1025.NCP_no_adv/1.45;
N_table.O2_glider.M1025.NCP_no_adv_error_C = N_table.O2_glider.M1025.NCP_no_adv_error/1.45;

% 10-03 Apr
check = N_table.O2_glider.days >= datenum(2016,03,10) & N_table.O2_glider.days <= datenum(2016,04,03);
N_table.O2_glider.M1003.days = N_table.O2_glider.days(check);
N_table.O2_glider.M1003.FADV = nanmean(N_table.O2_glider.ADV(check));
N_table.O2_glider.M1003.FADV_error = nanmean(N_table.O2_glider.ADV_error(check));
N_table.O2_glider.M1003.FASE = nanmean(N_table.O2_glider.ASE(check));
N_table.O2_glider.M1003.FASE_error = nanmean(N_table.O2_glider.ASE_error(check));
N_table.O2_glider.M1003.NCP = nanmean(N_table.O2_glider.NCP(check));
N_table.O2_glider.M1003.NCP_error = nanmean(N_table.O2_glider.NCP_error(check));
N_table.O2_glider.M1003.NCP_no_adv = nanmean(N_table.O2_glider.NCP_no_adv(check));
N_table.O2_glider.M1003.NCP_no_adv_error = nanmean(N_table.O2_glider.NCP_no_adv_error(check));
N_table.O2_glider.M1003.NCP_C = N_table.O2_glider.M1003.NCP/1.45;
N_table.O2_glider.M1003.NCP_error_C = N_table.O2_glider.M1003.NCP_error/1.45;
N_table.O2_glider.M1003.NCP_no_adv_C = N_table.O2_glider.M1003.NCP_no_adv/1.45;
N_table.O2_glider.M1003.NCP_no_adv_error_C = N_table.O2_glider.M1003.NCP_no_adv_error/1.45;

%% DIC

N_table.DIC_glider.ADV = [DIC_adv(2:end-1).adv];
N_table.DIC_glider.ASE = [DIC_ase.FDIC(2:end-1)];
N_table.DIC_glider.ADV_error = [errors.ADV_DIC(2:end-1).errors_adv];
N_table.DIC_glider.ASE_error = [errors.ASE.errors_ASE_DIC(2:end-1)];
N_table.DIC_glider.NCP = NCP_est_DIC;
N_table.DIC_glider.NCP_no_adv = NCP_est_no_adv_DIC;
N_table.DIC_glider.NCP_error = [errors.error_NCP_DIC];
N_table.DIC_glider.NCP_no_adv_error = [errors.error_NCP_no_adv_DIC];
N_table.DIC_glider.days = datenum(2016,03,10):1:datenum(2016,04,03);

% 10-19 March
check = N_table.DIC_glider.days >= datenum(2016,03,10) & N_table.DIC_glider.days <= datenum(2016,03,19);
N_table.DIC_glider.M1019.days = N_table.DIC_glider.days(check);
N_table.DIC_glider.M1019.FADV = nanmean(N_table.DIC_glider.ADV(check));
N_table.DIC_glider.M1019.FADV_error = nanmean(N_table.DIC_glider.ADV_error(check));
N_table.DIC_glider.M1019.FASE = nanmean(N_table.DIC_glider.ASE(check));
N_table.DIC_glider.M1019.FASE_error = nanmean(N_table.DIC_glider.ASE_error(check));
N_table.DIC_glider.M1019.NCP = nanmean(N_table.DIC_glider.NCP(check));
N_table.DIC_glider.M1019.NCP_error = nanmean(N_table.DIC_glider.NCP_error(check));
N_table.DIC_glider.M1019.NCP_no_adv = nanmean(N_table.DIC_glider.NCP_no_adv(check));
N_table.DIC_glider.M1019.NCP_no_adv_error = nanmean(N_table.DIC_glider.NCP_no_adv_error(check));
N_table.DIC_glider.M1019.NCP_O2 = N_table.DIC_glider.M1019.NCP*1.45;
N_table.DIC_glider.M1019.NCP_error_O2 = N_table.DIC_glider.M1019.NCP_error*1.45;
N_table.DIC_glider.M1019.NCP_no_adv_O2 = N_table.DIC_glider.M1019.NCP_no_adv*1.45;
N_table.DIC_glider.M1019.NCP_no_adv_error_O2 = N_table.DIC_glider.M1019.NCP_no_adv_error*1.45;

% 20-25 March
check = N_table.DIC_glider.days > datenum(2016,03,20) & N_table.DIC_glider.days < datenum(2016,03,25);
N_table.DIC_glider.M2025.days = N_table.DIC_glider.days(check);
N_table.DIC_glider.M2025.FADV = nanmean(N_table.DIC_glider.ADV(check));
N_table.DIC_glider.M2025.FADV_error = nanmean(N_table.DIC_glider.ADV_error(check));
N_table.DIC_glider.M2025.FASE = nanmean(N_table.DIC_glider.ASE(check));
N_table.DIC_glider.M2025.FASE_error = nanmean(N_table.DIC_glider.ASE_error(check));
N_table.DIC_glider.M2025.NCP = nanmean(N_table.DIC_glider.NCP(check));
N_table.DIC_glider.M2025.NCP_error = nanmean(N_table.DIC_glider.NCP_error(check));
N_table.DIC_glider.M2025.NCP_no_adv = nanmean(N_table.DIC_glider.NCP_no_adv(check));
N_table.DIC_glider.M2025.NCP_no_adv_error = nanmean(N_table.DIC_glider.NCP_no_adv_error(check));
N_table.DIC_glider.M2025.NCP_O2 = N_table.DIC_glider.M2025.NCP*1.45;
N_table.DIC_glider.M2025.NCP_error_O2 = N_table.DIC_glider.M2025.NCP_error*1.45;
N_table.DIC_glider.M2025.NCP_no_adv_O2 = N_table.DIC_glider.M2025.NCP_no_adv*1.45;
N_table.DIC_glider.M2025.NCP_no_adv_error_O2 = N_table.DIC_glider.M2025.NCP_no_adv_error*1.45;


% 10-25 March
check = N_table.DIC_glider.days >= datenum(2016,03,10) & N_table.DIC_glider.days <= datenum(2016,03,25);
N_table.DIC_glider.M1025.days = N_table.DIC_glider.days(check);
N_table.DIC_glider.M1025.FADV = nanmean(N_table.DIC_glider.ADV(check));
N_table.DIC_glider.M1025.FADV_error = nanmean(N_table.DIC_glider.ADV_error(check));
N_table.DIC_glider.M1025.FASE = nanmean(N_table.DIC_glider.ASE(check));
N_table.DIC_glider.M1025.FASE_error = nanmean(N_table.DIC_glider.ASE_error(check));
N_table.DIC_glider.M1025.NCP = nanmean(N_table.DIC_glider.NCP(check));
N_table.DIC_glider.M1025.NCP_error = nanmean(N_table.DIC_glider.NCP_error(check));
N_table.DIC_glider.M1025.NCP_no_adv = nanmean(N_table.DIC_glider.NCP_no_adv(check));
N_table.DIC_glider.M1025.NCP_no_adv_error = nanmean(N_table.DIC_glider.NCP_no_adv_error(check));
N_table.DIC_glider.M1025.NCP_O2 = N_table.DIC_glider.M1025.NCP*1.45;
N_table.DIC_glider.M1025.NCP_error_O2 = N_table.DIC_glider.M1025.NCP_error*1.45;
N_table.DIC_glider.M1025.NCP_no_adv_O2 = N_table.DIC_glider.M1025.NCP_no_adv*1.45;
N_table.DIC_glider.M1025.NCP_no_adv_error_O2 = N_table.DIC_glider.M1025.NCP_no_adv_error*1.45;

%% buoy NCP

%% DIC

% 20-25 March
N_table.DIC_buoy.M2025.ASE = nanmean(NCP_DIC.M2025.F);
N_table.DIC_buoy.M2025.ASE_error = buoy_error_DIC.ase_20_25;
N_table.DIC_buoy.M2025.ADV = N_table.DIC_glider.M2025.FADV;
N_table.DIC_buoy.M2025.ADV_error = buoy_error_DIC.adv_20_25;
N_table.DIC_buoy.M2025.NCP = NCP_DIC.M2025.NCP;
N_table.DIC_buoy.M2025.NCP_error = buoy_error_DIC.error_20_25;
N_table.DIC_buoy.M2025.NCP_ADV = NCP_DIC.M2025.NCP_ADV;
N_table.DIC_buoy.M2025.NCP_ADV_error = buoy_error_DIC.error_20_25_adv;
N_table.DIC_buoy.M2025.NCP_O2 = N_table.DIC_buoy.M2025.NCP*1.45;
N_table.DIC_buoy.M2025.NCP_error_O2 = N_table.DIC_buoy.M2025.NCP_error*1.45;
N_table.DIC_buoy.M2025.NCP_ADV_O2 = N_table.DIC_buoy.M2025.NCP_ADV*1.45;
N_table.DIC_buoy.M2025.NCP_ADV_error_O2 = N_table.DIC_buoy.M2025.NCP_ADV_error*1.45;

% 29-01 April
N_table.DIC_buoy.M2901.ASE = nanmean(NCP_DIC.M2901.F);
N_table.DIC_buoy.M2901.ASE_error = buoy_error_DIC.ase_29_01;
% N_table.DIC_buoy.ADV = N_table.DIC_glider.M2901.FADV;
% N_table.DIC_buoy.ADV_error = buoy_error_DIC.adv_29_01;
N_table.DIC_buoy.M2901.NCP = NCP_DIC.M2901.NCP;
N_table.DIC_buoy.M2901.NCP_error = buoy_error_DIC.error_29_01;
% N_table.DIC_buoy.NCP_ADV = NCP_DIC.M2901.NCP_ADV;
% N_table.DIC_buoy.NCP_ADV_error = buoy_error_DIC.error_29_25_01;
N_table.DIC_buoy.M2901.NCP_O2 = N_table.DIC_buoy.M2901.NCP*1.45;
N_table.DIC_buoy.M2901.NCP_error_O2 = N_table.DIC_buoy.M2901.NCP_error*1.45;
% N_table.DIC_buoy.M2025.NCP_ADV_O2 = N_table.DIC_buoy.M2025.NCP_ADV*1.45;
% N_table.DIC_buoy.M2025.NCP_ADV_error_O2 = N_table.DIC_buoy.M2025.NCP_ADV_error*1.45;

%% O2

% 20-25 March
N_table.O2_buoy.M2025.ASE = NCP_O2.M2025.ASE;
N_table.O2_buoy.M2025.ASE_error = buoy_error_O2.ase_20_25;
N_table.O2_buoy.M2025.ADV = N_table.O2_glider.M2025.FADV;
N_table.O2_buoy.M2025.ADV_error = buoy_error_O2.adv_20_25;
N_table.O2_buoy.M2025.NCP = NCP_O2.M2025.NCP;
N_table.O2_buoy.M2025.NCP_error = buoy_error_O2.error_20_25;
N_table.O2_buoy.M2025.NCP_ADV = NCP_O2.M2025.NCP_ADV;
N_table.O2_buoy.M2025.NCP_ADV_error = buoy_error_O2.error_20_25_adv;
N_table.O2_buoy.M2025.NCP_DIC = N_table.O2_buoy.M2025.NCP/1.45;
N_table.O2_buoy.M2025.NCP_error_DIC = N_table.O2_buoy.M2025.NCP_error/1.45;
N_table.O2_buoy.M2025.NCP_ADV_DIC = N_table.O2_buoy.M2025.NCP_ADV/1.45;
N_table.O2_buoy.M2025.NCP_ADV_error_DIC = N_table.O2_buoy.M2025.NCP_ADV_error/1.45;

% 29-01 April
N_table.O2_buoy.M2901.ASE = NCP_O2.M2901.ASE;
N_table.O2_buoy.M2901.ASE_error = buoy_error_O2.ase_29_01;
N_table.O2_buoy.M2901.ADV = N_table.O2_glider.M2901.FADV;
N_table.O2_buoy.M2901.ADV_error = buoy_error_O2.adv_29_01;
N_table.O2_buoy.M2901.NCP = NCP_O2.M2901.NCP;
N_table.O2_buoy.M2901.NCP_error = buoy_error_O2.error_29_01;
N_table.O2_buoy.M2901.NCP_ADV = NCP_O2.M2901.NCP_ADV;
N_table.O2_buoy.M2901.NCP_ADV_error = buoy_error_O2.error_29_01_adv;
N_table.O2_buoy.M2901.NCP_DIC = N_table.O2_buoy.M2901.NCP/1.45;
N_table.O2_buoy.M2901.NCP_error_DIC = N_table.O2_buoy.M2901.NCP_error/1.45;
N_table.O2_buoy.M2901.NCP_ADV_DIC = N_table.O2_buoy.M2901.NCP_ADV/1.45;
N_table.O2_buoy.M2901.NCP_ADV_error_DIC = N_table.O2_buoy.M2901.NCP_ADV_error/1.45;
