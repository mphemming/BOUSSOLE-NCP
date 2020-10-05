function [mooring] = get_mooring(filename,imos_time)
%
% function: get_mooring.m
%
% created by MPH @ UNSW Sydney, Australia on 27/02/2020
% last updated on 29/04/2020
% m.hemming@unsw.edu.au
% This script was created using MATLAB version 9.8.0.1323502 (R2020a)
%
% Input: 
%
% filename     |     string input including path and netCDF filename including extension
% imos_time    |     0/1 => no/yes to convert netCDF IMOS time (days since 1950-01-01) to MATLAB datenum 
%
% Output:
%
% mooring      |     structure file containing variables included in netCDF, plus file attributes and information

%% get information
mooring.file_info = ncinfo(filename);

%% obtain all variables in file

for n_var = 1:numel(mooring.file_info.Variables)
    % is datatype supported?
    if isempty(strmatch(mooring.file_info.Variables(n_var).Datatype,'UNSUPPORTED DATATYPE'))  
        eval(['mooring.',mooring.file_info.Variables(n_var).Name,' = ncread(filename,mooring.file_info.Variables(n_var).Name);']);
    else
        display(['UNSUPPORTED VARIABLE, n = ',num2str(n_var)])
        continue
    end
end

%% convert from IMOS time to MATLAB datenum if required

if imos_time == 1
    mooring.TIME = datenum(1950,01,01) + mooring.TIME;
end

end
