function orbit_data = ...
    read_GMAT_eph(file_directory, file_name, exact_time_step)
%% DESCRITION 
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Start date:           April 7, 2016
%       Last Modified:        April 7, 2016
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% This function reads in the orbit data from the NASA GMAT ephemeris files.
%
% -------------------------------------------------------------------------
% INPUT:
%
%   file_directory  = Directory of the GMAT ephemeris file. This is a
%                     string that ends with '.../'
%   file_name       = File name of the GMAT *.eph file.
%   exact_time_step = Boolean that flags whether the ephemeris file was
%                     produced at exact time intervals. Matlab's datenum
%                     function tends to be inaccurate. This rounds the time
%                     step to an exact value.
%
% -------------------------------------------------------------------------
% OUTPUT
%
%   orbit_data = Data structure of the form:
%               orbit_data.pos_m   = position vectors in meters.
%               orbit_data.vel_m_s = velocity vectors in meters / sec.
%               orbit_data.datestr = raw data time stamp.
%               orbit_data.datenum = matlab time stamp date num.
%               orbit_data.elapsed_time_sec = data elapsed time in sec.
%
%% IMPLEMENTATION

% Set deliminiters for reading the data file.
delim = ' ';
nheaderlines = 20;
orbit_data_raw = ...
    importdata([file_directory, file_name], delim, nheaderlines);

% Parse position and velocity data and set to the correct units.
orbit_data.pos_m = orbit_data_raw.data(:,1:3) * 1000;
orbit_data.vel_m_s = orbit_data_raw.data(:,4:6) * 1000;

% Format the date / time vectors.
orbit_data.datestr = orbit_data_raw.textdata(18:end,1);
orbit_data.datenum = datenum(orbit_data.datestr,'yyyy-mm-ddTHH:MM:SS.FFF');

% Make an elapsed time vector.
elapsed_time_sec = ...
    ( orbit_data.datenum - orbit_data.datenum(1) ) * 24 * 3600; % [sec]

% Make elapsed second vector exact since the matlab datenum function seems
% to round and loose precision.
if exact_time_step
    dt = round( median( diff(elapsed_time_sec) ) ); % [sec]

    elapsed_time_sec = ...
        [0:dt:( orbit_data.datenum(end) - ...
        orbit_data.datenum(1) ) * 24 * 3600]';
end

% Set the final elapsed time vector.
orbit_data.elapsed_time_sec = elapsed_time_sec; 

