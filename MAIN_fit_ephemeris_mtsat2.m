%% MAIN SCRIPT FOR THE SBAS L5 MOPS EPHEMERIS FIT ANALYSIS.
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Date:                 May  2, 2017
%       Last Modified:        February 25, 2018
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% This fits the L5 SBAS MOPS ephemeris message parameters to precision
% orbit data. It also performs fit error analysis and evaluates the message
% performance. Specificially, this looks at the corner cases that cause
% problems with the fitting algorithm convergence. For more info, please
% see Appendix B of:
%
%       T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite
%           Systems," Doctor of Philosophy, Aeronautics and Astronautics,
%           Stanford University, Stanford, CA, 2017.
%
% This is available from: https://purl.stanford.edu/dc409wn9227
%
% This particular version is for trials with MTSAT2.
%
%% SET UP WORKSPACE

clc
clear
close all

% Load physical constants file to enter them into the global workspace.
physical_constants_GPS

% Define global variables.
global omega_e mu R_e

% Turn rank deficient warning off.
warning('off','MATLAB:rankDeficientMatrix');

%% INITIAL CONDITIONS USED IN SIMULATION [For reference]

SMA = 42164.13866; % [km]
ecc = 0.00042061; % [-]
inc = 0.00195; % [deg]
RAAN = 0.391565; % [deg]
AOP = 327.12287; % [deg]
MA = 304.10548; % [deg]

%% IMPLEMENTATION

% Read in GMAT ephemeris file. 
file_directory_GMAT = [pwd, '/Results_Fitting_MTSAT2/'];
file_name_GMAT = 'EphemerisFile_mtsat2_1day.eph';

% Read data in the propagated orbit file.
exact_time_step = true;
orbit_data = ...
    read_GMAT_eph(file_directory_GMAT, file_name_GMAT, exact_time_step);

% Fit interval of interest.
fit_interval = 4 * 60; % [seconds]

% File names.
file_altitudes = SMA * 1000 - R_e; % [km]

% Directory for saving the data.
file_dir_save = [pwd, '/Results_Fitting_MTSAT2/'];
 
% Total time in file.
time_in_file = orbit_data.elapsed_time_sec(end); % [sec] 

% Start the clock.
tic;

% Get the analytical coefficients for the URE equations.
[coeff_A2, coeff_C2, coeff_R2, coeff_T2, coeff_RT, theta] = ...
    analytic_URE_eqn(  file_altitudes, R_e );

% Time between messages to be fit.
time_between_messages =  10; % [sec]

% Determine the message start times.
message_start_times = ...
    0:time_between_messages:(time_in_file - fit_interval);

% Determine the number of messages that can be made per file.
num_eph_per_file = length( message_start_times );

% Initialize variables.
rms_ure_Save = NaN(num_eph_per_file, 1);
rms_3D_Save = NaN(num_eph_per_file, 1);
max_3D_Save = NaN(num_eph_per_file, 1);
Num_Iter_Save = NaN(num_eph_per_file, 1);
convergence_crit_Save = NaN(num_eph_per_file, 1);
failure_flag_Save = zeros(num_eph_per_file, 1);
fit_type_Save = zeros(num_eph_per_file, 1);
eph_datenum_Save = NaN(num_eph_per_file, 1);

eph_Save(num_eph_per_file).Asqrt  = [];
eph_Save(num_eph_per_file).e      =[];
eph_Save(num_eph_per_file).i0     =[];
eph_Save(num_eph_per_file).Omega0 =[];
eph_Save(num_eph_per_file).Omega  =[];
eph_Save(num_eph_per_file).M0     =[];

eph_Save(num_eph_per_file).Cus    =[];
eph_Save(num_eph_per_file).Cuc    =[];
eph_Save(num_eph_per_file).Crs    =[];
eph_Save(num_eph_per_file).Crc    =[];
eph_Save(num_eph_per_file).Cis    =[];
eph_Save(num_eph_per_file).Cic    =[];

eph_Save(num_eph_per_file).IDOT      =[];
eph_Save(num_eph_per_file).Omega_dot =[];
eph_Save(num_eph_per_file).Delta_n   =[];
eph_Save(num_eph_per_file).Toe   =[];

for idx_message = 1:num_eph_per_file
    % Define the time vector for fitting.
    time = ...
        orbit_data.elapsed_time_sec(...
        orbit_data.elapsed_time_sec <= fit_interval );
    
    % Get the start / end index for the data to fit to.
    idx_start = find(orbit_data.elapsed_time_sec == ...
        message_start_times(idx_message));
    idx_end = find(orbit_data.elapsed_time_sec == ...
        message_start_times(idx_message) + fit_interval);
    
    % Get position and velocity vectors
    pos = ...
        orbit_data.pos_m(idx_start:idx_end,:);
    vel = ...
        orbit_data.vel_m_s(idx_start:idx_end,:);
    
    % Convert ECEF to ECI coordinates.
    % NOTE: We'll work in ECI coordinates for the purposes of this
    %       experiement but, ECEF is also possibele here with some
    %       small changes.
    % theta_g = utc2gmst( datevec(orbit_data.datenum(idx_start)) ); % [rad]
    % [R_test, V_test] = ECEF2ECI(...
    %     pos(1,:)',  vel(1,:)', theta_g)
    
    % Since we're dealing with ECI vs ECEF coordinates, we'll use a
    % zero offset between them (for the purposes of reusing other
    % code).
    theta_g = 0; % [rad]
    
    % Form initial guess with the 6 Keplerian elements.
    [coe, undefined, orbit_type] = ECI2COE(pos(1,:), vel(1,:));
    
    % Form the initial guess for the estimator.
    a = coe.a; % [m]
    n = sqrt(mu/a^3); % [rad/sec]
    ecc = coe.e; % [-]
    inc = coe.i * pi / 180; % [rad]
    RAAN = coe.RAAN * pi / 180; % [rad]
    omega = coe.omega * pi /180; % [rad]
    M0 = coe.M * pi / 180; % [rad]
    
    Cus = 0; % [rad]
    Cuc = 0; % [rad]
    Crs = 0; % [rad]
    Crc = 0; % [rad]
    Cis = 0; % [rad]
    Cic = 0; % [rad]
    
    IDOT = 0; % [rad/s]
    OMEGA_DOT = 0; % [rad/s]
    delta_n = 0; % [rad/s]
    
    % Form the initial guess.
    initial_guess(1) = a;
    initial_guess(2) = ecc;
    initial_guess(3) = inc;
    initial_guess(4) = RAAN;
    initial_guess(5) = omega;
    initial_guess(6) = M0;
    
    initial_guess(7)  = Cus;
    initial_guess(8)  = Cuc;
    initial_guess(9)  = Crs;
    initial_guess(10) = Crc;
    initial_guess(11) = Cis;
    initial_guess(12) = Cic;
    
    initial_guess(13) = IDOT;
    initial_guess(14) = OMEGA_DOT;
    initial_guess(15) = delta_n;
    
    % Define the convergence criteria.
    ConvCrit = 1e-11;
    
    % Fit L5 SBAS MOPS ephemeris parameters or subset.
    % Define the weighting matrix, use the identity matrix for now.
    Wmat = ones( size(time) );
    
    fit_parameters = zeros(1,15);
    
    % Keplerian Elements.
    fit_parameters(1) = 1; % a
    fit_parameters(2) = 1; % e
    fit_parameters(3) = 1; % inc
    fit_parameters(4) = 1; % RAAN
    fit_parameters(5) = 1; % omega
    fit_parameters(6) = 1; % M0
    
    % Corrections.
    fit_parameters(7) = 1; % Cus
    fit_parameters(8) = 1; % Cuc
    fit_parameters(13) = 1; % IDOT
    
    % Fit parameters.
    [a, ecc, inc, RAAN, omega, M0,...
        Cus, Cuc, Crs, Crc, Cis, Cic, ...
        IDOT, OMEGA_DOT, delta_n, flag, NumIter, fit_type] = ...
        COE15_estimator_wrapper(time, pos, vel, initial_guess, ...
        Wmat, ConvCrit, fit_parameters, ...
        theta_g, coeff_R2, coeff_A2, coeff_C2);
    
    % Quantize message parameters.
    [a, ecc, inc, RAAN, omega, M0, Cus, Cuc, IDOT, numbits] = ...
        bit_reduction(a, ecc, inc, RAAN, omega, M0, Cus, Cuc, IDOT);
    
    % Error analysis.
    [error_3d, rms_error, rms_ure, ...
        error_radial, error_along_track, error_cross_track, ...
        eph] = eph_error_analysis(sqrt(a), ecc, inc, RAAN, omega, M0, ...
        Cus, Cuc, Crc, Crs, Cic, Cis, ...
        IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
        coeff_R2, coeff_A2, coeff_C2);
    
    % Save results.
    rms_ure_Save(idx_message) = rms_ure;
    rms_3D_Save(idx_message) = rms_error;
    max_3D_Save(idx_message) = max(error_3d);
    Num_Iter_Save(idx_message) = NumIter;
    convergence_crit_Save(idx_message) =  ...
        ConvCrit;
    failure_flag_Save(idx_message) = flag;
    fit_type_Save(idx_message) = fit_type;
    eph_datenum_Save(idx_message) = orbit_data.datenum(1) + ...
        message_start_times(idx_message) / 24 / 3600;
    
    eph_Save(idx_message) = eph;
end % end idx_message

% Give brief summary of performance.
disp(['Total failures is: ', num2str(sum(failure_flag_Save(:)))])
disp(['Median RMS Error is: ', num2str(median(rms_3D_Save(:)))])
disp(['95th Percentile RMS Error is: ', num2str( prctile(rms_3D_Save(:),95)) ])

% Save a *.mat file.
file_name_save = [file_dir_save, 'Results_MTSAT2.mat'];
save(file_name_save, ...
    'rms_ure_Save', 'rms_3D_Save', 'max_3D_Save',...
    'Num_Iter_Save', 'convergence_crit_Save', 'failure_flag_Save', ...
    'eph_datenum_Save',...
    'eph_Save', 'fit_type_Save', 'orbit_data');

figure;
plot(rms_ure_Save)

% Clear saved variables to avoid confusion.
clear('rms_ure_Save', 'rms_3D_Save','max_3D_Save', ...
    'Num_Iter_Save', 'convergence_crit_Save', 'failure_flag_Save', ...
    'eph_datenum_Save',...
    'eph_Save', 'fit_type_Save', 'orbit_data');

% Turn rank deficient warning on.
warning('on','MATLAB:rankDeficientMatrix');

%% ANALYZE RESULTS

% Close all plot. 
close all

% Select data file. 
file_name_save = [file_dir_save, 'Results_MTSAT2.mat'];

% Set date limits for plotting. 
x_lower_plot = datenum('Feb 03, 2018 23:30:00');
x_upper_plot = datenum('Feb 04, 2018 23:30:00');

% Load data. 
load(file_name_save);

% Plot RMS URE as a function of date. 
figure; 
plot(eph_datenum_Save, max_3D_Save * 100, 'linewidth', 2)
grid on
datetick
xlim([x_lower_plot, x_upper_plot])
xlabel({'UTC',datestr(orbit_data.datenum(1), 'mmm-dd-yyyy')})
ylabel('Max 3D Error in Fit Interval [cm]')
title('Max 3D Representation Error as a Function of Time')
fileSave = [file_dir_save ,'error_max_3D_vs_time.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);

% Plot RMS 3D error as a function of date. 
figure; 
plot(eph_datenum_Save, rms_3D_Save * 100, 'linewidth', 2)
grid on
datetick
xlim([x_lower_plot, x_upper_plot])
xlabel({'UTC',datestr(orbit_data.datenum(1), 'mmm-dd-yyyy')})
ylabel('3D RMS Error [cm]')
title('3D RMS Representation Error as a Function of Time')
fileSave = [file_dir_save ,'error_rms_3D_vs_time.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);

% Plot RMS URE as a function of date. 
figure; 
plot(eph_datenum_Save, rms_ure_Save * 100, 'linewidth', 2)
grid on
datetick
xlim([x_lower_plot, x_upper_plot])
xlabel({'UTC',datestr(orbit_data.datenum(1), 'mmm-dd-yyyy')})
ylabel('RMS URE [cm]')
title('RMS Representation URE as a Function of Time')
fileSave = [file_dir_save ,'error_rms_ure_vs_time.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);

% Ouput performance statistics to file. 
file_summary = [file_dir_save, 'summary.txt'];
file_id = fopen(file_summary, 'w');

fprintf(file_id, ['Total failures is: ', num2str(sum(failure_flag_Save(:))), '\n']);
fprintf(file_id, '\n');
fprintf(file_id, ['Median RMS Error is: ', num2str(median(rms_3D_Save(:))), ' [m]\n']);
fprintf(file_id, ['95th Percentile RMS Error is: ', num2str( prctile(rms_3D_Save(:),95)), ' [m]\n']);
fprintf(file_id, ['Max RMS Error is: ', num2str( max(rms_3D_Save(:))), ' [m]\n']);
fprintf(file_id, '\n');
fprintf(file_id, ['Median RMS URE is: ', num2str(median(rms_ure_Save(:))), ' [m]\n']);
fprintf(file_id, ['95th Percentile RMS URE is: ', num2str( prctile(rms_ure_Save(:),95)), ' [m]\n']);
fprintf(file_id, ['Max RMS URE is: ', num2str( max(rms_ure_Save(:))) ,' [m]\n']);
fprintf(file_id, '\n');
fprintf(file_id, ['Max In-Fit-Interval Error is: ', num2str( prctile(max_3D_Save(:),95)), ' [m]\n' ]);

fclose(file_id);

%% PLOT RAAN FOR COMPARISON

% Load 7 day orbit propagation data set. 
file_directory_GMAT = [pwd, '/Results_Fitting_MTSAT2/'];
file_name_GMAT = 'EphemerisFile_mtsat2_7day.eph';

% Read data in the propagated orbit file.
exact_time_step = true;
orbit_data = ...
    read_GMAT_eph(file_directory_GMAT, file_name_GMAT, exact_time_step);

% Plot RAAN as a function of date based on Keplerian model of GMAT data. 
RAAN_propagated = length(orbit_data.pos_m);
for idx = 1:length(orbit_data.pos_m)
    [coe, ~, ~] = ...
        ECI2COE(orbit_data.pos_m(idx, :), orbit_data.vel_m_s(idx, :));
    
    % Save RAAN. 
    RAAN_propagated(idx) = coe.RAAN;
end

figure; 
plot(orbit_data.elapsed_time_sec/24/60/60, RAAN_propagated, 'linewidth', 2)
grid on
ylabel('RAAN [deg]')
xlabel(['Propagation Time Since ', orbit_data.datestr(1) ,' [days]'])
title('RAAN as a function of propation time')
fileSave = [file_dir_save ,'RAAN_vs_time.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);

