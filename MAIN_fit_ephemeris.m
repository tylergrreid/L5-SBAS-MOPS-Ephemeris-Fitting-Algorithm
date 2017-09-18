%% MAIN SCRIPT FOR THE SBAS L5 MOPS EPHEMERIS FIT ANALYSIS.
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Date:                 May  2, 2017
%       Updated:              May 24, 2017
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

% Load the info about the analysis.
load('orbit_numbering_GMAT.mat')

%% IMPLEMENTATION

% The RAANs to run through.
RAAN_cases = [0, 90, 180, 270];

% Fit interval of interest.
fit_interval = 4 * 60; % [seconds]

% File names.
file_altitudes = 35785; % [km]

% Directory for the propagated orbit data.
file_directory = [pwd, '/GMAT_precision_ephemeris_files/'];

% Directory for saving the data.
file_dir_save = [pwd, '/Results_Fitting/'];

% Total time in file.
time_in_file = 24 * 60 * 60; % [sec] This is 1 day.

% String for ease in parsing filenames.
preamble = '000';

% Max number of files.
max_file_num = num_cases;

% Start the clock.
tic;

% Quantize results flag. If set true, then the final output will be
% quantized as per the GPS message bits / scale factors.
quantize_flag = false;
if quantize_flag
    file_dir_save = file_dir_save_quantized;
end

% Get the analytical coefficients for the URE equations.
[coeff_A2, coeff_C2, coeff_R2, coeff_T2, coeff_RT, theta] = ...
    analytic_URE_eqn(  altitudes * 1000, R_e );

% Each RAAN will save it's own *.mat file.
for idx_RAAN = 1:length(RAAN_cases)
    
    % Time between messages to be fit.
    time_between_messages = 15 * 60; % [sec]
    
    % Determine the message start times.
    message_start_times = ...
        0:time_between_messages:(time_in_file - fit_interval);
    
    % Determine the number of messages that can be made per file.
    num_eph_per_file = length( message_start_times );
    
    % Determine the total number of messages.
    num_eph_total = num_eph_per_file * max_file_num;
    
    % Initialize variables.
    rms_ure_Save = NaN(num_eph_per_file, num_cases);
    rms_3D_Save = NaN(num_eph_per_file, num_cases);
    Num_Iter_Save = NaN(num_eph_per_file, num_cases);
    convergence_crit_Save = NaN(num_eph_per_file, num_cases);
    failure_flag_Save = zeros(num_eph_per_file, num_cases);
    fit_type_Save = zeros(num_eph_per_file, num_cases);
        
    eph_Save(num_eph_per_file, num_cases).Asqrt  = [];
    eph_Save(num_eph_per_file, num_cases).e      =[];
    eph_Save(num_eph_per_file, num_cases).i0     =[];
    eph_Save(num_eph_per_file, num_cases).Omega0 =[];
    eph_Save(num_eph_per_file, num_cases).Omega  =[];
    eph_Save(num_eph_per_file, num_cases).M0     =[];
    
    eph_Save(num_eph_per_file, num_cases).Cus    =[];
    eph_Save(num_eph_per_file, num_cases).Cuc    =[];
    eph_Save(num_eph_per_file, num_cases).Crs    =[];
    eph_Save(num_eph_per_file, num_cases).Crc    =[];
    eph_Save(num_eph_per_file, num_cases).Cis    =[];
    eph_Save(num_eph_per_file, num_cases).Cic    =[];
    
    eph_Save(num_eph_per_file, num_cases).IDOT      =[];
    eph_Save(num_eph_per_file, num_cases).Omega_dot =[];
    eph_Save(num_eph_per_file, num_cases).Delta_n   =[];
    eph_Save(num_eph_per_file, num_cases).Toe   =[];
    
    % Run through the orbit scenarios (eccentricities, inclinations)
    for idx_orbit_scenario = 1:num_cases
        % Give us a sense of where we are in processing (output).
        fprintf(['Starting RAAN = ',num2str(RAAN_cases(idx_RAAN)),...
            ' [deg] and scenario ', num2str( idx_orbit_scenario ),...
            ' out of ',...
            num2str(num_cases),...
            ', Elapsed time is ',num2str(toc/3600),' [hours] \n'])
        
        for idx_message = 1:num_eph_per_file
            % Get the propagated orbit data to be fit.
            % Get the file name.
            idx_zero = length(num2str(idx_orbit_scenario));
            file_name = horzcat('EphemerisFile_MOPSSat_alt_35785_km_RAAN_',...
                num2str(RAAN_cases(idx_RAAN)), '_', preamble(1:end-idx_zero),...
                num2str(idx_orbit_scenario),'.eph');
            
            % Read data in the file.
            exact_time_step = true;
            orbit_data = ...
                read_GMAT_eph(file_directory, file_name, exact_time_step);
            
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
            rms_ure_Save(idx_message, idx_orbit_scenario) = rms_ure;
            rms_3D_Save(idx_message, idx_orbit_scenario) = rms_error;
            Num_Iter_Save(idx_message, idx_orbit_scenario) = NumIter;
            convergence_crit_Save(idx_message, idx_orbit_scenario) =  ...
                ConvCrit;
            failure_flag_Save(idx_message, idx_orbit_scenario) = flag;
            fit_type_Save(idx_message, idx_orbit_scenario) = fit_type;

            eph_Save(idx_message, idx_orbit_scenario) = eph;
        end % end idx_message
    end % end idx_orbit_scenario
    
    % Give brief summary of performance.
    disp(['Total failures is: ', num2str(sum(failure_flag_Save(:)))])
    disp(['Median RMS Error is: ', num2str(median(rms_3D_Save(:)))])
    disp(['95th Percentile RMS Error is: ', num2str( prctile(rms_3D_Save(:),95)) ])
    
    % Save a *.mat file.
    file_name_save = [file_dir_save, 'Results_RAAN_',...
        num2str(RAAN_cases(idx_RAAN)),'.mat'];
    save(file_name_save, ...
        'rms_ure_Save', 'rms_3D_Save', ...
        'Num_Iter_Save', 'convergence_crit_Save', 'failure_flag_Save', ...
        'eph_Save', 'fit_type_Save');
    
    % Clear saved variables to avoid confusion.
    clear('rms_ure_Save', 'rms_3D_Save', ...
        'Num_Iter_Save', 'convergence_crit_Save', 'failure_flag_Save', ...
        'eph_Save', 'fit_type_Save');
    
end % end idx_RAAN

% Turn rank deficient warning on.
warning('on','MATLAB:rankDeficientMatrix');

%% ANALYZE RESULTS

close all

% Initialize arrays.
error_95_3D_array = NaN(length(ecc_cases), length(inc_cases));
error_95_ure_array = NaN(length(ecc_cases), length(inc_cases));
error_max_3D_array = NaN(length(ecc_cases), length(inc_cases));
error_max_ure_array = NaN(length(ecc_cases), length(inc_cases));
fit_type_array = zeros(length(ecc_cases), length(inc_cases));
fit_type_p_array = NaN(length(ecc_cases), length(inc_cases));

rms_ure_All = [];
rms_3D_All = [];
fit_type_All = [];

% Gather the data from each scenario.
for idx_RAAN = 1:length(RAAN_cases)
    % Load files and get data.
    file_name_save = [file_dir_save, 'Results_RAAN_',...
        num2str(RAAN_cases(idx_RAAN)),'.mat'];
    
    load(file_name_save);
    
    % Get all the data. 
    rms_ure_All = [rms_ure_All; rms_ure_Save];
    rms_3D_All = [rms_3D_All; rms_3D_Save];
    fit_type_All = [fit_type_All; fit_type_Save];
end

% Initialize case number
case_number = 1;

% Make RMS error into a grid corresponding to eccentricity and
% inclinations tested.
for ecc_idx = 1:length(ecc_cases)
    for inc_idx = 1:length(inc_cases)
        % Get the data for a given eccentricity and inclination. 

        % Get the 95% error.
        error_95_3D_array(inc_idx, ecc_idx) = ...
            prctile( rms_3D_All(:, case_number),95 );

        % Get 95th percentile ure.
        error_95_ure_array(inc_idx, ecc_idx) = ...
            prctile( rms_ure_All(:, case_number), 95 );

        % Get the max 3D error.
        error_max_3D_array(inc_idx, ecc_idx) = ...
            max( rms_3D_All(:, case_number) );
        
        % Get the max 3D error.
        error_max_ure_array(inc_idx, ecc_idx) = ...
            max( rms_ure_All(:, case_number) );
        
        % Get the fit method.
        fits = fit_type_All(:, case_number);

        fit_type_array(inc_idx, ecc_idx) = ...
            mode( fits(fits~=1) );

        % Get percentage of cases that need this.
        fit_type_p_array(inc_idx, ecc_idx) = ...
            length(fits(fits~=1))/length(fits);

        % Update count.
        case_number = case_number + 1;
    end
end
      
%% PLOT RESULTS

% Plot contour of results ure.
figure;
[cs, hc] = contourf(error_95_ure_array*100);
colorbar
xlim([1,13])
ylim([1,13])
% caxis([0.9, 1.8])
caxis([1.0, 3.0])
h = colorbar;
ylabel(h, '95% RMS URE [cm]')
set(hc, 'EdgeColor','none')
colormap parula;
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'RMS_URE_95_Error_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);


% Plot contour of results ure.
figure;
[cs, hc] = contourf(error_95_3D_array*100);
colorbar
xlim([1,13])
ylim([1,13])
% caxis([0.9, 1.8])
h = colorbar;
ylabel(h, '95% RMS 3D Error [cm]')
set(hc, 'EdgeColor','none')
colormap parula;
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'RMS_3D_95_Error_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);


% Plot contour of results 3D max error.
figure;
[cs, hc] = contourf(error_max_3D_array*100);
colorbar
xlim([1,13])
ylim([1,13])
% caxis([1, 6.0])
caxis([1, 10.0])
h = colorbar;
ylabel(h, 'Max 3D RMS Error [cm]')
set(hc, 'EdgeColor','none')
colormap parula(10);
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'RMS_3D_max_Error_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);


% Plot contour of results ure max error.
figure;
[cs, hc] = contourf(error_max_ure_array*100);
colorbar
xlim([1,13])
ylim([1,13])
caxis([1, 6.0])
h = colorbar;
ylabel(h, 'Max RMS URE [cm]')
set(hc, 'EdgeColor','none')
colormap parula(10);
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'RMS_ure_max_Error_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);


% Plot contour of results fit type.
colors = [hsv(11)];
figure;
hold all;
% [cs, hc] = contourf(fit_type_array);
for ecc_idx = 1:length(ecc_cases)
    for inc_idx = 1:length(inc_cases)
        type = fit_type_array(inc_idx, ecc_idx)-1;
        
        if ~isnan(type)
            plot(ecc_idx, inc_idx, 'ks' ,...
                'MarkerFaceColor', colors(type,:),...
                'MarkerSize', 20, 'linewidth', 2)
        end
    end
end
% shading flat
colorbar
xlim([1,13])
ylim([1,13])
caxis([1.5, 12.5])
h = colorbar;
ylabel(h, 'Fit Type')
set(h,'YTick',[2:12])
set(hc, 'EdgeColor','none')
colormap(colors)
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])
grid on

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'fit_type_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);


% Plot contour of results percentage.
figure;
[cs, hc] = contourf(fit_type_p_array * 100);
colorbar
xlim([1,13])
ylim([1,13])
caxis([0, 100])
h = colorbar;
ylabel(h, '% Alternate Fit Cases')
set(hc, 'EdgeColor','none')
colormap jet(10);
xlabel('Eccentricity')
ylabel('Inclination [deg]')
set(gca,'TickDir','out')
xticks([1:13])
yticks([1:13])

% Eccentricity.
xticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'})

% Inclination.
yticklabels({'0', '10^{-16}', '10^{-12}', '10^{-8}', '10^{-4}', ...
    '10^{-3}', '10^{-2}', '0.1','0.5', '1', '2', '5', '10'})

fileSave = [file_dir_save ,'fit_percentage_All.tiff'];
exportfig(gcf,fileSave,'height',9,'width',12,'fontsize',22,'resolution',220);
