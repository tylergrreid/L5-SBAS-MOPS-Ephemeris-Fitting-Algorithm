function [a, ecc, inc, RAAN, omega, M0,...
    Cus, Cuc, Crs, Crc, Cis, Cic, ...
    IDOT, OMEGA_DOT, delta_n, flag, NumIter, fit_type] = ...
    COE15_estimator_wrapper(time, pos, vel,...
    initial_guess, Wmat, ConvCrit, ...
    fit_parameters, theta_g, coeff_R2, coeff_A2, coeff_C2)
%% DESCRIPTION
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Date:                 May   , 2017
%       Updated:              May 24, 2017
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% This function fits the GPS / L5 SBAS MOPS ephemeris message to 
% precision orbit data.
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% This function fits the GPS / L5 SBAS MOPS ephemeris message to 
% precision orbit data. It attempts combinations of removals of certain
% orbital elements in the event of convergence failure and chooses the best
% option. 
%
% -------------------------------------------------------------------------
% INPUT:
%       time          = Input vector of times starting from 0 [sec] which
%                       are consistent with the pos / vel (assumes a fixed
%                       time step).
%       pos           = Position of the spacecraft in an ECEF
%                       coordinate frame.
%       vel           = Velocity of the spacecraft in an ECEF
%                       coordinate frame.
%       initial_guess = Initial guess of orbital elements in the form:
%                       initial_guess = [a; e; inc; RAAN; omega; M0;
%                           Cus; Cuc; Crs; Crc; Cis; Cic;
%                           IDOT; OMEGA_DOT; delta_n];
%       Wmat          = Values to form the weighting matrix for least
%                       squares, there should be one value for each point.
%       ConvCrit      = Value of the convergence criteria used in
%                       evaluating the newton decrement. 
%       fit_parameters= Parameters to be fit of the 15 orbital elements. 
%       theta_g       = Greenwich mean sidereal time [rad], see: utc2gmst.m
%       coeff_R2, coeff_A2, coeff_C2 = 
%                       Analytical coefficients for the URE equations.
%
% -------------------------------------------------------------------------
% OUTPUT:
%       a            = Best fit semi-major axis
%       e            = Best fit eccentricity
%       inc          = Best fit inclination
%       RAAN         = Best fit right ascension of the ascending node
%       omega        = Best fit argument of perigee
%       M0           = Best fit mean anomaly at epoch
%       Cus,Cuc      = Best fit along-track harmonic correction terms
%       Crs,Crc      = Best fit radial harmonic correction terms
%       Cis,Cic      = Best fit cross-track harmonic correction terms
%       IDOT         = Best fit inclination correction rate
%       OMEGA_DOT    = Best fit rate in right ascention
%       delta_n      = Best fit orbital rate offset
%       flag         = 0 if convergence fails, 1 if successfull
%       NumIter      = Number of Iterations
%
%% MAIN

% Define the fit parameter scenarios to be tested.
fit_scenarios = repmat(fit_parameters, 10, 1);

% No RAAN.
fit_scenarios(2, 4) = 0;

% No omega.
fit_scenarios(3, 5) = 0;

% No inclination.
fit_scenarios(4, 3) = 0;

% No RAAN / inc.
fit_scenarios(5, 3) = 0;
fit_scenarios(5, 4) = 0;

% No omega / inc.
fit_scenarios(6, 3) = 0;
fit_scenarios(6, 5) = 0;

% No RAAN / ecc.
fit_scenarios(7, 2) = 0;
fit_scenarios(7, 4) = 0;

% No omega / ecc.
fit_scenarios(8, 2) = 0;
fit_scenarios(8, 4) = 0;

% No RAAN / omega.
fit_scenarios(9, 4) = 0;
fit_scenarios(9, 5) = 0;

% No RAAN / omega / inc.
fit_scenarios(10, 3) = 0;
fit_scenarios(10, 4) = 0;
fit_scenarios(10, 5) = 0;

% No RAAN / omega / ecc.
fit_scenarios(11, 2) = 0;
fit_scenarios(11, 4) = 0;
fit_scenarios(11, 5) = 0;

% No RAAN / omega / ecc / inc.
fit_scenarios(12, 2) = 0;
fit_scenarios(12, 3) = 0;
fit_scenarios(12, 4) = 0;
fit_scenarios(12, 5) = 0;

% Get the number of scenarios.
[num_scenarios,~] = size(fit_scenarios);

% Initialize the rms ure (this is a dummy value).
rms_ure_best = 999;

for i = 1:num_scenarios
    % Try fitting without removing elements.
    [a_test, ecc_test, inc_test, RAAN_test, omega_test, M0_test,...
        Cus_test, Cuc_test, Crs_test, Crc_test, Cis_test, Cic_test, ...
        IDOT_test, OMEGA_DOT_test, delta_n_test, ...
        flag_test, NumIter_test] = ...
        COE15_estimator(time, pos, initial_guess, Wmat, ConvCrit, ...
        fit_scenarios(i,:));
    
    % If the inclination is less than zero, we have failed.
    if inc_test < 0
        flag_test = 1;
    end
    
    % Evaluate the message performance.
    if flag_test == 0
        [~, ~, rms_ure, ~, ~, ~, ~] = ...
            eph_error_analysis(sqrt(a_test), ecc_test, inc_test,...
            RAAN_test, omega_test, M0_test, ...
            Cus_test, Cuc_test, Crc_test, Crs_test, ...
            Cic_test, Cis_test, ...
            IDOT_test, OMEGA_DOT_test, delta_n_test, ...
            time, pos, vel, theta_g, ...
            coeff_R2, coeff_A2, coeff_C2);
    end

    if flag_test == 0
        if rms_ure < rms_ure_best            
            % Assign the best rms_ure seen so far.
            rms_ure_best = rms_ure;
            
            % Assign the fit type.
            fit_type = i;
            
            % Assign the orbital elements as the output.
            a = a_test; % [m]
            ecc = ecc_test; % [-]
            inc = inc_test; % [rad]
            RAAN = RAAN_test; % [rad]
            omega = omega_test; % [rad]
            M0 = M0_test; % [rad]
            
            Cus   = Cus_test;    % [rad]
            Cuc   = Cuc_test;    % [rad]
            Crs   = Crs_test;    % [m]
            Crc   = Crc_test;   % [m]
            Cis   = Cis_test;   % [rad]
            Cic   = Cic_test;   % [rad]
            
            IDOT      = IDOT_test; % [rad/s]
            OMEGA_DOT = OMEGA_DOT_test; % [rad/s]
            delta_n   = delta_n_test; % [rad/s]
            
            % Number of iterations.
            NumIter = NumIter_test;
            
            % Set the flag.
            flag = flag_test;
            
            % If we succeed without having to remove any elements, stop. 
            if i == 0
               break; 
            end
        end
    end
end
