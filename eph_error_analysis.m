function [error_3d, rms_error, rms_ure, ...
    error_radial, error_along_track, error_cross_track, ...
    eph] = eph_error_analysis(sqrt_a, ecc, inc, RAAN, omega, M0, ...
    Cus, Cuc, Crc, Crs, Cic, Cis, ...
    IDOT, OMEGA_DOT, delta_n, time, pos, vel, theta_g, ...
    coeff_R2, coeff_A2, coeff_C2)
%% DESCRIPTION
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Start date:           April 17, 2016 
%       Last modified:        April 17, 2016

% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
%  Compute the ephemeris message error summary. 
%
% -------------------------------------------------------------------------
% INPUT:
% 
% - Ephemeris message parameters in meters / radians: 
%     sqrt_a, ecc, inc, RAAN, omega, M0, ...
%         Cus, Cuc, Crc, Crs, Cic, Cis, ...
%         IDOT, OMEGA_DOT, delta_n
% - time, pos, vel: time, position, and velocity of truth data orbit data.
% - theta_g: sidereal time and time(1).
% - coeff_R2, coeff_A2, coeff_C2: coefficients of the analytical URE
%   equation.
%
% -------------------------------------------------------------------------
% OUTPUT:
%
% - error_3d:          Vec of 3D errors.
% - rms_error:         RMS 3D error.
% - rms_ure:           RMS URE error.
% - error_radial:      Vector of radial errors.
% - error_along_track: Vector of along-track errors.
% - error_cross_track: Vector of cross-track errors.
% - eph:               Data structure of the ephemeris message.
%
%% IMPLEMENTATION

% global omega_e
omega_e = 0; % Set to zero in this context because we are working in ECI 
             % instead of ECEF coordinates. To work in ECEF, add omega_e to
             % the list of globals.

% Make an ephemeris structure to pass to the eph2xyz function.
eph.Asqrt  = sqrt_a;
eph.e      = ecc;
eph.i0     = inc;
eph.Omega0 = RAAN;
eph.Omega  = omega;
eph.M0     = M0;

eph.Cus    = Cus;
eph.Cuc    = Cuc;
eph.Crs    = Crs;
eph.Crc    = Crc;
eph.Cis    = Cis;
eph.Cic    = Cic;

eph.IDOT      = IDOT;
eph.Omega_dot = OMEGA_DOT;
eph.Delta_n   = delta_n;

eph.Toe = time(1);

% Number of data points.
n_data = length(time);

% Initialize vectors.
pos_theo = NaN(n_data, 3);
error_3d = NaN(n_data, 1);
error_radial = NaN(n_data, 1);
error_along_track = NaN(n_data, 1);
error_cross_track = NaN(n_data, 1);

% Compute message error.
for i = 1:n_data
    
    % Transmission time.
    ttx = time(i);
    
    % Compute message theorectical position.
    [x_ecef_eph, ~, ~, ~, ~, ~, ~, ~, ~, ~] ...
        = eph2xyz(eph, ttx);
    pos_theo(i,:) = x_ecef_eph;
    
    % Get the sidereal time.
    theta_g_epoch = theta_g + omega_e * time(i);
    
    % Convert ECEF pos / vel to ECI
    [x_eci, v_eci] = ECEF2ECI( pos(i,:)', vel(i,:)', theta_g_epoch );
    [x_eci_eph, ~] = ECEF2ECI( x_ecef_eph', NaN(3,1), theta_g_epoch );
    
    % Compute transformation matrix ECI_2_RIC (radial / along track / 
    % cross track).
    ECI_2_RIC_Mat = ECI2RIC( x_eci, v_eci );
    x_ric = ECI_2_RIC_Mat * x_eci;
    x_ric_eph = ECI_2_RIC_Mat * x_eci_eph;
    
    % Compute the radial / along-track / cross-track errors.
    error_radial(i) = x_ric_eph(1) - x_ric(1);
    error_along_track(i) = x_ric_eph(2) - x_ric(2);
    error_cross_track(i) = x_ric_eph(3) - x_ric(3);
    
    % Compute 3D error.
    error_3d(i) = norm(x_ecef_eph - pos(i,:));
    
end

% RMS 3D error.
rms_error = ( sum( ...
    error_radial.^2 + error_along_track.^2 + error_cross_track.^2 ...
    ) / n_data ) ^ (1/2);

% Calculate the RMS error component errors.
RMS_radial = rms( error_radial );
RMS_cross_track = rms( error_cross_track );
RMS_along_track = rms( error_along_track );

% Calculate the RMS orbit only ure.
rms_ure = sqrt( coeff_A2 * RMS_along_track ^ 2 + ...
    coeff_C2 * RMS_cross_track ^ 2 + ...
    coeff_R2 * RMS_radial ^ 2 );
