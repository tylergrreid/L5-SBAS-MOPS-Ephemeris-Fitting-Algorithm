%% PHYSICAL CONSTANTS
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Start Date:           April 7, 2016
%       Last updated:         April 7, 2016
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% Loading this file enters all of the physical constants needed to perform
% all of the neccessary calculations into the global workspace.  All
% constants are given in SI units. This are the physical constants as
% defined by the GPS interface spec.
%
%% DEFINE CONSTANTS

global mu omega_e R_e 

% Earth's Mean Equitorial Radius
R_e =      6.37813649e6;        % [m]

% Earth Gravitational Parameter mu = G*M_earth
mu =       3.986005e14;      % [m^3/s^2]

% Mean Angular Velocity of the Earth
omega_e =  7.2921151467e-5;    % [rad/s]

% J_2 - second zonal harmonic
% J_2 =      1.0826300e-3;        % [-]

% Earth's shape - eccentricity^2
% Earth_E2 = 0.006694385000;      % [-]


