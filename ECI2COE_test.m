%% ECI2COE TEST SCRIPT
% 
%  This is to test the ECI to Classical Orbital Elements function. 
%
%  Date: May 3, 2017
%
%% WORKSPACE

clc
clear
close all

% load physical constants file to enter them into the global workspace
physical_constants_GPS

% Define global variables.
global omega_e mu R_e

%% TEST 1 - ARBITRARY ORBIT

% Orbital elements. 
a = R_e + 800e3; % [m]
e = 0.1;
inc = 51.0 * pi / 180; % [rad]
RAAN = 192.0 * pi / 180; % [rad]
omega = 19.0 * pi / 180; % [rad]
M = 45.0 * pi / 180; % [deg]

% Convert to ECI position and velocity. 
[X_truth,V_truth] = COE2RV(a,e,inc,RAAN,omega,M);

% Get orbital elemtents. 
[coe, undefined, orbit_type] = ECI2COE(X_truth,V_truth);
% coe

% Convert orbital elements back to ECI position and velcity. 
[X_test,V_test] = COE2RV(coe.a,coe.e,coe.i*pi/180,...
    coe.RAAN*pi/180,coe.omega*pi/180,coe.M*pi/180);

% Assess error. 
position_error = norm(X_truth - X_test);
velocity_error = norm(V_truth - V_test);

disp('Test 1 - Inclined, Elliptical')
disp(['Position error = ',num2str(position_error*100), ' cm'])
disp(['Velocity error = ',num2str(velocity_error*100), ' cm/second'])
coe
orbit_type
disp(' ')

%% TEST 2 - ELLIPTICAL EQUITORIAL

% Orbital elements. 
a = R_e + 800e3; % [m]
e = 0.1;
inc = 0.0 * pi / 180; % [rad]
RAAN = 270.0 * pi / 180; % [rad]
omega = 19.0 * pi / 180; % [rad]
M = 45.0 * pi / 180; % [deg]

% Convert to ECI position and velocity. 
[X_truth,V_truth] = COE2RV(a,e,inc,RAAN,omega,M);

% Get orbital elemtents. 
[coe, undefined, orbit_type] = ECI2COE(X_truth,V_truth);
% coe

% Convert orbital elements back to ECI position and velcity. 
[X_test,V_test] = COE2RV(coe.a,coe.e,coe.i*pi/180,...
    coe.RAAN*pi/180,coe.omega*pi/180,coe.M*pi/180);

% Assess error. 
position_error = norm(X_truth - X_test);
velocity_error = norm(V_truth - V_test);

disp('Test 2 - Elliptical, Equitorial')
disp(['Position error = ',num2str(position_error*100), ' cm'])
disp(['Velocity error = ',num2str(velocity_error*100), ' cm/second'])
coe
orbit_type
disp(' ')

%% TEST 3 - CIRCULAR INCLINED

% Orbital elements. 
a = R_e + 800e3; % [m]
e = 0.0;
inc = 2 * pi / 180; % [rad]
RAAN = 192.0 * pi / 180; % [rad]
omega = 145.0 * pi / 180; % [rad]
M = 45.0 * pi / 180; % [deg]

% Convert to ECI position and velocity. 
[X_truth,V_truth] = COE2RV(a,e,inc,RAAN,omega,M);

% Get orbital elemtents. 
[coe, undefined, orbit_type] = ECI2COE(X_truth,V_truth);

% Convert orbital elements back to ECI position and velcity. 
[X_test,V_test] = COE2RV(coe.a,coe.e,coe.i*pi/180,...
    coe.RAAN*pi/180,coe.omega*pi/180,coe.M*pi/180);

% Assess error. 
position_error = norm(X_truth - X_test);
velocity_error = norm(V_truth - V_test);

disp('Test 3 - Circular Inclined')
disp(['Position error = ',num2str(position_error*100), ' cm'])
disp(['Velocity error = ',num2str(velocity_error*100), ' cm/second'])
coe
orbit_type
disp(' ')

%% TEST 4 - CIRCULAR EQUITORIAL

% Orbital elements. 
a = R_e + 800e3; % [m]
e = 0;
inc = 1e-9; % [rad]
RAAN = 192.0 * pi / 180; % [rad]
omega = 19.0 * pi / 180; % [rad]
M = 45.0 * pi / 180; % [deg]

% Convert to ECI position and velocity. 
[X_truth,V_truth] = COE2RV(a,e,inc,RAAN,omega,M);

% Get orbital elemtents. 
[coe, undefined, orbit_type] = ECI2COE(X_truth,V_truth);

% Convert orbital elements back to ECI position and velcity. 
[X_test,V_test] = COE2RV(coe.a,coe.e,coe.i*pi/180,...
    coe.RAAN*pi/180,coe.omega*pi/180,coe.M*pi/180);

% Assess error. 
position_error = norm(X_truth - X_test);
velocity_error = norm(V_truth - V_test);

disp('Test 4 - Circular, Equitorial')
disp(['Position error = ',num2str(position_error*100), ' cm'])
disp(['Velocity error = ',num2str(velocity_error*100), ' cm/second'])
coe
orbit_type
