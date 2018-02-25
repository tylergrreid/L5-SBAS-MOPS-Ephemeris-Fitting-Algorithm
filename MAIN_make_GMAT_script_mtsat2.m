%% MAKE GMAT SCRIPT TO GET ALL DESIRED ORBITS
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Project Title:        L5 SBAS MOPS
%       Start Date:           April 28, 2017
%       Last Modified:        February 25, 2018
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% This produces a script that configures and runs the NASA General Mission
% Analysis Tool (GMAT) for a variety of orbits and time steps. 
% For more info on GMAT, please see: https://gmat.gsfc.nasa.gov/
%
% This particular version creates the script for MTSAT-2 for evaluation of
% the L5 MOPS ephemeris message. 
%
%% SET UP WORKSPACE

clear
clc
close all

%% SET UP ORBITS / PHYSICAL PARAMETERS

% Bring physical constants into the workspace. 
physical_constants_GPS;

% Define initial conditions of the orbit. 
SMA = 42164.13866; % [km]
ecc = 0.00042061; % [-]
inc = 0.00195; % [deg]
RAAN = 0.391565; % [deg]
AOP = 327.12287; % [deg]
MA = 304.10548; % [deg]

EA = Keplers_Eqn(MA*pi/180, ecc); % [rad]
cosTA = (cos(EA) - ecc) / (1 - ecc * cos(EA));
sinTA = sin(EA) * sqrt(1 - ecc^2) / (1 - ecc * cos(EA)); 

TA = atan2(sinTA, cosTA) * 180 / pi; 
TA = wrapTo360(TA);

% Typical GPS mass / areas. 
area_MEOplus = 21; % [m^2]
mass_MEOplus = 1665; % [kg]

% Start and end date of the simulation. 
start_date = '03 Feb 2018 23:30:00.000';

% Propagation time.
prop_time = 86400 * 7; % [sec]

% Spacecraft name. 
sv_name = 'mtsat2'; 

%% MAKE THE ORBIT EVERY X MINUTES FOR A YEAR (MAKE 24 HOUR TRACKS)

% Open text file.
fileID = fopen(horzcat(pwd,'/Results_Fitting_MTSAT2/RUN_GMAT_mtsat2.script'),'w');

% Create spacecraft.
fprintf(fileID,['Create Spacecraft ',sv_name,';\n']);

% Define coordinates / timing. 
fprintf(fileID,['GMAT ',sv_name,'.DateFormat = UTCGregorian;\n']);
fprintf(fileID,['GMAT ',sv_name,'.Epoch = ''',datestr(start_date,'dd mmm yyyy HH:MM:SS.FFF'),''';\n']);
fprintf(fileID,['GMAT ',sv_name,'.CoordinateSystem = EarthMJ2000Eq;\n']);

% Define orbital elements. 
fprintf(fileID,['GMAT ',sv_name,'.SMA = ',num2str(SMA),';\n']); % [km]
fprintf(fileID,['GMAT ',sv_name,'.ECC = ',num2str(ecc),';\n']); % [-]
fprintf(fileID,['GMAT ',sv_name,'.RAAN = ',num2str(RAAN),';\n']); % [deg]
fprintf(fileID,['GMAT ',sv_name,'.AOP = ',num2str(AOP),';\n']); % [deg]
fprintf(fileID,['GMAT ',sv_name,'.TA = ',num2str(TA),';\n']); % [deg]
fprintf(fileID,['GMAT ',sv_name,'.INC = ',num2str(inc),';\n']); % [deg]

% Define physical parameters. 
fprintf(fileID,['GMAT ',sv_name,'.DryMass = ',num2str(mass_MEOplus),';\n']); % [kg]
fprintf(fileID,['GMAT ',sv_name,'.SRPArea = ',num2str(area_MEOplus),';\n']); % [m^2]
fprintf(fileID,['GMAT ',sv_name,'.DragArea = ',num2str(area_MEOplus),';\n']); % [m^2]

fprintf(fileID,['GMAT ',sv_name,'.Cr = 1.8;\n']);
fprintf(fileID,['GMAT ',sv_name,'.Cd = 2.2;\n']);

%% FORCE MODEL

% Create force model.
fprintf(fileID,'Create ForceModel fm;\n');

% Earth gravity + tides.
fprintf(fileID,'GMAT fm.CentralBody = Earth;\n');
fprintf(fileID,'GMAT fm.PrimaryBodies = {Earth};\n');
fprintf(fileID,'GMAT fm.GravityField.Earth.Degree = 70;\n');
fprintf(fileID,'GMAT fm.GravityField.Earth.Order  = 70;\n');
fprintf(fileID,'GMAT fm.GravityField.Earth.PotentialFile = ''EGM96.cof'';\n'); % GPS currently uses the EGM96 model which has degree / order 70
fprintf(fileID,'GMAT fm.GravityField.Earth.EarthTideModel = ''SolidAndPole'';\n');
fprintf(fileID,'\n');

% Third body gravity.
% GPS currently only models the Moon and Sun
fprintf(fileID,'GMAT fm.PointMasses = {Luna, Sun};\n');
fprintf(fileID,'\n');

% Drag. - This is negligable for GPS but will keep it on for LEO.
fprintf(fileID,'GMAT fm.Drag.AtmosphereModel = MSISE90;\n');
fprintf(fileID,'GMAT fm.Drag.HistoricWeatherSource = ''ConstantFluxAndGeoMag'';\n'); % Default
fprintf(fileID,'GMAT fm.Drag.PredictedWeatherSource = ''ConstantFluxAndGeoMag'';\n'); % Default
fprintf(fileID,'GMAT fm.Drag.F107 = 150;\n'); % Default
fprintf(fileID,'GMAT fm.Drag.F107A = 150;\n'); % Default
fprintf(fileID,'GMAT fm.Drag.MagneticIndex = 3;\n'); % Default
fprintf(fileID,'\n');

% Solar Radiation Pressure.
% GPS has its own specific model we won't get into the details of here.
fprintf(fileID,'GMAT fm.SRP = On;\n');
fprintf(fileID,'GMAT fm.SRP.Flux = 1367;\n'); % Default
fprintf(fileID,'GMAT fm.SRP.SRPModel = Spherical;\n'); % Default
fprintf(fileID,'GMAT fm.SRP.Nominal_Sun = 149597870.691;\n'); % Default
fprintf(fileID,'\n');

% Relativity.
fprintf(fileID,'GMAT fm.RelativisticCorrection = On;\n'); % GPS currently includes this.
fprintf(fileID,'\n');

% Error control.
fprintf(fileID,'GMAT fm.ErrorControl = RSSStep;\n');
fprintf(fileID,'\n');

%% SET UP PROPAGATOR

% Create propagator.
fprintf(fileID,'Create Propagator prop;\n');

fprintf(fileID,'GMAT prop.FM = fm;\n');
fprintf(fileID,'GMAT prop.Type = RungeKutta89;\n'); % Good performance in LEO according to documentation
fprintf(fileID,'GMAT prop.InitialStepSize = 30;\n');
fprintf(fileID,'GMAT prop.Accuracy = 9.999999999999999e-12;\n');
fprintf(fileID,'GMAT prop.MinStep = 0.001;\n');
fprintf(fileID,'GMAT prop.MaxStep = 30;\n');
fprintf(fileID,'GMAT prop.MaxStepAttempts = 50;\n');
fprintf(fileID,'GMAT prop.StopIfAccuracyIsViolated = true;\n');
fprintf(fileID,'\n');

%% SET UP EPHEMERIS FILES

% Create / configure ephemeris files.

% Create ephemeris file.
fprintf(fileID,...
    ['Create EphemerisFile EphmerisFile_',sv_name,';\n']);

% Set spacecraft.
fprintf(fileID,...
    ['EphmerisFile_',sv_name,'.Spacecraft = ',sv_name,';\n']);

% Set file name output.
fprintf(fileID,...
    ['EphmerisFile_',sv_name,'.Filename = ''','EphemerisFile_',sv_name,'.eph'';\n']);

% Set coordinate system, we'll used ECEF.
% fprintf(fileID,...
%    ['EphmerisFile_',sv_name{alt,orb_param},'.CoordinateSystem = EarthFixed;\n']);

% This is the code for inertial coordinates.
fprintf(fileID,...
    ['EphmerisFile_',sv_name,'.CoordinateSystem = EarthMJ2000Eq;\n']);

% Set the step size.
fprintf(fileID,...
    ['EphmerisFile_',sv_name,'.StepSize = 10;\n']);

% New line.
fprintf(fileID,'\n');

%% RUN THE MISSION

fprintf(fileID,'BeginMissionSequence;\n');

% Propagate for 24 hours.
fprintf(fileID,...
    ['Propagate prop(',sv_name,') {',sv_name,'.ElapsedSecs = ',num2str(prop_time),'.0};\n']);

% Close file.
fclose(fileID);
