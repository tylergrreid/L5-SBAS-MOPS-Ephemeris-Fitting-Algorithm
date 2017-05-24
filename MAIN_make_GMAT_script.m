%% MAKE GMAT SCRIPT TO GET ALL DESIRED ORBITS
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Project Title:        L5 SBAS MOPS
%       Start Date:           April 28, 2017
%       Last Modified:        May 24, 2017
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% This produces a script that configures and runs the NASA General Mission
% Analysis Tool (GMAT) for a variety of orbits and time steps. 
% For more info on GMAT, please see: https://gmat.gsfc.nasa.gov/
%
%% SET UP WORKSPACE

clear
clc
close all

%% SET UP ORBITS / PHYSICAL PARAMETERS

% Bring physical constants into the workspace. 
physical_constants_GPS;

% Altitudes.
altitudes =  42164-R_e/1000; % [km]

% Create Inclinations and eccentricities to test. 
ecc = [0, 1e-16, 1e-12, 1e-8, 1e-4, 1e-3, ...
    1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6]';
inc = [0, 1e-16, 1e-12, 1e-8, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1, 2, 5, 10]';

% Define RAAN. 
RAAN = 270;
AOP = 0;
TA = 0;

% The number of orbits being examined.
num_cases = length(ecc) * length(inc) * ...
    length(RAAN) * length(AOP) * length(TA);
disp(['Number of cases: ',num2str(num_cases)])

% Set MEOs/IGSOs to be have typical GPS satellite physical parameters.
area_MEOplus = 21; % [m^2]
mass_MEOplus = 1665; % [kg]
inc_MEOplus = 55; % [deg]

% Start and end date of the simulation. 
start_date = '01 Jan 2017 00:00:00.000';
end_date   = '01 Jan 2018 00:00:00.000';

% Propagation time.
prop_time = 86400; % [sec]

% Time between successive propagations.
next_prop = 366 * 24 * 60 * 60; % [sec]

% Start times for propagations.
times = 0:next_prop:(datenum(end_date) - datenum(start_date))*24*3600;

disp(['Number of orbit props: ',num2str(length(times))])

disp(['Number of total cases: ',num2str(num_cases * length(times))])

%% MAKE SPACECRAFT NAMES

max_digits = length(num2str(length(times) * num_cases)); % The maximum number of digits.
preamble = '';
for i = 1:max_digits-1
    preamble = [preamble,'0'];
end

for alt = 1:length(altitudes)
    for orb_param = 1:num_cases
        idx = length(num2str(orb_param))-1;
        
        sv_name{alt,orb_param} = ['MOPSSat_alt_',...
            num2str(floor(altitudes(alt))),'_km_RAAN_',num2str(RAAN),'_',preamble(1:end-idx),num2str(orb_param)];
    end
end

%% MAKE THE ORBIT EVERY X MINUTES FOR A YEAR (MAKE 24 HOUR TRACKS)

% Open text file.
fileID = fopen(horzcat('RUN_ALL_GMAT_ORBITS_RAAN_',num2str(RAAN),'.script'),'w');

% initialize case number.
case_number = 1;

% Create all spacecraft orbits.
for alt = 1:length(altitudes)
    for ecc_idx = 1:length(ecc)
        for inc_idx = 1:length(inc)
            for RAAN_idx = 1:length(RAAN)
                % Create spacecraft.
                fprintf(fileID,['Create Spacecraft ',sv_name{alt,case_number},';\n']);
                
                % Update date / time.
                new_date = addtodate(datenum(start_date), times(1), 'second');
                
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.DateFormat = UTCGregorian;\n']);
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.Epoch = ''',datestr(new_date,'dd mmm yyyy HH:MM:SS.FFF'),''';\n']);
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.CoordinateSystem = EarthMJ2000Eq;\n']);
                
                % Update orbital elements.
                semi_major_axis = R_e/1000 + altitudes(alt); % [km]
                %n = sqrt(mu/semi_major_axis^3); % [rad/s] - Mean orbital rate.
                %true_anomaly = TA + n*times(t)*180/pi; % [deg]
                %true_anomaly = mod(true_anomaly, 360); % put between 0 and 360.
                
                % Choose true anomaly from uniform random
                % true_anomaly = unifrnd(0, 360); % [deg]
                true_anomaly = TA; 
                
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.SMA = ',num2str(semi_major_axis),';\n']); % [km]
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.ECC = ',num2str(ecc(ecc_idx)),';\n']); % [-]
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.RAAN = ',num2str(RAAN(RAAN_idx)),';\n']); % [deg]
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.AOP = ',num2str(AOP),';\n']); % [deg]
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.TA = ',num2str(true_anomaly),';\n']); % [deg]
                
                % Inclination
                if altitudes(alt)>20000 % MEO and above.
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.INC = ',num2str(inc(inc_idx)),';\n']); % [deg]
                else % LEO
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.INC = ',num2str(inc_LEO),';\n']); % [deg]
                end % end if.
                
                % Physical parameters.
                if altitudes(alt)>20000 % MEO and above
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.DryMass = ',num2str(mass_MEOplus),';\n']); % [kg]
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.SRPArea = ',num2str(area_MEOplus),';\n']); % [m^2]
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.DragArea = ',num2str(area_MEOplus),';\n']); % [m^2]
                else % LEO
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.DryMass = ',num2str(mass_LEO),';\n']); % [kg]
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.SRPArea = ',num2str(area_LEO),';\n']); % [m^2]
                    fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.DragArea = ',num2str(area_LEO),';\n']); % [m^2]
                end % end if.
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.Cr = 1.8;\n']);
                fprintf(fileID,['GMAT ',sv_name{alt,case_number},'.Cd = 2.2;\n']);
                
                % Save the eccentricity, inclination, as a function of case
                % number. 
                ecc_orb_scenario(case_number) = ecc(ecc_idx);
                inc_orb_scenario(case_number) = inc(inc_idx);
                
                % March forward. 
                case_number = case_number + 1;
                
                % New line.
                fprintf(fileID,'\n');
            end % end RAAN.
        end % end inc.
    end % end ecc.
end % end alt.

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
for alt = 1:length(altitudes)
    for orb_param = 1:num_cases
        % Create ephemeris file.
        fprintf(fileID,...
            ['Create EphemerisFile EphmerisFile_',sv_name{alt,orb_param},';\n']);
        
        % Set spacecraft.
        fprintf(fileID,...
            ['EphmerisFile_',sv_name{alt,orb_param},'.Spacecraft = ',sv_name{alt,orb_param},';\n']);
        
        % Set file name output.
        fprintf(fileID,...
            ['EphmerisFile_',sv_name{alt,orb_param},'.Filename = ''','EphemerisFile_',sv_name{alt,orb_param},'.eph'';\n']);
        
        % Set coordinate system, we'll used ECEF.
%         fprintf(fileID,...
%             ['EphmerisFile_',sv_name{alt,orb_param},'.CoordinateSystem = EarthFixed;\n']);
    
        % This is the code for inertial coordinates.
        fprintf(fileID,...
            ['EphmerisFile_',sv_name{alt,orb_param},'.CoordinateSystem = EarthMJ2000Eq;\n']);
        
        % Set the step size.
        fprintf(fileID,...
            ['EphmerisFile_',sv_name{alt,orb_param},'.StepSize = 10;\n']);
        
        % New line.
        fprintf(fileID,'\n');
    end % end orb_param.
end % end alt.

%% RUN THE MISSION

fprintf(fileID,'BeginMissionSequence;\n');

% Run all orbits.
for alt = 1:length(altitudes)
    for orb_param = 1:num_cases
        % Propagate for 24 hours.
        fprintf(fileID,...
            ['Propagate prop(',sv_name{alt,orb_param},') {',sv_name{alt,orb_param},'.ElapsedSecs = 86400.0};\n']);
    end % end t.
end % end alt.

% Close file.
fclose(fileID);

%% SAVE NUMBERING AND PARAMETERS

ecc_cases = ecc;
inc_cases = inc;

save('orbit_numbering_GMAT.mat', 'sv_name' , 'ecc_cases', 'inc_cases',...
    'ecc_orb_scenario', 'inc_orb_scenario', ...
    'AOP', 'TA', 'altitudes', 'num_cases') 
