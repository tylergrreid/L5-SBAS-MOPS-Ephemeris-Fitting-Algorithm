function [a, ecc, inc, RAAN, omega, M0,...
    Cus, Cuc, Crs, Crc, Cis, Cic, ...
    IDOT, OMEGA_DOT, delta_n, flag,NumIter] = ...
    COE15_estimator(time, pos, initial_guess, Wmat, ConvCrit, ...
    fit_parameters)
%% DESCRIPTION:
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Date:                 April 6, 2016
%       Updated:              April 6, 2016
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% This function fits the GPS / L5 SBAS MOPS ephemeris message to 
% precision orbit data.
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
%% GLOBAL VARIABLES

global mu R_e omega_e

%% IMPLEMENTATION

% Length of data set.
m = length(time);

% Rearrange the data - stack position vectors on top of each other.
count = 1;
for i = 1:m
    X_data_rearrange(count:count+2) = pos(i,:)';
    count = count + 3;
end

% Extract initial guess info.
a     = initial_guess(1);    % [m]
ecc   = initial_guess(2);    % [-]
inc   = initial_guess(3);    % [rad]
RAAN  = initial_guess(4);    % [rad]
omega = initial_guess(5);    % [rad]
M0    = initial_guess(6);    % [rad]

Cus   = initial_guess(7);    % [rad]
Cuc   = initial_guess(8);    % [rad]
Crs   = initial_guess(9);    % [m]
Crc   = initial_guess(10);   % [m]
Cis   = initial_guess(11);   % [rad]
Cic   = initial_guess(12);   % [rad]

IDOT      = initial_guess(13); % [rad/s]
OMEGA_DOT = initial_guess(14); % [rad/s]
delta_n   = initial_guess(15); % [rad/s]

% Scale to normalize lengths.
scale_meters = a;

% Non-dimensionalize a for numerical stability.
a = a / scale_meters;
Crs = Crs / scale_meters;
Crc = Crc / scale_meters;

% Form intial guess state vector.
COEvec = [a, ecc, inc, RAAN, omega, M0, ...
            Cus, Cuc, Crs, Crc, Cis, Cic, ...
            IDOT, OMEGA_DOT, delta_n];
% COEvec = initial_guess;

% Logic variables.
true  = 1;
false = 0;
done  = false;

% Keep track of the number of iterations.
NumIter = 0;

% Set the failure flag to null. This will ultimately tell us whether or not
% we have failed in producing a message.
flag = 0;

% Form weighting matrix based on diagonal input weights.
count = 1;
W = zeros(length(Wmat)*3);
for i = 1:length(Wmat)
    W(count:count+2,count:count+2) = eye(3)*Wmat(i);
    count = count + 3;
end

% set the maximum number of iterations
MaxIter = 100;

% Define the error smoothing parameters.
forward_diff_coeff = [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6];
D_hat = zeros( 3 * m );
num_neighbours = 6;
D_hat = D_hat + diag( ones( 3 * m, 1 ) ) * forward_diff_coeff(1) ;
for i = 1:num_neighbours
    D_hat = D_hat + diag( ones( 3 * (m - i), 1 ), 3 * i ) * ...
        forward_diff_coeff(i+1);
end
D_hat = D_hat(1:end-3*num_neighbours, 1:end);

u_hat = 0;

% START ITERATION SCHEME FOR 15 ELEMENT ESTIMATE.
while done == false && NumIter <= MaxIter
    % Initialize count variable.
    count = 1;

    for i = 1:length(time)
        % Make an ephemeris structure to pass to the eph2xyz function.
        eph.Asqrt  = sqrt(a * scale_meters);
        eph.e      = ecc;
        eph.i0     = inc;
        eph.Omega0 = RAAN; % TODO this may need to be adjusted???
        eph.Omega  = omega;
        eph.M0     = M0;

        eph.Cus    = Cus;
        eph.Cuc    = Cuc;
        eph.Crs    = Crs * scale_meters;
        eph.Crc    = Crc * scale_meters;
        eph.Cis    = Cis;
        eph.Cic    = Cic;

        eph.IDOT      = IDOT;
        eph.Omega_dot = OMEGA_DOT;
        eph.Delta_n   = delta_n;
        
        eph.Toe = time(1);
        
        % Define the time of transmission.
        ttx = time(i);
        
        % Compute the ECEF position as well as the other needed parameters.
        [X_ECEF, t_k, A_k, e_n, E_k, u_k, i_k, OMEGA_k, PHI_k, r_k] = ...
            eph2xyz(eph, ttx);
        
        A_k = A_k / scale_meters;
        
        
        % Assign the theoretical position vector.
        X_theo(count:count+2) = X_ECEF' / scale_meters;
        
        % Compute vectors for partial derivatives.
        r_k_hat = X_ECEF' / r_k;
        r_k = r_k / scale_meters;
        
        dr_k_hat_du_k = [
            -sin(u_k) * cos(OMEGA_k) - cos(u_k) * cos(i_k) * sin(OMEGA_k);
            -sin(u_k) * sin(OMEGA_k) + cos(u_k) * cos(i_k) * cos(OMEGA_k);
             cos(u_k) * sin(i_k)];
         
        dr_k_hat_du_k2 = [
            -sin(u_k) * sin(OMEGA_k) - cos(u_k) * cos(i_k) * sin(OMEGA_k);
            -sin(u_k) * sin(OMEGA_k) + cos(u_k) * cos(i_k) * cos(OMEGA_k);
             cos(u_k) * sin(i_k)];
        dr_k_hat_du_k2 = dr_k_hat_du_k; 
            
        dr_k_hat_di_k = sin(u_k) * [
             sin(i_k) * sin(OMEGA_k);
            -sin(i_k) * cos(OMEGA_k);
             cos(i_k)];
            
        
        % Compute derivatives.
        drda = (1 - e_n * cos(E_k)) * r_k_hat;
        
        drde = A_k * (e_n - cos(E_k)) / (1 - e_n * cos(E_k)) * r_k_hat +...
            r_k * (2 * sin(E_k) - e_n * sin(E_k) * cos(E_k) - e_n ^ 2 * sin(E_k)) / ...
            sqrt(1 - e_n ^ 2) / (1 - e_n * cos(E_k)) ^ 2 * dr_k_hat_du_k;
        
        drdi0 = r_k * dr_k_hat_di_k;
        
        drdOMEGA0 = r_k * [
            -cos(u_k) * sin(OMEGA_k) - sin(u_k) * cos(i_k) * cos(OMEGA_k); 
             cos(u_k) * cos(OMEGA_k) - sin(u_k) * cos(i_k) * sin(OMEGA_k);
             0];
         
         drdomega = r_k * dr_k_hat_du_k2;
         
         drdM0 = A_k * e_n * sin(E_k) / (1 - e_n * cos(E_k)) * r_k_hat + ...
             r_k * sqrt(1 - e_n^2) / (1 - e_n * cos(E_k)) ^2 * dr_k_hat_du_k;
         
         drdCus = sin(2*PHI_k) * r_k * dr_k_hat_du_k2;
        
         drdCuc = cos(2*PHI_k) * r_k * dr_k_hat_du_k2;
        
         drdCrs = sin(2*PHI_k) * r_k_hat;
         
         drdCrc = cos(2*PHI_k) * r_k_hat;
         
         drdCis = sin(2*PHI_k) * r_k * dr_k_hat_di_k;
         
         drdCic = cos(2*PHI_k) * r_k * dr_k_hat_di_k;
         
         drdIDOT = t_k * drdi0;
         
         drdOMEGA_DOT = t_k * drdOMEGA0;
         
         drddelta_n = t_k * drdM0;
         
        % Form Jacobian matrix.
        A(count:count+2,:) = [
            drda, drde, drdi0, drdOMEGA0, drdomega, drdM0, ...
            drdCus, drdCuc, drdCrs, drdCrc, drdCis, drdCic, ...
            drdIDOT, drdOMEGA_DOT, drddelta_n].*...
            [fit_parameters;fit_parameters;fit_parameters];            
        
        % Update the counter.
        count = count + 3;  
    end
    
    % Solve for the update using Matlab's matrix divide. 
    y = X_data_rearrange' / scale_meters - X_theo';
    dCOEvec = A \ y;

    % If we want weighted least squares.
    % dCOEvec = (W * A) \ ( W * y);
    
    % Update orbital element vector.
    COEvec = COEvec + dCOEvec';
    
    % Assign the new values of the orbital elements.
    a     = COEvec(1);    % [m]
    ecc   = COEvec(2);    % [-]
    inc   = COEvec(3);    % [rad]
    RAAN  = COEvec(4);    % [rad]
    omega = COEvec(5);    % [rad]
    M0    = COEvec(6);    % [rad]
    
    Cus   = COEvec(7);    % [rad]
    Cuc   = COEvec(8);    % [rad]
    Crs   = COEvec(9);    % [m]
    Crc   = COEvec(10);   % [m]
    Cis   = COEvec(11);   % [rad]
    Cic   = COEvec(12);   % [rad]
    
    IDOT      = COEvec(13); % [rad/s]
    OMEGA_DOT = COEvec(14); % [rad/s]
    delta_n   = COEvec(15); % [rad/s]
    
    % Mitigate negative eccetricity.
    if ecc < 0
        ecc = abs( ecc );
        COEvec(2) = ecc;
        M0 = M0 + pi;
        COEvec(6) = M0;
        omega = omega + pi;
        COEvec(5) = omega;
    end
    
    % If we get nonsensical things, we have failed. 
    if a < 0 || a > 1.5 || ecc > 1
        done = true;
        NumIter = MaxIter;
        flag = 1;
    end
    
    % Update number of iterations.
    NumIter = NumIter + 1;
    
    % Check for convergence.
    newton_decrement = norm( A * dCOEvec );
    if newton_decrement < ConvCrit
        done = true;
    end 
end

%% DETERMINE PERFORMANCE AND OUTPUT

% Re-dimensionalize.
a = a * scale_meters; % [m]
Crs = Crs * scale_meters; % [m]
Crc = Crc * scale_meters; % [m]

% Put the angular quantities in the correct range (between -pi and pi)
inc = wrapToPi(inc);
RAAN = wrapToPi(RAAN);
omega = wrapToPi(omega);
M0 = wrapToPi(M0);

% Final check for failure.
if NumIter >= MaxIter
    % Set flag.
    flag = 1;
    
    % Set output to NaN.
    a     = NaN;
    ecc   = NaN;
    inc   = NaN;
    RAAN  = NaN;
    omega = NaN;
    M0    = NaN;
    Cus   = NaN;
    Cuc   = NaN;
    Crs   = NaN;
    Crc   = NaN;
    Cis   = NaN;
    Cic   = NaN;
    IDOT  = NaN;
    OMEGA_DOT = NaN;
    delta_n   = NaN;
    NumIter   = NaN;
end
