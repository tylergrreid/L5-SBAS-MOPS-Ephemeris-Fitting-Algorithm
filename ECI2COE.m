function [coe, undefined, orbit_type] = ECI2COE(X_ECI,V_ECI)
%% DESCRIPTION:
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       AA 279:               Problem Set 4
%       Date:                 April 22, 2011
%       Date modified:        May 3, 2016
%
%       Modified to handle special cases of circular and equitorial.
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Based on Vallado (2007) Algorithm 10
%
% This algorithm will compute the Earth Centered Inertial (ECI) position
% and velocity vectors which are equivalent to the given classical orbital
% elements
%
% -------------------------------------------------------------------------
% INPUT:
%
%       X = ECI position vector of the spacecraft    [length] (1)
%       V = ECI veloicty vector of the spacecraft    [length / time] (1)
%
% -------------------------------------------------------------------------
% OUTPUT:
%
%       a = semi-major axis                          [length] (1)
%       e = eccentrity                               [deg]
%       i = inclination                              [deg]
%   OMEGA = right ascension of the ascending node    [deg]
%       w = argument of perigee                      [deg]
%       f = true anomaly                             [deg]
%
% -------------------------------------------------------------------------
% NOTES:
%
% (1) This quantity can be expressed in either m or km or etc as long
%     as the global value of mu (the Earth's gravitational parameter) is in
%     consitant units.
%
% (2) This algorithm can handle special cases such as circular inclined,
%     elliptical equitorial, circular equitorial and outputs the
%     corresponding relevant parameters.
%
%% DEFINE GLOBAL VARIABLES USED

global mu

%% IMPLEMTENTATION

% Compute the magnitude of position and velocity.
r = norm(X_ECI);
v = norm(V_ECI);

% Compute angular momentum vector.
h_vec = cross(X_ECI,V_ECI);
h     = norm(h_vec);

% Compute the line of nodes vector.
n_vec = cross([0 0 1]',h_vec);
n     = norm(n_vec);

% Compute the specific mechanical energy.
energy = (v^2)/2 - mu/r;

% Compute the eccentricity vector.
e_vec = ((v^2 - mu/r)*X_ECI - dot(X_ECI,V_ECI)*V_ECI ) / mu;


% -------------------------------------------------------------------------
% eccentricity e
% -------------------------------------------------------------------------
e = norm(e_vec);

tol = 1e-15;

if abs(e-1)<tol
    e = 1;
end

if abs(e)<tol
    e = 0;
end

% -------------------------------------------------------------------------
% semi parameter p
% -------------------------------------------------------------------------

% Filter parabolic case.
if e ~= 1
    a = -mu/2/energy;
    p = a*(1-e^2);
else % parabolic case
    a = Inf;
    p = h^2/mu;
end


% -------------------------------------------------------------------------
% inclination i
% -------------------------------------------------------------------------

i = acos(h_vec(3)/h); % [rad]

% -------------------------------------------------------------------------
% RAAN
% -------------------------------------------------------------------------

% Compute the right ascention of the ascending node (RAAN).
RAAN = acos( n_vec(1)/n );

% Determine the quadrant.
if n_vec(2)<0
    RAAN = 2*pi-RAAN;
end


% -------------------------------------------------------------------------
% argument of perigee w
% -------------------------------------------------------------------------

% Compute the argument of perigee (w).
w = real(acos( dot(n_vec,e_vec) / e / n )); % [rad]

% Determine the quadrant.
if e_vec(3)<0
    w = 2*pi-w;
end


% -------------------------------------------------------------------------
% true anomaly f
% -------------------------------------------------------------------------

% True anomaly (f).
f = real(acos( dot(e_vec,X_ECI) / r / e )); % [rad]

% Determine the quadrant.
if dot(X_ECI,V_ECI)<0
    f = 2*pi-f;
end

% -------------------------------------------------------------------------
% mean anomaly f
% -------------------------------------------------------------------------

% Find the eccentric anomaly.
cosf = cos(f);
cosE = (e+cosf)/(1+e*cosf);
sinE = sin(f)*sqrt(1-e^2)/(1+e*cosf);

E = atan2(sinE,cosE);

% Compute mean anomaly.
M = E - e*sin(E);

%% HANDLE SPECIAL CASES AND FORMAT OUTPUT

% Set special case = 0, if there are no special cases, output the above 6
% COE's

% Determine which special case we have.
orbit_case = 1; % Default.

% Case I - Elliptical, Equitorial.
if (e<1 && e~=0) && (i == 0 || i == pi)
    orbit_case = 2;
end

% Case II - Circular, Inclined. 
if e == 0 && (i ~= 0 || i ~= pi)
    orbit_case = 3;
end

% Case III - Circular, Equitorial. 
if e == 0 && (i == 0 || i == pi)
    orbit_case = 4;
end

switch orbit_case
    
    case 1 % Default.
        
        coe.p = p;
        coe.a = a;
        coe.e = e;
        coe.i = i*180/pi;
        coe.omega = w*180/pi;
        coe.RAAN = RAAN*180/pi;
        coe.f = f*180/pi;
        coe.M = M*180/pi;
        
        undefined.p = 0;
        undefined.a = 0;
        undefined.e = 0;
        undefined.i = 0;
        undefined.omega = 0;
        undefined.RAAN = 0;
        undefined.f = 0;
        undefined.M = 0;
        
        orbit_type = 'elliptical inclined';
        
    case 2 % Elliptical equitorial.
        
        % Compute the true longitude of periapsis.
        w_true = acos(e_vec(1)/e); % [rad]
        
        % Determine the quadrant.
        if e_vec(2)<0
            w_true =  2*pi-w_true;
        end
        
        % Output structure.
        coe.p = p;
        coe.a = a;
        coe.e = e;
        coe.i = i*180/pi;
        coe.omega = w_true*180/pi;
        coe.RAAN = 0;
        coe.f = f*180/pi;
        coe.M = M*180/pi;
        
        undefined.p = 0;
        undefined.a = 0;
        undefined.e = 0;
        undefined.i = 0;
        undefined.omega = 1;
        undefined.RAAN = 1;
        undefined.f = 0;
        undefined.M = 0;
        
        orbit_type = 'elliptical equitorial';
        
    case 3 % Circular Inclined.
        
        % Compute the argument of latitude.
        u = acos( dot(n_vec,X_ECI) / n / r );
        
        % Determine the quadrant.
        if X_ECI(3)<0
            u = 2*pi-u;
        end
        
        % Output structure.
        coe.p = p;
        coe.a = a;
        coe.e = e;
        coe.i = i*180/pi;
        coe.omega = 0;
        coe.RAAN = RAAN*180/pi;
        coe.f = u*180/pi;
        coe.M = coe.f;
        
        undefined.p = 0;
        undefined.a = 0;
        undefined.e = 0;
        undefined.i = 0;
        undefined.omega = 1;
        undefined.RAAN = 0;
        undefined.f = 1;
        undefined.M = 1;
        
        orbit_type = 'circular inclined';
        
    case 4 % Circular Equitorial.
        
        % Compute the true longitude.
        lambda_true = acos(X_ECI(1) / r);
        
        % Determine the quadrant.
        if X_ECI(2)<0
            lambda_true = 2*pi-lambda_true;
        end
        
        % Output structure.
        coe.p = p;
        coe.a = a;
        coe.e = e;
        coe.i = i*180/pi;
        coe.omega = 0;
        coe.RAAN = 0;
        coe.f = lambda_true*180/pi;
        coe.M = coe.f;
        
        undefined.p = 0;
        undefined.a = 0;
        undefined.e = 0;
        undefined.i = 0;
        undefined.omega = 1;
        undefined.RAAN = 1;
        undefined.f = 1;
        undefined.M = 1;
        
        orbit_type = 'circular equitorial';
end



