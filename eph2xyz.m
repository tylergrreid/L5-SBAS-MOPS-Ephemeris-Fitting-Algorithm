function [X_ECEF, t_k, A_k, e_n, E_k, u_k, i_k, OMEGA_k, PHI_k, r_k] ...
    = eph2xyz(eph,ttx)
%% DESCRIPTION
%
%       Written by: Tyler Reid (tyreid@stanford.edu)
%       Date:       Feb, 2010
%       Course:     AA 272C
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION:
%
% Given the scaled ephemeris data output the satellite position in the ECEF
% coordinate frame.
%
% -------------------------------------------------------------------------
% INPUTS:
%
%       eph       = structured ephemeris date
%       gps_week  = current GPS week
%       ttx       = time of transmission
%
% -------------------------------------------------------------------------
% OUTPUTS:
%
%       X_ECEF    = satellite vehicle (SV) in the ECEF frame 
%
%% GLOBAL CONSTANTS

global mu
omega_e = 0; % Set to zero in this context because we are working in ECI 
             % instead of ECEF coordinates. To work in ECEF, add omega_e to
             % the list of globals.

%% IMPLEMENTATION

% -------------------------------------------------------------------------
% EXTRACT KEPLERIAN ORBITAL ELEMENTS
% -------------------------------------------------------------------------

% Find t_k.
t_k = ttx - eph.Toe;

% Account for week crossovers. Not needed in this analysis. 
% if t_k  > 302400
%     t_k = t_k - 604800;
% elseif t_k < -302400
%     t_k = t_k + 604800;
% end

% Semi-major axis at reference time.
A_0 = eph.Asqrt^2; % [m]
A_k = A_0;

% Mean motion.
n_0 = sqrt(mu / A_0^3); % [rad/s]
n_a = n_0 + eph.Delta_n;

% Eccentricity.
e_n = eph.e;

% Corrected mean anomaly.
M_k = eph.M0 + n_a*t_k;

% Eccentric anomaly.
E_k = Keplers_Eqn(M_k,e_n);

% True anomaly.
tan_a = sqrt(1-e_n^2)*sin(E_k)/(1-e_n*cos(E_k));
tan_b = (cos(E_k)-e_n)/(1-e_n*cos(E_k));

v_k = atan2(tan_a,tan_b);

% Harmonic correction terms.
PHI_k = v_k + eph.Omega;

du_k = eph.Cus*sin(2*PHI_k) + eph.Cuc*cos(2*PHI_k);
dr_k = eph.Crs*sin(2*PHI_k) + eph.Crc*cos(2*PHI_k);
di_k = eph.Cis*sin(2*PHI_k) + eph.Cic*cos(2*PHI_k);

% Argument of latitude.
u_k = PHI_k + du_k;

% Radius.
r_k = A_k*(1-e_n*cos(E_k)) + dr_k;

% Inclination.
i_k = eph.i0 + eph.IDOT*t_k + di_k;

% Longitude of the ascending node.
OMEGA_dot = eph.Omega_dot;
OMEGA_k = eph.Omega0 + (OMEGA_dot - omega_e)*t_k  - omega_e*eph.Toe;


% -------------------------------------------------------------------------
% POSITION OF SATELLITE
% -------------------------------------------------------------------------

% SV position in the satellite orbit.
xk_p = r_k*cos(u_k);
yk_p = r_k*sin(u_k);

% SV position in the ECEF frame. 
xk = xk_p*cos(OMEGA_k) - yk_p*cos(i_k)*sin(OMEGA_k);
yk = xk_p*sin(OMEGA_k) + yk_p*cos(i_k)*cos(OMEGA_k);
zk = yk_p*sin(i_k);

X_ECEF = [xk yk zk];


% -------------------------------------------------------------------------
% CLOCK CORRECTION
% -------------------------------------------------------------------------
% Not used here.
% t = ttx - eph.Toc;
% dtsve = eph.a_f0 + eph.a_f1*t + eph.a_f2*t^2 + F*e_n*eph.Asqrt*sin(E_k);

