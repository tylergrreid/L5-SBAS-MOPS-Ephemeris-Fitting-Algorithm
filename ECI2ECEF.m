function [X_ECEF,V_ECEF] = ECI2ECEF(X_ECI,V_ECI,GMST)
%% DESCRIPTION
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Project Title:        Arctic Navigation / WAAS
%       Project Start Date:   March 28, 2011
%       Last updated:         April 18, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Given the X,Y,Z coordinates of a spacecraft in ECI coordintes and the
% Greenwich Mean Sidereal Time (GMST), compute the position vector of the
% spacecraft in the ECEF frame.
%
% -------------------------------------------------------------------------
% INPUT:
%   
%       X_ECI = ECI position vector of the spacecraft    [length]*
%       V_ECI = ECI velocity vector of the spacecraft    [length/time]*
%        GMST = Greenwich mean sidereal time             [rad]
%
% ------------------------------------------------------------------------- 
%
% OUTPUT:
%      
%       X_ECEF = ECEF position vector of the spacecraft   [length]*
%       V_ECI = ECI velocity vector of the spacecraft     [length/time]*
%
% -------------------------------------------------------------------------
%
% NOTES:
%
% * this quantity can be expressed in any length/time unit, the output 
%   will be consistant with the input
%
%% DEFINE GLOBAL VARIABLES USED

global omega_e

%% IMPLEMENTATION

% define the tranformation matrix 
ECEF_C_ECI = [cos(GMST) sin(GMST) 0;-sin(GMST) cos(GMST) 0;0 0 1] ;

% determine if input arguments are column or row vectors
[m n] = size(X_ECI);

if m>n % input arguments are column vectors
    X_ECEF = ECEF_C_ECI*X_ECI;
    V_ECEF = ECEF_C_ECI*V_ECI;
elseif n>m % input arguments are row vectors 
    X_ECEF = X_ECI*ECEF_C_ECI';
    V_ECEF = V_ECI*ECEF_C_ECI';
else
    fprintf('Error in ECI2ECEF - matrix dimensions\n')
end

