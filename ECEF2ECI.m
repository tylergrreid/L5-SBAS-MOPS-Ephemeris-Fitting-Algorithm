function [X_ECI,V_ECI] = ECEF2ECI(X_ECEF,V_ECEF,GMST)
%% DESCRIPTION
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       Lab:                  Stanford GPS Lab
%       Project Title:        Arctic Navigation / WAAS
%       Project Start Date:   March 28, 2011
%       Last updated:         April 23, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Given the X,Y,Z coordinates of a spacecraft in ECEF coordintes and the
% Greenwich Mean Sidereal Time (GMST), compute the position vector of the
% spacecraft in the ECI frame.
%
% -------------------------------------------------------------------------
% INPUT:
%   
%       X_ECEF = ECEF position vector of the spacecraft    [length]*
%       V_ECEF = ECEF velocity vector of the spacecraft    [length/time]*
%         GMST = Greenwich mean sidereal time              [rad]
%
% ------------------------------------------------------------------------- 
% OUTPUT:
%      
%       X_ECI = ECI position vector of the spacecraft      [length]*
%       V_ECI = ECI velocity vector of the spacecraft      [length/time]*
%
% -------------------------------------------------------------------------
% NOTES:
%
% (1) this quantity can be expressed in any length/time unit, the output 
%     will be consistant with the input
% (2) position and velocity inputs are assumed to be column vectors
%
%% DEFINE GLOBAL VARIABLES USED

global omega_e

%% IMPLEMENTATION

ECI_C_ECEF = [cos(GMST) -sin(GMST) 0;sin(GMST) cos(GMST) 0;0 0 1] ;

X_ECI = ECI_C_ECEF*X_ECEF;
V_ECI = ECI_C_ECEF*V_ECEF + cross([0 0 omega_e]',X_ECI);

