function [coeff_A2, coeff_C2, coeff_R2, coeff_T2, coeff_RT, theta] = ...
    analytic_URE_eqn( altitude , R_e )
%% DESCRIPTION
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
%   This function calculates the coefficients based on altitude of the RMS
%   Global URE statistical calculation as outlines in the GPS SPS in
%   formula (A-1). This formulation is based on the method described in: 
%
%   L. Chen, W. Jiao, X. Huang, C. Geng, L. Ai, L. Lu, et al.,
%   "Study on Signal-In-Space Errors Calculation Method and Statistical
%   Characterization of BeiDou Navigation Satellite System," 
%   in China Satellite Navigation Conference (CSNC) 2013 Proceedings: 
%   BeiDou/GNSS Navigation Applications ? Test & Assessment Technology ? 
%   User Terminal Technology, J. Sun, W. Jiao, H. Wu, and C. Shi, Eds., ed 
%   Berlin, Heidelberg: Springer Berlin Heidelberg, 2013, pp. 423-434.
%
% -------------------------------------------------------------------------
% INPUTS
%  
%   altitude - spacecraft altitude [length]
%   R_e      - radius of Earh [length] (consistent units).
%
% -------------------------------------------------------------------------
% OUTPUTS
%
%   These are the coefficients of the URE equation given in the GPS SPS,
%   equation A-1 on page A-19 (4th edition):
%       Global Average URE = ( (T)^2 + (0.980*R)^2 + (0.141*A)^2 + 
%                               +(0.141*C)^2  - 1.960*T*R ) ^ (1/2)
%       where T, R, A, C are the RMS time, radial, along-track, and
%       cross-track error components in [m].
%
%       This is for GPS specifically, more generally this is:
%       Global Average URE = ( coeff_T2 * (cxT)^2 + 
%                              coeff_R2 * R^2 + 
%                              coeff_A2 * A^2 +
%                              coeff_C2 * C^2 - 
%                              ceoff_RT * (c*T*R)
%
%% IMPLEMENTATION

% Set up for decimal conversion.
digits(5);

% Compute the orbital radius.
R_s = ( R_e + altitude ); % []

% Compute theta.
theta = asin(R_e / R_s); % [rad]
% disp(['theta = ', num2str(theta*180/pi),' [deg]'])

% Define symbolic variables.
alpha = sym('alpha');
beta = sym('beta');
T = sym('T');
C = sym('C');
A = sym('A');
R = sym('R');

% Add constraints on these symbolic variables.
assume(alpha,'real')
assume(beta,'real')
assume(T,'real')
assume(C,'real')
assume(A,'real')
assume(R,'real')
assume(R,'positive')

% Define the error vector.
E = [-C, A, -R]';

% Define the l vector. (line of sight vector)
l = [R_e*sin(alpha)*cos(beta), R_e*sin(alpha)*sin(beta), ...
    R_e*cos(alpha) - R_s]' / sqrt(R_e^2 - 2*R_e*R_s*cos(alpha) + R_s^2);

% Compute the RMS SISRE relationship.
% Integrate with respect to beta.
integrand1 = ( E'*l - T ) ^ 2 * sin(alpha) / 2 / pi / (1 - sin(theta));
integrand2 = int(integrand1, beta, 0, 2*pi);

% Integrate with respect to alpha.
F = int(integrand2, alpha, 0, pi/2 - theta);

% Display the equation.
% disp(' ')
% disp('Equation is:')
% vpa(F)

% Get the coefficients of the polynomial.
[coeff, ~] = coeffs(vpa(F),A);
coeff_A2 = double( coeff(1) );

[coeff, ~] = coeffs(vpa(F),C);
coeff_C2 = double( coeff(1) );

[coeff, ~] = coeffs(vpa(F),R);
coeff_R2 = double( coeff(1) );

[coeff, ~] = coeffs(vpa(F),T);
coeff_T2 = double( coeff(1) );

[coeff, ~] = coeffs(coeff(2),R);
coeff_RT = double( coeff(1) );

