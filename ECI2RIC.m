function ECI_2_RIC_Mat = ECI2RIC( x_eci, v_eci )
%% DESCRIPTION
%       
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       PI:                   Todd Walter, Per Enge
%       Lab:                  Stanford University GPS Lab
%       Start Date:           April 12, 2016
%       Last Modified:        April 12, 2016
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Give the position and velocity column vectors of the spacecraft, output
% the transformation matrix from ECI to radial, along-track, and
% cross-track coordinates.
%
% -------------------------------------------------------------------------
% INPUT:
%   
%       x_eci = ECI position vector of the spacecraft    [length]*
%       v_eci = ECI velocity vector of the spacecraft    [length/time]*
%
% ------------------------------------------------------------------------- 
% OUTPUT:
%      
%       ECI_2_RIC_Mat = Transformation matrix from ECI to local
%       coordinates (radial, along-track, and cross-track). Where we have
%       x_ric = ECI_2_RIC_Mat * x_eci
%
% -------------------------------------------------------------------------
% NOTES:
%
% (1) This quantity can be expressed in any length/time unit, the output 
%     will be consistant with the input.
% (2) Position and velocity inputs are assumed to be column vectors.
%
%% IMPLEMENTATION

% Radial.
r_hat = x_eci / norm(x_eci);

% Cross Track.
c_hat = cross( x_eci, v_eci );
c_hat = c_hat / norm( c_hat );

% In / Along Track.
i_hat = cross( c_hat, r_hat );

% Return the transformation matrix.
ECI_2_RIC_Mat = [ r_hat'; i_hat'; c_hat'];
