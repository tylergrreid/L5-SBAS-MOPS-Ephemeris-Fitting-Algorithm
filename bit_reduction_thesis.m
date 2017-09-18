function [ab,eb,incb,RAANb,omegab,M0b,Cusb,Cucb,IDOTb,numbits] = bit_reduction(a,e,inc,RAAN,omega,M0,Cus,Cuc,IDOT)
%% BIT REDUCTION
% -------------------------------------------------------------------------
% DESCRIPTION
%
%  Written By: Tyler Reid (tyreid@stanford.edu)
%  PI: Todd Walter, Per Enge
%  Lab: Stanford University GPS Lab
%  Date Start: Summer 2013
%  Modified to match IWG draft ICD March 3, 2015
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% This takes the L5 SBAS MOPS ephemeris message parameters and applies
% quantization. This reduces the number of bits used in matlab to represent
% these parameters to those specified in the ICD. 
% 
% -------------------------------------------------------------------------
% INPUT:
%   
%       a         = Best fit semi-major axis [m]
%       e         = Best fit eccentricity [-]
%       inc       = Best fit inclination [rad]
%       RAAN      = Best fit right ascension of the ascending node [rad]
%       omega     = Best fit argument of perigee [rad]
%       M0        = Best fit mean anomaly at epoch [rad]
%       Cus, Cuc  = Best fit along-track harmonic correction terms [rad/s]
%       IDOT      = Best fit inclination correction rate [rad/s]
%
% ------------------------------------------------------------------------- 
% OUTPUT:
%      
%       ab         = Quantized (truncated) a
%       eb         = Quantized (truncated) e
%       incb       = Quantized (truncated) inc
%       RAANb      = Quantized (truncated) RAAN
%       omegab     = Quantized (truncated) omega
%       M0b        = Quantized (truncated) M0
%       Cusb, Cucb = Quantized (truncated) Cus, Cuc
%       IDOTb      = Quantized (truncated) IDOT
%       num_bits   = Total number of message bits (combined)
%
%% IMPLEMENTATION

% Parameter bit assignment.
% Based on the IWG draft MOPS ICD 16/12/2014 
% This was different for ION GNSS 2013.
a_bits        = 32;
e_bits        = 31;
i_bits        = 34; 
RAAN_bits     = 35; % allows 1 sign bit
omega_bits    = 35; % allows 1 sign bit
M0_bits       = 35; % allows 1 sign bit
Cuc_bits      = 22; % allows 1 sign bit
Cus_bits      = 22; % allows 1 sign bit
IDOT_bits     = 22; % allows 1 sign bit

% Compute the number of message bits. 
numbits = a_bits+e_bits+i_bits+RAAN_bits+omega_bits+M0_bits+IDOT_bits+Cuc_bits+Cus_bits;

% Define scale factors (SF).
% SF's are as per the IWG draft MOPS ICD 16/12/2014 
% This was different for ION GNSS 2013. 
scale_a     = 0.01; % 32 bits (0 sign)
scale_e     = 2^-31; % 33 bits (0 sign)
scale_i     = pi*2^-34; % 34 bits (0 sign)
scale_RAAN  = pi*2^-34; % 35 bits (1 sign)
scale_omega = pi*2^-34; % 35 bits (1 sign)
scale_M0    = pi*2^-34; % 35 bits (1 sign)
scale_Cuc   = ((pi/2)*1e-4)*(2^-21); % 22 bits (1 sign)
scale_Cus   = ((pi/2)*1e-4)*(2^-21); % 22 bits (1 sign)
scale_IDOT  = ((7*pi/6)*1e-6)*(2^-21); % 22 bits (1 sign)

% -------------------------------------------------------------------------
% Semi major axis
% -------------------------------------------------------------------------

% Apply scale factor.
temp     = a/scale_a;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within bit limit.
if length(tempbin)>a_bits
    error('Semi major axis bit overflow')
    return;
else
    ab = bin2dec(tempbin)*scale_a;
    a_bits = length(tempbin);
end

% -------------------------------------------------------------------------
% Eccentricity
% -------------------------------------------------------------------------

% Apply scale factor.
temp     = e/scale_e;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within bit limit.
if length(tempbin)>e_bits
    error('Eccentricity bit overflow')
    return;
else
    eb = bin2dec(tempbin)*scale_e;
    e_bits = length(tempbin);
end


% -------------------------------------------------------------------------
% Inclination
% -------------------------------------------------------------------------

% Make sure inclination is in acceptable limits.
if inc < 0 || inc > pi
    fprintf('Inclination must be a number between 0 and pi [rad]\n');
    return;
end

% Apply scale factor.
temp     = inc/scale_i;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within the bit limit.
if length(tempbin)>i_bits
    error('Inclination bit overflow')
    return;
else
    incb = bin2dec(tempbin)*scale_i;
    i_bits = length(tempbin);
end


% -------------------------------------------------------------------------
% RAAN
% -------------------------------------------------------------------------

% Make sure RAAN is within +/- pi bounds.
if RAAN > pi || RAAN < -pi
    % Put in the correct range. 
    RAAN = mod(RAAN,2*pi);
    
    if RAAN > pi
        RAAN = RAAN-2*pi;
    end
end

% Deal with negative numbers.
if RAAN < 0
    RAAN   = -RAAN;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = RAAN/scale_RAAN;

% Convert to binary string.
tempbin  = dec2bin(floor(temp));

% Determine if the string is within the bit limit.
if length(tempbin)>RAAN_bits-1 % Allows for theoretical sign bit
    error('RAAN bit overflow')
    return;
else
    RAANb = bin2dec(tempbin)*scale_RAAN*(-1)^negate;
    RAAN_bits = length(tempbin);
end

% -------------------------------------------------------------------------
% omega
% -------------------------------------------------------------------------

% Make sure omega is within +/- pi bounds
if omega >= pi || omega <= -pi
    % Put in the correct range. 
    omega = mod(omega, 2*pi);
    
    if omega >= pi
        omega = omega-2*pi;
    end
end

% Deal with negative numbers.
if omega <= 0
    omega   = -omega;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = omega/scale_omega;
if omega == pi
    temp = temp - 1;
end

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within bit limit.
if length(tempbin)>omega_bits-1 % allows for theoretical sign bit
    error('omega bit overflow')
    return;
else
    omegab = bin2dec(tempbin)*scale_omega*(-1)^negate;
    omega_bits = length(tempbin);
end


% -------------------------------------------------------------------------
% M0
% -------------------------------------------------------------------------

% Make sure omega is within +/- pi bounds.
if M0 > pi || M0 < -pi
    % Put in the correct range.
    M0 = mod(M0,2*pi);
    
    if M0 > pi
        M0 = M0-2*pi;
    end
end

% Deal with negative numbers.
if M0 < 0
    M0   = -M0;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = M0/scale_M0;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within bit limit.
if length(tempbin)>M0_bits-1 % allows for theoretical sign bit
    error('M0 bit overflow')
    return;
else
    M0b = bin2dec(tempbin)*scale_M0*(-1)^negate;
    M0_bits = length(tempbin);
end


% -------------------------------------------------------------------------
% IDOT
% -------------------------------------------------------------------------

% Deal with negative numbers.
if IDOT < 0
    IDOT   = -IDOT;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = IDOT/scale_IDOT;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is within bit limit.
if length(tempbin)>IDOT_bits-1 % allows for theoretical sign bit
    error('IDOT bit overflow')
    return;
else
    IDOTb = bin2dec(tempbin)*scale_IDOT*(-1)^negate;
    IDOT_bits = length(tempbin);
end


% -------------------------------------------------------------------------
% Cus
% -------------------------------------------------------------------------

% Deal with negative numbers.
if Cus < 0
    Cus   = -Cus;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = Cus/scale_Cus;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is withing bit limit.
if length(tempbin)>Cus_bits-1 % allows for theoretical sign bit
    error('Cus bit overflow')
    return;
else
    Cusb = bin2dec(tempbin)*scale_Cus*(-1)^negate;
    Cus_bits = length(tempbin);
end



% -------------------------------------------------------------------------
% Cuc
% -------------------------------------------------------------------------

% Deal with negative numbers.
if Cuc < 0
    Cuc   = -Cuc;
    negate = 1;
else
    negate = 0;
end

% Apply scale factor.
temp     = Cuc/scale_Cuc;

% Convert to binary string.
tempbin  = dec2bin(round(temp));

% Determine if the string is withing bit limit.
if length(tempbin)>Cuc_bits-1 % allows for theoretical sign bit
    error('Cuc bit overflow')
    return;
else
    Cucb = bin2dec(tempbin)*scale_Cuc*(-1)^negate;
    Cuc_bits = length(tempbin);
end

% Count the total number of bits of the message.
numbits = a_bits + e_bits + i_bits + RAAN_bits + 1 + omega_bits +...
    1 + M0_bits + 1 + IDOT_bits + 1 + Cuc_bits + 1 + Cus_bits;
