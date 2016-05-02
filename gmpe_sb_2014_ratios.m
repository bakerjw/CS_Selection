function [ ratio, sigma, phi, tau ] = gmpeSB_2014_ratios( T )
% Created by Jack Baker, August 21, 2015
% Updated 12 October, 2015, to include Phi and Tau values
%
% Compute Sa_RotD100/Sa_RotD50 ratios, from the following model:
%
%   Shahi, S. K., and Baker, J. W. (2014). "NGA-West2 models for ground-
%   motion directionality." Earthquake Spectra, 30(3), 1285-1300.
%
% INPUTS:
%
%   T         = Period(s) of interest (sec) 
%
% OUTPUTS:
%
%   ratio         = geometric mean of Sa_RotD100/Sa_RotD50
%
%   sigma         = standard deviation of log(Sa_RotD100/Sa_RotD50)
%
%   phi           = within-event standard deviation
%
%   tau           = between-event standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Model coefficient values from Table 1 of the above-reference paper
periods_orig =  [0.0100000000000000,0.0200000000000000, 0.0300000000000000, 0.0500000000000000, 0.0750000000000000, 0.100000000000000,  0.150000000000000,  0.200000000000000,  0.250000000000000,  0.300000000000000,  0.400000000000000,  0.500000000000000,  0.750000000000000,  1,                  1.50000000000000,   2,                  3,                  4,                  5,              7.50000000000000,   10];
ratios_orig =   [1.19243805900000,  1.19124621700000,   1.18767783300000,   1.18649074900000,   1.18767783300000,   1.18767783300000,   1.19961419400000,   1.20562728500000,   1.21652690500000,   1.21896239400000,   1.22875320400000,   1.22875320400000,   1.23738465100000,   1.24110237900000,   1.24234410200000,   1.24358706800000,   1.24732343100000,   1.25985923900000,   1.264908769000, 1.28531008400000,   1.29433881900000];
sigma_orig =    [0.08,              0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,           0.08,               0.08];
phi_orig =      [0.08,              0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,               0.08,           0.08,               0.07];
tau_orig =      [0.01,              0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.01,               0.02,           0.02,               0.03];


% Interpolate to compute values for the user-specified periods
ratio = interp1(log(periods_orig), ratios_orig, log(T));
sigma = interp1(log(periods_orig), sigma_orig, log(T));
phi = interp1(log(periods_orig), phi_orig, log(T));
tau = interp1(log(periods_orig), tau_orig, log(T));


end

