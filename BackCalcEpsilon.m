% This script back-calculates epsilon values based on M, R, Sa inputs
% from USGS deaggregations

% Ting Lin
% 6/22/2010

clear all
close all
clc

% Variable definitions
% M_bar         = Moment Magnitude
% T1            = Period (sec); Use Period = -1 for PGV computation
%                 Use 1000 to output the array of median with original period
%                 (no interpolation)
% Rjb           = Joyner-Boore distance (km)
% Fault_Type    = 0 for unspecified fault
%               = 1 for strike-slip fault
%               = 2 for normal fault
%               = 3 for reverse fault
% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China or Turkey
%               = 4 for Italy
% z1            = Basin depth (km); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = 999 if unknown
% Vs30          = shear wave velocity averaged over top 30 m in m/s


% Ground Motion Prediction Model inputs: 
% Fault_Type, region, z1, Vs30 values for the site
Fault_Type = 1;
region = 1;
z1 = 999;
Vs30 = 863;

% USGS hazard and deaggregation inputs:
T1 = 0.5; % Period of interest
M_bar = 6.6;
R_bar = 17.7;
Sa_star = 0.5;

% Assume Rrup and Rjb are the same as R
Rjb = R_bar;

% Obtain median and log standard deviation predictions of Sa given the M/R scenario 
[Sa_1, sigma_1] = BSSA_2014_nga(M_bar, T1, Rjb, Fault_Type, region, z1, Vs30);
% Back-calculate epsilon value, eps_bar, corresponding to target Sa_star
eps_bar = (log(Sa_star) - log(Sa_1))/sigma_1
 
