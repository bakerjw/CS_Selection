function [rho] = baker_jayaram_correlation(T1, T2);

% Created by Jack Baker, 2/28/07 (updated 6/25/2007)
% Compute the correlation of epsilons for the NGA ground motion models
%
% The function is strictly emperical, fitted over the range the range 0.01s <= T1, T2 <= 10s
%
% Documentation is provided in the following document:
% Baker, J.W. and Jayaram, N., "Correlation of spectral acceleration values from NGA ground 
% motion models," Earthquake Spectra, (in review).

% INPUT
%
%   T1, T2      = The two periods of interest. The periods may be equal,
%                 and there is no restriction on which one is larger.
%
% INPUT
%
%   rho         = The predicted correlation coefficient


 
        T_min = min(T1, T2);
        T_max = max(T1, T2);
        
        C1 = (1-cos(pi/2 - log(T_max/max(T_min, 0.109)) * 0.366 ));
        if T_max < 0.2
            C2 = 1 - 0.105*(1 - 1./(1+exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099);
        end
        if T_max < 0.109
            C3 = C2;
        else
            C3 = C1;
        end
        C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi*(T_min)/(0.109)));

        if T_max <= 0.109
            rho = C2;
        elseif T_min > 0.109
            rho = C1;
        elseif T_max < 0.2
            rho = min(C2, C4);
        else
            rho = C4;
        end

        
