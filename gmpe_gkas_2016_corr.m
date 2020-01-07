%% Function implementing correlation model by Gulerce, Kamai, Abrahamson, and Silva 2017
function [rho_total, rho_between, rho_within] = gmpe_gkas_2016_corr(Ta, Tb, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj)
%{
Gulerce et al 2016 ground motion correlation model; citation:

Zeynep Gülerce, Ronnie Kamai, Norman A. Abrahamson, 
and Walter J. Silva. "Ground Motion Prediction Equations for
the Vertical Ground Motion Component 
Based on the NGA-W2 Database." Earthquake Spectra 33(2) (2017): 499-528.

Provides correlations between epsilons for vertical (V) component of ground
motion at periods Ta and Tb

By N. Simon Kwong; nealsimonkwong@berkeley.edu

Input Variables
Ta          = 1st Period for vertical component (sec); from 0.01 to 10 sec
Tb          = 2nd Period for vertical component (sec); from 0.01 to 10 sec

M       = Moment magnitude
Rrup    = Closest distance to rupture plane (km)
Rjb     = Closest distance to surface projection of rupture plane (km)
Rx      = Horizontal distance from surface projection of top edge of 
        rupture plane to site, measured perpendicular to average strike, 
        and is negative for footwall but positive for hanging wall (km)
FRV     = Indicator variable representing reverse and reverse-oblique
        faulting (=1 when rake is within 30 to 150 deg)
FNM     = Indicator variable representing normal and normal-oblique
        faulting (=1 when rake is within -30 to -150 deg)
dip     = Average dip angle of rupture plane (deg)
Vs30    = Time-averaged shear wave velocity in top 30m of site (m/s)
region  = 0 for Global
        = 1 for California
        = 2 for Japan and Italy
        = 3 for eastern China
Sj      = Indicator variable representing Japan's site effects

Output Variables
rho_total    = Estimated correlation
%}

%% Specify periods used in GMPM development
T_GKAS2016 = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7 8 10]; % Periods that were used in developing correlation model

%% Execute GMPE for desired T
% Check input period
if Ta>=min(T_GKAS2016) && Ta<=max(T_GKAS2016) && Tb>=min(T_GKAS2016) && Tb<=max(T_GKAS2016)
    if Ta==Tb % Correlation should be unity
        % Get exact values for known periods
        rhoKnown_total = zeros(size(T_GKAS2016)); rhoKnown_between = zeros(size(T_GKAS2016)); rhoKnown_within = zeros(size(T_GKAS2016));
        for ii=1:length(T_GKAS2016)
            [rhoKnown_total(ii), rhoKnown_between(ii), rhoKnown_within(ii)] = GKAS_2016_corr_fixedT(T_GKAS2016(ii), T_GKAS2016(ii), M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);
        end
        % Linear interpolation
        rho_total = interp1(log(T_GKAS2016),rhoKnown_total,log(Ta)); 
        rho_between = interp1(log(T_GKAS2016),rhoKnown_between,log(Ta)); 
        rho_within = interp1(log(T_GKAS2016),rhoKnown_within,log(Ta)); 
    else
        % Determine neighboring periods
        Ta_lo = max(T_GKAS2016(T_GKAS2016<=Ta));
        Ta_up = min(T_GKAS2016(T_GKAS2016>=Ta));
        Tb_lo = max(T_GKAS2016(T_GKAS2016<=Tb));
        Tb_up = min(T_GKAS2016(T_GKAS2016>=Tb));
        
        if Ta_lo==Ta_up && Tb_lo==Tb_up % Unique Ta and Tb
            % Data already available for input periods
            [rho_total, rho_between, rho_within] = GKAS_2016_corr_fixedT(Ta, Tb, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);         
        elseif Ta_lo==Ta_up && Tb_lo~=Tb_up % Unique Ta but non-unique Tb
            [rho_total_lo,  rho_between_lo,  rho_within_lo] = GKAS_2016_corr_fixedT(Ta, Tb_lo, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);  
            [rho_total_up,  rho_between_up,  rho_within_up] = GKAS_2016_corr_fixedT(Ta, Tb_up, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);
            % Linearly interpolate
            rho_total = interp1( log([Tb_lo Tb_up]), [rho_total_lo rho_total_up], log(Tb)); % Impt to interpolate on semilog scale
            rho_between = interp1( log([Tb_lo Tb_up]), [rho_between_lo rho_between_up], log(Tb)); % Impt to interpolate on semilog scale
            rho_within = interp1( log([Tb_lo Tb_up]), [rho_within_lo rho_within_up], log(Tb)); % Impt to interpolate on semilog scale                                   
        elseif Ta_lo~=Ta_up && Tb_lo==Tb_up % Unique Tb but non-unique Ta
            [rho_total_lo,  rho_between_lo,  rho_within_lo] = GKAS_2016_corr_fixedT(Ta_lo, Tb, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);
            [rho_total_up,  rho_between_up,  rho_within_up] = GKAS_2016_corr_fixedT(Ta_up, Tb, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);
            % Linearly interpolate
            rho_total = interp1( log([Ta_lo Ta_up]), [rho_total_lo rho_total_up], log(Ta)); % Impt to interpolate on semilog scale
            rho_between = interp1( log([Ta_lo Ta_up]), [rho_between_lo rho_between_up], log(Ta)); % Impt to interpolate on semilog scale
            rho_within = interp1( log([Ta_lo Ta_up]), [rho_within_lo rho_within_up], log(Ta)); % Impt to interpolate on semilog scale
        else
            % Determine output for neighboring periods
            [X,Y] = meshgrid( [Ta_lo Ta_up], [Tb_lo Tb_up] ); % Create combinations of input periods
            rho_total_approx = zeros(size(X));
            rho_between_approx = zeros(size(X));
            rho_within_approx = zeros(size(X));            
            for ii=1:2
                for jj=1:2
                    xx = X(ii,jj);
                    yy = Y(ii,jj);
                    [rho_total_approx(ii,jj), rho_between_approx(ii,jj), rho_within_approx(ii,jj)] = GKAS_2016_corr_fixedT(xx, yy, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj);
                end
            end
            % Linearly interpolate
            rho_total = interp2(log(X),log(Y),rho_total_approx,log(Ta),log(Tb),'linear'); % Impt to interpolate on semilog scale
            rho_between = interp2(log(X),log(Y),rho_between_approx,log(Ta),log(Tb),'linear'); % Impt to interpolate on semilog scale
            rho_within = interp2(log(X),log(Y),rho_within_approx,log(Ta),log(Tb),'linear'); % Impt to interpolate on semilog scale            
        end
    end
else
    disp(['Input Ta is ' num2str(Ta,4) ' and input Tb is ' num2str(Tb,4)]);
    error('Invalid input periods T');
end
end



function [rho_total, rho_between, rho_within] = GKAS_2016_corr_fixedT(Ta, Tb, M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj)
%% Load previously saved correlation data
% Intra-event data
Table4 = [0.01 1.000 0.995 0.970 0.916 0.900 0.891 0.854 0.805 0.757 0.713 0.636 0.581 0.484 0.423 0.352 0.333 0.328 0.322 0.298 0.264 0.242 0.228 0.249;
0.02 0.995 1.000 0.980 0.926 0.900 0.885 0.841 0.788 0.740 0.696 0.618 0.564 0.467 0.408 0.337 0.319 0.320 0.317 0.294 0.261 0.239 0.226 0.247;
0.03 0.970 0.980 1.000 0.942 0.895 0.861 0.797 0.737 0.686 0.643 0.566 0.513 0.421 0.368 0.302 0.288 0.299 0.305 0.284 0.251 0.228 0.216 0.235;
0.05 0.916 0.926 0.942 1.000 0.911 0.843 0.742 0.670 0.612 0.563 0.483 0.432 0.349 0.303 0.243 0.236 0.262 0.271 0.254 0.222 0.199 0.188 0.207;
0.075 0.900 0.900 0.895 0.911 1.000 0.906 0.776 0.688 0.622 0.569 0.489 0.434 0.351 0.305 0.248 0.238 0.256 0.265 0.247 0.213 0.190 0.175 0.194;
0.1 0.891 0.885 0.861 0.843 0.906 1.000 0.853 0.754 0.682 0.622 0.533 0.473 0.378 0.325 0.261 0.242 0.252 0.253 0.234 0.200 0.176 0.159 0.180;
0.15 0.854 0.841 0.797 0.742 0.776 0.853 1.000 0.876 0.788 0.720 0.626 0.560 0.450 0.383 0.308 0.278 0.268 0.259 0.242 0.209 0.183 0.163 0.185;
0.2 0.805 0.788 0.737 0.670 0.688 0.754 0.876 1.000 0.896 0.821 0.722 0.652 0.530 0.453 0.369 0.325 0.296 0.282 0.263 0.228 0.202 0.177 0.196;
0.25 0.757 0.740 0.686 0.612 0.622 0.682 0.788 0.896 1.000 0.912 0.802 0.727 0.596 0.509 0.410 0.358 0.309 0.291 0.267 0.232 0.207 0.187 0.213;
0.3 0.713 0.696 0.643 0.563 0.569 0.622 0.720 0.821 0.912 1.000 0.864 0.780 0.646 0.554 0.446 0.389 0.328 0.301 0.279 0.247 0.226 0.207 0.233;
0.4 0.636 0.618 0.566 0.483 0.489 0.533 0.626 0.722 0.802 0.864 1.000 0.886 0.736 0.636 0.515 0.448 0.368 0.332 0.306 0.270 0.251 0.228 0.251;
0.5 0.581 0.564 0.513 0.432 0.434 0.473 0.560 0.652 0.727 0.780 0.886 1.000 0.814 0.707 0.576 0.502 0.413 0.375 0.340 0.305 0.286 0.269 0.282;
0.75 0.484 0.467 0.421 0.349 0.351 0.378 0.450 0.530 0.596 0.646 0.736 0.814 1.000 0.853 0.708 0.624 0.502 0.442 0.402 0.364 0.348 0.328 0.342;
1 0.423 0.408 0.368 0.303 0.305 0.325 0.383 0.453 0.509 0.554 0.636 0.707 0.853 1.000 0.807 0.712 0.572 0.495 0.447 0.397 0.375 0.361 0.362;
1.5 0.352 0.337 0.302 0.243 0.248 0.261 0.308 0.369 0.410 0.446 0.515 0.576 0.708 0.807 1.000 0.858 0.685 0.585 0.522 0.462 0.436 0.422 0.406;
2 0.333 0.319 0.288 0.236 0.238 0.242 0.278 0.325 0.358 0.389 0.448 0.502 0.624 0.712 0.858 1.000 0.803 0.692 0.610 0.543 0.510 0.483 0.463;
3 0.328 0.320 0.299 0.262 0.256 0.252 0.268 0.296 0.309 0.328 0.368 0.413 0.502 0.572 0.685 0.803 1.000 0.861 0.763 0.681 0.621 0.573 0.536;
4 0.322 0.317 0.305 0.271 0.265 0.253 0.259 0.282 0.291 0.301 0.332 0.375 0.442 0.495 0.585 0.692 0.861 1.000 0.892 0.796 0.721 0.666 0.612;
5 0.298 0.294 0.284 0.254 0.247 0.234 0.242 0.263 0.267 0.279 0.306 0.340 0.402 0.447 0.522 0.610 0.763 0.892 1.000 0.912 0.826 0.758 0.682;
6 0.264 0.261 0.251 0.222 0.213 0.200 0.209 0.228 0.232 0.247 0.270 0.305 0.364 0.397 0.462 0.543 0.681 0.796 0.912 1.000 0.923 0.844 0.751;
7 0.242 0.239 0.228 0.199 0.190 0.176 0.183 0.202 0.207 0.226 0.251 0.286 0.348 0.375 0.436 0.510 0.621 0.721 0.826 0.923 1.000 0.939 0.826;
8 0.228 0.226 0.216 0.188 0.175 0.159 0.163 0.177 0.187 0.207 0.228 0.269 0.328 0.361 0.422 0.483 0.573 0.666 0.758 0.844 0.939 1.000 0.896;
10 0.249 0.247 0.235 0.207 0.194 0.180 0.185 0.196 0.213 0.233 0.251 0.282 0.342 0.362 0.406 0.463 0.536 0.612 0.682 0.751 0.826 0.896 1.000];
T_GKAS2016 = Table4(:,1); 
corrCoef_intra = Table4(:,2:end);

% Inter-event data
Table5 = [0.01 1.000 0.997 0.984 0.958 0.963 0.960 0.925 0.853 0.755 0.664 0.513 0.395 0.198 0.121 0.044 -0.008 0.028 -0.014 0.034 -0.007 -0.037 -0.067 -0.011;
0.02 0.997 1.000 0.992 0.967 0.963 0.954 0.910 0.831 0.730 0.637 0.483 0.364 0.168 0.091 0.020 -0.030 0.010 -0.024 0.028 -0.007 -0.038 -0.068 -0.012;
0.03 0.984 0.992 1.000 0.983 0.964 0.942 0.873 0.779 0.669 0.569 0.409 0.289 0.090 0.016 -0.047 -0.086 -0.035 -0.053 0.003 -0.029 -0.057 -0.083 -0.030;
0.05 0.958 0.967 0.983 1.000 0.973 0.938 0.843 0.727 0.603 0.496 0.331 0.209 0.012 -0.061 -0.120 -0.155 -0.106 -0.112 -0.046 -0.078 -0.106 -0.133 -0.089;
0.075 0.963 0.963 0.964 0.973 1.000 0.975 0.886 0.776 0.651 0.544 0.377 0.250 0.054 -0.019 -0.083 -0.129 -0.103 -0.133 -0.076 -0.103 -0.125 -0.136 -0.099;
0.1 0.960 0.954 0.942 0.938 0.975 1.000 0.936 0.835 0.718 0.615 0.450 0.329 0.130 0.057 -0.022 -0.079 -0.069 -0.113 -0.058 -0.091 -0.119 -0.148 -0.100;
0.15 0.925 0.910 0.873 0.843 0.886 0.936 1.000 0.947 0.861 0.776 0.626 0.508 0.305 0.226 0.132 0.053 0.041 -0.026 0.015 -0.028 -0.058 -0.094 -0.034;
0.2 0.853 0.831 0.779 0.727 0.776 0.835 0.947 1.000 0.959 0.900 0.769 0.657 0.464 0.384 0.292 0.207 0.199 0.120 0.136 0.083 0.052 0.005 0.065;
0.25 0.755 0.730 0.669 0.603 0.651 0.718 0.861 0.959 1.000 0.968 0.867 0.772 0.602 0.524 0.424 0.347 0.321 0.205 0.203 0.141 0.099 0.051 0.109;
0.3 0.664 0.637 0.569 0.496 0.544 0.615 0.776 0.900 0.968 1.000 0.937 0.861 0.707 0.637 0.538 0.468 0.424 0.287 0.260 0.195 0.158 0.103 0.156;
0.4 0.513 0.483 0.409 0.331 0.377 0.450 0.626 0.769 0.867 0.937 1.000 0.962 0.845 0.777 0.670 0.597 0.528 0.370 0.309 0.240 0.195 0.146 0.185;
0.5 0.395 0.364 0.289 0.209 0.250 0.329 0.508 0.657 0.772 0.861 0.962 1.000 0.925 0.868 0.764 0.692 0.630 0.476 0.402 0.332 0.300 0.228 0.262;
0.75 0.198 0.168 0.090 0.012 0.054 0.130 0.305 0.464 0.602 0.707 0.845 0.925 1.000 0.967 0.880 0.816 0.746 0.597 0.505 0.441 0.417 0.350 0.359;
1 0.121 0.091 0.016 -0.061 -0.019 0.057 0.226 0.384 0.524 0.637 0.777 0.868 0.967 1.000 0.941 0.884 0.817 0.664 0.553 0.479 0.453 0.369 0.369;
1.5 0.044 0.020 -0.047 -0.120 -0.083 -0.022 0.132 0.292 0.424 0.538 0.670 0.764 0.880 0.941 1.000 0.960 0.881 0.755 0.633 0.567 0.554 0.478 0.441;
2 -0.008 -0.030 -0.086 -0.155 -0.129 -0.079 0.053 0.207 0.347 0.468 0.597 0.692 0.816 0.884 0.960 1.000 0.929 0.798 0.685 0.612 0.590 0.523 0.484;
3 0.028 0.010 -0.035 -0.106 -0.103 -0.069 0.041 0.199 0.321 0.424 0.528 0.630 0.746 0.817 0.881 0.929 1.000 0.923 0.833 0.747 0.697 0.610 0.569;
4 -0.014 -0.024 -0.053 -0.112 -0.133 -0.113 -0.026 0.120 0.205 0.287 0.370 0.476 0.597 0.664 0.755 0.798 0.923 1.000 0.944 0.866 0.795 0.717 0.669;
5 0.034 0.028 0.003 -0.046 -0.076 -0.058 0.015 0.136 0.203 0.260 0.309 0.402 0.505 0.553 0.633 0.685 0.833 0.944 1.000 0.955 0.881 0.806 0.747;
6 -0.007 -0.007 -0.029 -0.078 -0.103 -0.091 -0.028 0.083 0.141 0.195 0.240 0.332 0.441 0.479 0.567 0.612 0.747 0.866 0.955 1.000 0.966 0.918 0.842;
7 -0.037 -0.038 -0.057 -0.106 -0.125 -0.119 -0.058 0.052 0.099 0.158 0.195 0.300 0.417 0.453 0.554 0.590 0.697 0.795 0.881 0.966 1.000 0.975 0.915;
8 -0.067 -0.068 -0.083 -0.133 -0.136 -0.148 -0.094 0.005 0.051 0.103 0.146 0.228 0.350 0.369 0.478 0.523 0.610 0.717 0.806 0.918 0.975 1.000 0.952;
10 -0.011 -0.012 -0.030 -0.089 -0.099 -0.100 -0.034 0.065 0.109 0.156 0.185 0.262 0.359 0.369 0.441 0.484 0.569 0.669 0.747 0.842 0.915 0.952 1.000];
corrCoef_inter = Table5(:,2:end);

%% Sigmas from GMPM for V component
% Std dev for V at Ta
[~, sig_Va, tau_Va, phi_Va] = gmpeV_bc_2016(M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj, Ta);

% Std dev for V at Tb
[~, sig_Vb, tau_Vb, phi_Vb] = gmpeV_bc_2016(M, Rrup, Rjb, Rx, FRV, FNM, dip, Vs30, region, Sj, Tb);

%% Compute correlation
id_Ta = find( T_GKAS2016 == Ta );
id_Tb = find( T_GKAS2016 == Tb );
rho_within = corrCoef_intra( id_Tb, id_Ta );
rho_between = corrCoef_inter( id_Tb, id_Ta );
rho_total = (phi_Va*phi_Vb)/(sig_Va*sig_Vb)*rho_within + (tau_Va*tau_Vb)/(sig_Va*sig_Vb)*rho_between;

end
