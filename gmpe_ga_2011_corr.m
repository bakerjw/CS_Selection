%% Function implementing correlation model by Gulerce and Abrahamson 2011
function [rho_total, rho_between, rho_within] = gmpe_ga_2011_corr(T_H, T_VH, M, Rrup, Vs30, FRV, FNM)
%{
Gulerce and Abrahamson 2011 ground motion prediction model; citation:

Gülerce, Zeynep, and Norman A. Abrahamson. "Site-specific design spectra 
for vertical ground motion." Earthquake Spectra 27(4) (2011): 1023-1047.

Provides correlations between epsilon at T_{H} for horizontal (H) component
of ground motion and epsilon at T_{V/H} for vertical-to-horizontal (V/H)
ratio

By N. Simon Kwong; nealsimonkwong@berkeley.edu

Input Variables
T_H          = Period for H component (sec); from 0.01 to 10 sec
T_VH         = Period for V/H ratio (sec); from 0.01 to 10 sec
M            = Magnitude
Rrup         = Closest distance to plane of rupture (km)
Vs30         = Average shear wave velocity (m/s) in upper 30m of soil
FRV          = Flag for reverse faulting
FNM          = Flag for normal faulting

Output Variables
rho_total    = Estimated correlation
%}

%% Specify periods used in GMPM development
T_GA2011 = [0.01 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1 2 3 4 5 10]; % Periods that were used in developing correlation model

%% Execute GMPE for desired T
% Check input period
if T_H>=min(T_GA2011) && T_H<=max(T_GA2011) && T_VH>=min(T_GA2011) && T_VH<=max(T_GA2011)
    if T_H==T_VH
        % Get exact values for known periods
        rhoKnown_total = zeros(size(T_GA2011)); rhoKnown_between = zeros(size(T_GA2011)); rhoKnown_within = zeros(size(T_GA2011));
        for ii=1:length(T_GA2011)
            [rhoKnown_total(ii), rhoKnown_between(ii), rhoKnown_within(ii)] = GA_2011_corr_fixedT(T_GA2011(ii),T_GA2011(ii),M,Rrup,Vs30,FRV,FNM);
        end
        % Linear interpolation
        rho_total = interp1(log(T_GA2011),rhoKnown_total,log(T_H)); % T_H = T_V and impt to interpolate on semilog scale
        rho_between = interp1(log(T_GA2011),rhoKnown_between,log(T_H)); % T_H = T_V and impt to interpolate on semilog scale
        rho_within = interp1(log(T_GA2011),rhoKnown_within,log(T_H)); % T_H = T_V and impt to interpolate on semilog scale
    else
        % Determine neighboring periods
        T_H_lo = max(T_GA2011(T_GA2011<=T_H));
        T_H_up = min(T_GA2011(T_GA2011>=T_H));
        T_VH_lo = max(T_GA2011(T_GA2011<=T_VH));
        T_VH_up = min(T_GA2011(T_GA2011>=T_VH));
        
        if T_H_lo==T_H_up && T_VH_lo==T_VH_up % Unique T_H and T_VH
            % Data already available for input periods
            [rho_total, rho_between, rho_within] = GA_2011_corr_fixedT(T_H,T_VH,M,Rrup,Vs30,FRV,FNM);            
        elseif T_H_lo==T_H_up && T_VH_lo~=T_VH_up % Unique T_H but non-unique T_VH
            [rho_total_lo,  rho_between_lo,  rho_within_lo] = GA_2011_corr_fixedT(T_H,T_VH_lo,M,Rrup,Vs30,FRV,FNM);
            [rho_total_up,  rho_between_up,  rho_within_up] = GA_2011_corr_fixedT(T_H,T_VH_up,M,Rrup,Vs30,FRV,FNM);
            % Linearly interpolate
            rho_total = interp1( log([T_VH_lo T_VH_up]), [rho_total_lo rho_total_up], log(T_VH)); % Impt to interpolate on semilog scale
            rho_between = interp1( log([T_VH_lo T_VH_up]), [rho_between_lo rho_between_up], log(T_VH)); % Impt to interpolate on semilog scale
            rho_within = interp1( log([T_VH_lo T_VH_up]), [rho_within_lo rho_within_up], log(T_VH)); % Impt to interpolate on semilog scale                                   
        elseif T_H_lo~=T_H_up && T_VH_lo==T_VH_up % Unique T_VH but non-unique T_H
            [rho_total_lo,  rho_between_lo,  rho_within_lo] = GA_2011_corr_fixedT(T_H_lo,T_VH,M,Rrup,Vs30,FRV,FNM);
            [rho_total_up,  rho_between_up,  rho_within_up] = GA_2011_corr_fixedT(T_H_up,T_VH,M,Rrup,Vs30,FRV,FNM);
            % Linearly interpolate
            rho_total = interp1( log([T_H_lo T_H_up]), [rho_total_lo rho_total_up], log(T_H)); % Impt to interpolate on semilog scale
            rho_between = interp1( log([T_H_lo T_H_up]), [rho_between_lo rho_between_up], log(T_H)); % Impt to interpolate on semilog scale
            rho_within = interp1( log([T_H_lo T_H_up]), [rho_within_lo rho_within_up], log(T_H)); % Impt to interpolate on semilog scale
        else
            % Determine output for neighboring periods
            [X,Y] = meshgrid( [T_H_lo T_H_up], [T_VH_lo T_VH_up] ); % Create combinations of input periods
            rho_total_approx = zeros(size(X));
            rho_between_approx = zeros(size(X));
            rho_within_approx = zeros(size(X));            
            for ii=1:2
                for jj=1:2
                    xx = X(ii,jj);
                    yy = Y(ii,jj);
                    [rho_total_approx(ii,jj), rho_between_approx(ii,jj), rho_within_approx(ii,jj)] = GA_2011_corr_fixedT(xx,yy,M,Rrup,Vs30,FRV,FNM);
                end
            end
            % Linearly interpolate
            rho_total = interp2(log(X),log(Y),rho_total_approx,log(T_H),log(T_VH),'linear'); % Impt to interpolate on semilog scale
            rho_between = interp2(log(X),log(Y),rho_between_approx,log(T_H),log(T_VH),'linear'); % Impt to interpolate on semilog scale
            rho_within = interp2(log(X),log(Y),rho_within_approx,log(T_H),log(T_VH),'linear'); % Impt to interpolate on semilog scale            
        end
    end
else
    disp(['Input T_H is ' num2str(T_H,4) ' and input T_VH is ' num2str(T_VH,4)]);
    error('Invalid input periods T');
end
end



function [rho_total, rho_between, rho_within] = GA_2011_corr_fixedT(T_H,T_VH,M,Rrup,Vs30,FRV,FNM)
%% Load previously saved correlation data
% Intra-event data
Table4 = [0.01 -0.389 -0.312 -0.326 -0.317 -0.331 -0.332 -0.325 -0.353 -0.311 -0.259 -0.130 -0.094 -0.116 -0.152 0.030
0.05 -0.324 -0.273 -0.328 -0.310 -0.318 -0.301 -0.287 -0.303 -0.271 -0.216 -0.109 -0.073 -0.099 -0.139 0.067
0.1 -0.264 -0.173 -0.358 -0.344 -0.330 -0.274 -0.233 -0.243 -0.205 -0.165 -0.072 -0.027 -0.043 -0.079 0.063
0.15 -0.286 -0.180 -0.291 -0.417 -0.380 -0.311 -0.273 -0.269 -0.208 -0.173 -0.070 -0.020 -0.047 -0.067 0.043
0.2 -0.311 -0.216 -0.258 -0.325 -0.439 -0.348 -0.299 -0.300 -0.233 -0.180 -0.074 -0.038 -0.056 -0.077 0.033
0.3 -0.347 -0.269 -0.234 -0.240 -0.289 -0.450 -0.377 -0.370 -0.296 -0.236 -0.097 -0.065 -0.077 -0.102 0.014
0.4 -0.352 -0.291 -0.221 -0.200 -0.219 -0.333 -0.466 -0.431 -0.338 -0.277 -0.137 -0.106 -0.104 -0.144 0.011
0.5 -0.346 -0.296 -0.204 -0.162 -0.183 -0.269 -0.361 -0.493 -0.398 -0.327 -0.166 -0.141 -0.147 -0.185 -0.005
0.75 -0.287 -0.250 -0.160 -0.105 -0.119 -0.176 -0.232 -0.329 -0.483 -0.396 -0.217 -0.168 -0.170 -0.204 -0.020
1 -0.247 -0.216 -0.126 -0.086 -0.095 -0.140 -0.175 -0.258 -0.350 -0.443 -0.236 -0.191 -0.194 -0.206 -0.058
2 -0.162 -0.108 -0.049 -0.057 -0.072 -0.098 -0.131 -0.193 -0.226 -0.254 -0.318 -0.220 -0.247 -0.324 -0.260
3 -0.115 -0.060 -0.004 -0.018 -0.017 -0.055 -0.066 -0.151 -0.176 -0.208 -0.210 -0.254 -0.255 -0.324 -0.330
4 -0.080 -0.019 0.046 -0.003 -0.020 -0.044 -0.041 -0.138 -0.152 -0.173 -0.151 -0.125 -0.295 -0.343 -0.407
5 -0.072 -0.009 0.061 0.009 -0.009 -0.009 -0.032 -0.112 -0.098 -0.122 -0.119 -0.060 -0.224 -0.396 -0.511
10 -0.083 -0.020 0.042 -0.013 -0.061 -0.085 -0.086 -0.199 -0.196 -0.179 -0.163 -0.071 -0.205 -0.309 -0.577];
T_GA2011 = Table4(:,1); % Periods for V/H component same as those for H comp
corrCoef_intra = Table4(:,2:end);

% Inter-event data
Table5 = [0.01 -0.193 -0.156 -0.210 -0.320 -0.339 -0.314 -0.257 -0.147 -0.162 -0.061 -0.087 -0.016 0.036 -0.018 -0.249
0.05 -0.087 -0.068 -0.192 -0.306 -0.259 -0.238 -0.150 -0.036 -0.077 0.043 -0.001 0.045 0.106 0.063 -0.002
0.1 -0.072 -0.016 -0.186 -0.305 -0.325 -0.216 -0.144 -0.047 -0.091 0.044 -0.028 0.058 0.100 0.017 -0.265
0.15 -0.110 -0.048 -0.159 -0.314 -0.382 -0.292 -0.253 -0.157 -0.168 -0.033 -0.080 0.030 0.101 0.027 -0.341
0.2 -0.145 -0.101 -0.145 -0.287 -0.413 -0.321 -0.282 -0.173 -0.175 -0.072 -0.117 -0.021 0.035 -0.065 -0.390
0.3 -0.216 -0.165 -0.163 -0.249 -0.369 -0.418 -0.379 -0.239 -0.233 -0.163 -0.179 -0.073 -0.015 -0.070 -0.315
0.4 -0.237 -0.157 -0.108 -0.196 -0.331 -0.350 -0.416 -0.308 -0.300 -0.218 -0.231 -0.143 -0.088 -0.127 -0.413
0.5 -0.205 -0.131 -0.061 -0.166 -0.289 -0.300 -0.382 -0.351 -0.358 -0.264 -0.286 -0.224 -0.126 -0.148 -0.420
0.75 -0.187 -0.188 -0.061 -0.114 -0.235 -0.223 -0.242 -0.231 -0.370 -0.276 -0.338 -0.303 -0.219 -0.171 -0.365
1 -0.167 -0.247 0.018 -0.070 -0.176 -0.204 -0.210 -0.224 -0.258 -0.245 -0.286 -0.261 -0.226 -0.226 -0.375
2 -0.230 -0.325 -0.080 0.006 -0.061 -0.147 -0.200 -0.167 -0.254 -0.241 -0.361 -0.339 -0.304 -0.232 -0.437
3 -0.262 -0.302 -0.168 -0.166 -0.199 -0.169 -0.242 -0.223 -0.348 -0.302 -0.349 -0.358 -0.252 -0.204 -0.292
4 -0.209 -0.268 -0.129 -0.118 -0.175 -0.160 -0.297 -0.258 -0.303 -0.253 -0.346 -0.303 -0.181 -0.161 -0.266
5 -0.345 -0.362 -0.247 -0.240 -0.354 -0.287 -0.367 -0.309 -0.323 -0.275 -0.402 -0.324 -0.198 -0.200 -0.366
10 -0.268 -0.301 -0.289 -0.250 -0.268 -0.192 0.019 0.062 0.014 0.039 -0.182 -0.074 -0.037 -0.087 -0.171];
corrCoef_inter = Table5(:,2:end);

%% Infer other seismological parameters for AS2008 GMPM
% Assume vertical strike-slip faulting mechanism
Rjb = Rrup;
Rx = Rrup;
dip = 90; 
Ztor = 0; 
Z10 = 0; 
W = 15; 
FAS = 0;
FHW = 0;
FVS30 = 0;

%% Sigmas from GMPMs
% Std dev for H
[~, sig_H, ~, ~, phi_H, tau_H] = gmpe_as_2008_mod_wCD(M, Vs30, T_H, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30);

% Std dev for V/H
[~, sig_VH, phi_VH, tau_VH] = gmpeVH_ga_2011(M,Rrup,Vs30,FRV,FNM,T_VH);

%% Compute correlation
id_Th = find( T_GA2011 == T_H );
id_Tvh = find( T_GA2011 == T_VH );
rho_within = corrCoef_intra( id_Tvh, id_Th );
rho_between = corrCoef_inter( id_Tvh, id_Th );
rho_total = (phi_H*phi_VH)/(sig_H*sig_VH)*rho_within + (tau_H*tau_VH)/(sig_H*sig_VH)*rho_between;

end
