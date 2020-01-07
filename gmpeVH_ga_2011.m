%% Function implementing GMPM by Gulerce and Abrahamson 2011
function [Y_vec, sig_lnY_vec, sig_0_vec, tau_0_vec] = gmpeVH_ga_2011(M, Rrup, Vs30, FRV, FNM, T_vec)
%{
Gulerce and Abrahamson 2011 ground motion prediction model; citation:

Gülerce, Zeynep, and Norman A. Abrahamson. "Site-specific design spectra 
for vertical ground motion." Earthquake Spectra 27(4) (2011): 1023-1047.

Provides GM predictions for medians and logarithmic standard deviations of 
V/H ratios of PGA, PGV, and 5% damped spectral acceleration

By N. Simon Kwong; nealsimonkwong@berkeley.edu

Input Variables
M            = Magnitude
Rrup         = Closest distance to plane of rupture (km)
Vs30         = Average shear wave velocity (m/s) in upper 30m of soil
FRV          = Flag for reverse faulting
FNM          = Flag for normal faulting
T            = Period (sec); from 0.01 to 10 sec (use 0 for PGA and -1 for PGV)

Output Variables
Y             = Median V/H ratio of PGA, PGV, or 5% damped SA
sig_lnY       = Total standard deviation of log of V/H ratio
sig_0         = Intra-event (within-event) standard deviation of log of V/H ratio
tau_0         = Inter-event (between-event) standard deviation of log of V/H ratio
%}

%% Specify periods used in GMPM development
T_GA2011 = [0.01	0.02	0.029	0.04	0.05	0.075	0.1	0.15	0.2	0.26	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10];

%% Execute GMPE for desired T
Y_vec = zeros(size(T_vec));
sig_lnY_vec = zeros(size(T_vec));
sig_0_vec = zeros(size(T_vec));
tau_0_vec = zeros(size(T_vec));
for ii = 1:length(T_vec)
    T = T_vec(ii);
    % Check input period
    if ismember(T,[0 -1])
        % Return output for PGA (g), PGV (cm/sec)
        [Y, sig_lnY, sig_0, tau_0] = GA_2011_fixedT(M,Rrup,Vs30,FRV,FNM,T);
    elseif T>=min(T_GA2011) && T<=max(T_GA2011)
        % Determine neighboring periods
        T_lo = T_GA2011(find(T_GA2011<=T,1,'last'));
        T_up = T_GA2011(find(T_GA2011>=T,1,'first'));
        if T_lo==T_up
            % Input T is an element of T_GA2011
            [Y, sig_lnY, sig_0, tau_0] = GA_2011_fixedT(M,Rrup,Vs30,FRV,FNM,T);
        else
            % Determine output for neighboring periods
            [Y_lo, sig_lnY_lo, sig_0_lo, tau_0_lo] = GA_2011_fixedT(M,Rrup,Vs30,FRV,FNM,T_lo);
            [Y_up, sig_lnY_up, sig_0_up, tau_0_up] = GA_2011_fixedT(M,Rrup,Vs30,FRV,FNM,T_up);
            % Linearly interpolate on log scale
            Y = exp(interp1(log([T_lo T_up]),log([Y_lo Y_up]),log(T)));
            sig_lnY = interp1(log([T_lo T_up]),[sig_lnY_lo sig_lnY_up],log(T));
            sig_0 = interp1(log([T_lo T_up]),[sig_0_lo sig_0_up],log(T));
            tau_0 = interp1(log([T_lo T_up]),[tau_0_lo tau_0_up],log(T));
        end
    else
        disp(['Input period is ' num2str(T,4)]);
        error('Invalid input period T');
    end
    Y_vec(ii) = Y;
    sig_lnY_vec(ii) = sig_lnY;
    sig_0_vec(ii) = sig_0;
    tau_0_vec(ii) = tau_0;
end
end



function [Y, sig_lnY, sig_0, tau_0] = GA_2011_fixedT(M,Rrup,Vs30,FRV,FNM,T_fixed)
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

%% Estimate median PGA for rock
Vs30_rock = 1100;
FVS30_rock = 0;
T_rock = 0;
PGAhat1100 = gmpe_as_2008_mod_wCD(M, Vs30_rock, T_rock, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30_rock);
% PGAhat1100 = 0.1; % Used only for reproducing Figs 9-10 in EQS paper


%% Load previously saved regression output
% Fixed inputs
c1 = 6.75;
c4 = 10;
a3 = 0.0147;
a4 = 0.0334;  
a5 = -0.034;
n = 1.18;
c = 1.88;

% Period-dependent inputs
regCoef = [0 865.1 -1.186 0.140 -0.160 -0.105 0.000 0.003 -1.230 0.422 0.333 0.213 0.161
-1 400 -1.955 -1.200 0.090 0.050 0.150 0.022 -1.847 0.373 0.369 0.234 0.170
0.010 865.1 -1.186 0.140 -0.160 -0.105 0.000 0.003 -1.230 0.450 0.330 0.230 0.150
0.020 865.1 -1.219 0.140 -0.160 -0.105 0.000 0.003 -1.268 0.450 0.330 0.230 0.150
0.029 898.6 -1.269 0.335 -0.185 -0.140 0.000 0.003 -1.366 0.450 0.330 0.230 0.150
0.040 994.5 -1.308 0.562 -0.238 -0.160 0.000 0.003 -1.457 0.450 0.341 0.230 0.150
0.050 1053.5 -1.346 0.720 -0.275 -0.136 0.000 -0.001 -1.533 0.450 0.351 0.230 0.150
0.075 1085.7 -1.471 0.552 -0.240 -0.019 0.000 -0.007 -1.706 0.450 0.370 0.230 0.150
0.10 1032.5 -1.624 0.214 -0.169 0.000 0.017 -0.010 -1.831 0.450 0.384 0.230 0.150
0.15 877.6 -1.931 -0.262 -0.069 0.000 0.040 -0.008 -2.114 0.450 0.403 0.230 0.150
0.20 748.2 -2.188 -0.600 0.002 0.000 0.057 -0.003 -2.362 0.450 0.416 0.230 0.150
0.26 639 -2.412 -0.769 0.023 0.000 0.072 0.001 -2.527 0.450 0.429 0.230 0.150
0.30 587.1 -2.518 -0.861 0.034 0.000 0.080 0.006 -2.598 0.450 0.436 0.230 0.150
0.40 503 -2.657 -1.045 0.057 0.000 0.097 0.015 -2.685 0.450 0.449 0.230 0.150
0.50 456.6 -2.669 -1.189 0.075 0.000 0.110 0.022 -2.657 0.450 0.460 0.230 0.150
0.75 410.5 -2.401 -1.250 0.090 0.000 0.133 0.022 -2.265 0.450 0.479 0.237 0.150
1.0 400 -1.955 -1.209 0.090 0.000 0.150 0.022 -1.685 0.450 0.492 0.266 0.150
1.5 400 -1.025 -1.152 0.090 0.029 0.150 0.022 -0.570 0.450 0.511 0.307 0.150
2.0 400 -0.299 -1.111 0.090 0.050 0.150 0.022 0.250 0.532 0.520 0.337 0.150
3.0 400 0.000 -1.054 0.090 0.079 0.150 0.022 0.460 0.648 0.520 0.378 0.213
4.0 400 0.000 -1.014 0.090 0.100 0.150 0.022 0.460 0.700 0.520 0.407 0.258
5.0 400 0.000 -1.000 0.090 0.100 0.150 0.022 0.460 0.700 0.520 0.430 0.292
7.5 400 0.000 -1.000 0.090 0.100 0.150 0.022 0.460 0.700 0.520 0.471 0.355
10.0 400 0.000 -1.000 0.090 0.100 0.150 0.022 0.460 0.700 0.520 0.500 0.400];

periodList = regCoef(:,1);
% % T_GA2011 = periodList(3:end,1);
V_LIN = regCoef(:,2);
b = regCoef(:,3);
a1 = regCoef(:,4);
a2 = regCoef(:,5);
a6 = regCoef(:,6);
a7 = regCoef(:,7);
a8 = regCoef(:,8);
a10 = regCoef(:,9);
s1 = regCoef(:,10);
s2 = regCoef(:,11);
s3 = regCoef(:,12);
s4 = regCoef(:,13);

%% Extract coefficients for given vibration period
ip = find(periodList == T_fixed);

%% Magnitude term
R = sqrt(Rrup^2 + c4^2);
if M<=c1
    f1 = a1(ip) + a4*(M-c1) + a8(ip)*(8.5-M)^2 + (a2(ip)+a3*(M-c1)) * log(R);
else
    f1 = a1(ip) + a5*(M-c1) + a8(ip)*(8.5-M)^2 + (a2(ip)+a3*(M-c1)) * log(R);
end

%% Site amplification term
% Define V1
if T_fixed==-1 % If PGV
    V1 = 862;
elseif T_fixed<=0.5
    V1 = 1500;
elseif 0.5<T_fixed && T_fixed<=1
    V1 = exp( 8-0.795*log(T_fixed/0.21) );
elseif 1<T_fixed && T_fixed<2
    V1 = exp( 6.76-0.297*log(T_fixed) );
elseif T_fixed>=2
    V1 = 700;
else
    error('Invalid input period T entered');
end

% Define Vs30star
if Vs30<V1
    Vs30star = Vs30;
else
    Vs30star = V1;
end

% Determine f5 term
if Vs30star<V_LIN(ip)
    f5 = a10(ip)*log(Vs30star/V_LIN(ip)) - (  -b(ip)*log(PGAhat1100+c) + b(ip)*log( PGAhat1100+c*(Vs30star/V_LIN(ip))^n )  );
else
    f5 = a10(ip)*log(Vs30star/V_LIN(ip)) - (  (b(ip)*n)*log(Vs30star/V_LIN(ip))  );
end

%% Median prediction
lnYhat = f1 + a6(ip)*FRV + a7(ip)*FNM + f5;
Y = exp(lnYhat);

%% Aleatory variability
% Intra-event variability
if M<5
    sig_0 = s1(ip);
elseif M>7
    sig_0 = s2(ip);
else
    sig_0 = s1(ip) + (s2(ip)-s1(ip))/2*(M-5);
end

% Inter-event variability
if M<5
    tau_0 = s3(ip);
elseif M>7
    tau_0 = s4(ip);
else
    tau_0 = s3(ip) + (s4(ip)-s3(ip))/2*(M-5);
end

% Total variability
sig_lnY = sqrt( sig_0^2 + tau_0^2 );

end
