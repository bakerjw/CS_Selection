% coded by Yoshifumi Yamamoto, 5/6/10
%               Stanford University
%               yama4423@stanford.edu
% based on AS1997 by Jack W. Baker, 5/5/05
%               Stanford University
%               bakerjw@stanford.edu
%
%   updated 2011/01/19
%          Hanging Wall term
%          from the Errata for “AS NGA” model (http://peer.berkeley.edu/products/abrahamson-silva_nga_report_files/AS08_NGA_errata.pdf)
%          This version includes constant displacement effect
%   updated 2010/05/20
%   updated 2009/05/05
%
% Summary of the Abrahamson & Silva Ground-Motion Relations
% Norman Abrahamson, and Walter Silva
% Earthquake Spectra, Volume 24, No.1, pages 67-97, February 2008
%
% This script has been modified to correct an error based on the openSHA
% about constant displacement model(equation 22) and standard
% deviation(equation 24 and 26).


function [Sa, tsigma, period1, pga_rock, sigma, tau] = gmpe_as_2008_mod_wCD(M, Vs30, T, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30)
% Initialize struct V = constants
V = struct('period',0,'lin',0,...
    'a1',0,'a2',0,'a3',0,'a4',0,'a5',0,'a8',0,'a10',0,'a12',0,'a13',0,'a14',0,'a15',0,'a16',0,'a18',0,...
    'b',0,'c',0,'c1',0,'c2',0,'c4',0,'n',0,'s1',0,'s2',0,'s3',0,'s4',0,...
    'ro',0,'v1',0,'e2',0,'a22',0,'Vs30s',0);  

% Get info for rock conditions
Vrock = get_abrahamson_silva_constants(1,1100,FVS30,V); % Update constants for rock
pga_rock = exp(calc_val(M, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, 0, 1100, Vrock));
Vupdated = get_abrahamson_silva_constants(1,Vs30,FVS30,V);
[~, ~, ~, pga_sigmaB, pga_tauB] = abrahamson_silva_sigma(M, pga_rock, Vs30,0, 0, Vupdated);

% Get data for Td
Td = 10^( -1.25 + 0.3*M );
[SaTd] = AS_2008_nga_sub(M, 1100, Td, Rrup, Rjb, Rx, dip, Ztor, 0, W, FRV, FNM, FAS, FHW, FVS30, 0, 0, 1, V, pga_rock, pga_sigmaB, pga_tauB);

% Actual GMPM
[Sa, tsigma, period1, pga_rock, sigma, tau] = AS_2008_nga_sub(M, Vs30, T, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30, Td, SaTd, 0, V, pga_rock, pga_sigmaB, pga_tauB);


function [Sa, tsigma, period1, pga_rock, sigma, tau] = AS_2008_nga_sub(M, Vs30, T, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30, Td, SaTd, irock, V, pga_rock, pga_sigmaB, pga_tauB)
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   M               = moment magnitude
%   T               = period of vibration
%                 Use 1000 for output the array of Sa with period
%   Rrup            = closest distance to fault rupture (rupture distance)(km)
%   Rjb             = Joyner-Boore distance(km)
%   Rx              = Horizontal distance (km) from top edge of rupture
%   dip             = dip angle of fault in degree
%   Ztor            = Depth-to-top of rupture (km)
%   Z10             = Depth to Vs=1.0km/s at the site(m)
%   W               = Down-dip rupture width (km)
%   FRV             = Flag for reverse faulting earthquakes
%                   = 1 for reverse and reverse/oblique earthquakes defined
%                           by rake angles between 30 and 150 degrees
%                   = 0 otherwise
% 
%   FNM             = Flag for normal faulting earthquakes
%                   = 1 for normal earthquakes defined by rake angles
%                           between -60 and -120 degrees
%                   = 0 otherwise
% 
%   FAS             = Flag for aftershocks
%                   = 1 for aftershocks
%                   = 0 for mainshocks, foreshocks, and swarms
% 
%   FHW             = 1 for Hanging Wall sites
%                   = 0 otherwise
% 
%   FVS30           = 1 for estimated Vs30
%                   = 0 for measured  Vs30
% 
%   ircok           = 1 rock
%                   = 0 not rock
%
% OUTPUT   
%
%   Sa              = median spectral acceleration prediction
%   tsigma           = logarithmic standard deviation of spectral acceleration
%                     prediction FOR AN ARBITRARY OR AVERAGE COMPONENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the given period T, get the index for the constants
period = [    0,      -1,    0.01,    0.02,    0.03,    0.04,    0.05,   0.075,     0.1,    0.15,     0.2,    0.25,     0.3,     0.4,     0.5,    0.75,       1,     1.5,       2,       3,       4,       5,     7.5,      10];
    
nT=length(T); % Length of input vector of vibration periods T
iflg=0;
if (nT==1 && T==1000); % If input period is scalar and if input period is 1000 sec
    iflg=1;
    nperi=length(period);
    Sa=zeros(1,nperi-2);
    tsigma=zeros(1,nperi-2);
    period1=period(3:end);
    
    for index=3:1:nperi;
        % get constants for the given index value
        V = get_abrahamson_silva_constants(index,Vs30,FVS30,V);
        if period(index)<=Td || SaTd==0;
            Sa(index-2) = exp(calc_val(M, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, pga_rock, Vs30, V) + (1-irock)* f_10(Z10, Vs30, V));
        else
            Sa(index-2) = exp(calc_val2(SaTd,Td, period(index), Z10, pga_rock, Vs30, V));
        end;
        [tsigma(index-2), sigma, tau, sigmaB, tauB] = abrahamson_silva_sigma(M, pga_rock, Vs30, pga_sigmaB, pga_tauB, V);
    end;
end;

if(iflg==0);
    Sa=zeros(1,nT);
    tsigma=zeros(1,nT);
    sigma = zeros(1,nT);
    tau = zeros(1,nT);
    period1=T;
    for it=1:1:nT;
        Teach=T(it); % Current value of period within input period vector
        if Teach>period(end); Teach = period(end); end;
    % interpolate between periods if neccesary  
        if all(Teach ~= period) % If current period not one of those considered in GMPM; use logical indexing 
            T_low = max(period(period<Teach)); % Avoid using find()
            T_hi = min(period(period>Teach));
                 
            % Reduce calls to AS_2008_nga_sub
            T_sub = [T_low T_hi];
            [sa_sub, tsigma_sub, ~, ~, sigma_sub, tau_sub] = AS_2008_nga_sub(M, Vs30, T_sub, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, FVS30, Td, SaTd, irock, V, pga_rock, pga_sigmaB, pga_tauB);
            sa_low = sa_sub(1,1); sa_hi = sa_sub(1,2);
            tsigma_low = tsigma_sub(1,1); tsigma_hi = tsigma_sub(1,2);
            sigma_low = sigma_sub(1,1); sigma_hi = sigma_sub(1,2);
            tau_low = tau_sub(1,1); tau_hi = tau_sub(1,2);
                        
            x = [log(T_low) log(T_hi)];
            Y_sa = [log(sa_low) log(sa_hi)];
            Y_tsigma = [tsigma_low tsigma_hi];
            Y_sigma = [sigma_low sigma_hi];
            Y_tau = [tau_low tau_hi];
            Sa(it) = exp(interp1(x,Y_sa,log(Teach)));
            tsigma(it) = interp1(x,Y_tsigma,log(Teach));
            sigma(it) = interp1(x,Y_sigma,log(Teach));
            tau(it) = interp1(x,Y_tau,log(Teach));
        else
            index = find(abs((period - Teach)) < 0.0001); % Identify the period

            % get constants for the given index value
            V = get_abrahamson_silva_constants(index,Vs30,FVS30,V);
            if period(index)<=Td || SaTd==0;
                Sa(it) = exp(calc_val(M, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, pga_rock, Vs30, V) + (1-irock)* f_10(Z10, Vs30, V));
            else
                Sa(it) = exp(calc_val2(SaTd,Td, period(index), Z10, pga_rock, Vs30, V));
            end;

            [tsigma(it), sigma(it), tau(it), sigmaB, tauB] = abrahamson_silva_sigma(M, pga_rock, Vs30, pga_sigmaB, pga_tauB, V);

        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1] = f_1(M, R, V)
% value of f1
	if (M <= V.c1) 
        f1 = V.a1 + V.a4 * (M - V.c1) + V.a8 * (8.5 - M)^2 + (V.a2 + V.a3 * (M - V.c1)) * log(R);
	else 
        f1 = V.a1 + V.a5 * (M - V.c1) + V.a8 * (8.5 - M)^2 + (V.a2 + V.a3 * (M - V.c1)) * log(R);
	end

function [f4] = f_4(Rjb, Rx, dip, Ztor, M, W, V)
% value of f_4
    if Rjb<30
        T1=1-Rjb/30;
    else
        T1=0;
    end;
    W1=W*cos(dip);
    if Rx<=W1;
        T2=0.5+Rx/(2*W1);
    elseif Rx>W1 || dip==90;
        T2=1;
    end;
    if Rx>=Ztor;
        T3=1;
    else
        T3=Rx/Ztor;
    end;
    if M<=6;
        T4=0;
    elseif M<7;
        T4=M-6;
    else
        T4=1;
    end;
% from paper AS08
%     if dip>=70;
%         T5=1-(dip-70)/20;
%     else
%         T5=1;
%     end;
% from the Errata for “AS NGA” model (http://peer.berkeley.edu/products/abrahamson-silva_nga_report_files/AS08_NGA_errata.pdf)
    if dip>=30;
        T5=1-(dip-30)/60;
    else
        T5=1;
    end;
    
    
    f4 = V.a14*T1*T2*T3*T4*T5;

function [f5] = f_5(pga_rock, Vs30, V)
% value of f_5
if Vs30 < V.lin;
    f5 = V.a10 * log(V.Vs30s/V.lin) - V.b*log(pga_rock+V.c) + V.b*log(pga_rock+V.c*((V.Vs30s/V.lin)^V.n));
else
    f5 = (V.a10 + V.b*V.n) * log(V.Vs30s/V.lin);
end;

function [f6] = f_6(Ztor, V)
% value of f_6
    if Ztor<10
        f6=V.a16*Ztor/10;
    else
        f6=V.a16;
    end;
    
function [f8] = f_8(Rrup, M, V)
% value of f_8
    if M<5.5;
        T6=1;
    elseif M<=6.5;
        T6=0.5*(6.5-M)+0.5;
    else
        T6=0.5;
    end;
    if Rrup<100
        f8=0;
    else
        f8=V.a18*(Rrup-100)*T6;
    end;

function [f10] = f_10(Z10, Vs30, V)
% value of f_10
    if Vs30<180;
        Z10h=exp(6.745);
    elseif Vs30<=500;
        Z10h=exp(6.745-1.35*log(Vs30/180));
    else
        Z10h=exp(5.394-4.48*log(Vs30/500));
    end;
    a211=(V.a10+V.b*V.n)*log(V.Vs30s/min(V.v1,1000));
    a212=log((Z10+V.c2)/(Z10h+V.c2));
    if Vs30>=1000;
        a21=0;
    elseif a211+V.e2*a212<0;
        a21=-a211/a212;
    else
        a21=V.e2;
    end;
    
    f10=a21*a212;
    if Z10>=200
        f10=f10+V.a22*log(Z10/200);
    end;

   
    
function [X] = calc_val(M, Rrup, Rjb, Rx, dip, Ztor, Z10, W, FRV, FNM, FAS, FHW, pga_rock, Vs30, constants)
% calculate predicted value

    R = sqrt(Rrup^2 + constants.c4^2);
	X = f_1(M, R, constants) +  constants.a12*FRV + constants.a13*FNM + constants.a15*FAS ...
        + f_5(pga_rock, Vs30, constants) + FHW*f_4(Rjb, Rx, dip, Ztor, M, W, constants) + f_6(Ztor, constants) ...
        + f_8(Rrup, M, constants);
%         + f_8(Rrup, M, constants) + f_10(Z10, Vs30, constants);

function [X] = calc_val2(SaTd, Td, T, Z10, pga_rock, Vs30, constants)
% calculate predicted value
	X = log((SaTd) * Td^2 / T^2) - f_5(pga_rock, 1100, constants) + f_5(pga_rock, Vs30, constants) + f_10(Z10, Vs30, constants);
  
function [constants] = get_abrahamson_silva_constants(index,Vs30,FVS30,constants)
% get relevant constants

% arrays with values by index
period = [      0,      -1,    0.01,    0.02,    0.03,    0.04,    0.05,   0.075,     0.1,    0.15,     0.2,    0.25,     0.3,     0.4,     0.5,    0.75,       1,     1.5,       2,       3,       4,       5,     7.5,      10];
lin    = [  865.1,   400.0,   865.1,   865.1,   907.8,   994.5,  1053.5,  1085.7,  1032.5,   877.6,   748.2,   654.3,   587.1,   503.0,   456.6,   410.5,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0];
b      = [ -1.186,  -1.955,  -1.186,  -1.219,  -1.273,  -1.308,  -1.346,  -1.471,  -1.624,  -1.931,  -2.188,  -2.381,  -2.518,  -2.657,  -2.669,  -2.401,  -1.955,  -1.025,  -0.299,     0.0,     0.0,     0.0,     0.0,     0.0];
a1     = [  0.804,  5.7578,   0.811,   0.855,   0.962,   1.037,   1.133,   1.375,   1.563,   1.716,   1.687,   1.646,   1.601,   1.511,   1.397,   1.137,   0.915,   0.510,   0.192,  -0.280,  -0.639,  -0.936,  -1.527,  -1.993];
a2     = [-0.9679, -0.9046, -0.9679, -0.9774, -1.0024, -1.0289, -1.0508, -1.0810, -1.0833, -1.0357, -0.9700, -0.9202, -0.8974, -0.8677, -0.8475, -0.8206, -0.8088, -0.7995, -0.7960, -0.7960, -0.7960, -0.7960, -0.7960, -0.7960];
a8     = [-0.0372,   -0.12, -0.0372, -0.0372, -0.0372, -0.0315, -0.0271, -0.0191, -0.0166, -0.0254, -0.0396, -0.0539, -0.0656, -0.0807, -0.0924, -0.1137, -0.1289, -0.1534, -0.1708, -0.1954, -0.2128, -0.2263, -0.2509, -0.2683];
a10    = [ 0.9445,  1.5390,  0.9445,  0.9834,  1.0471,  1.0884,  1.1333,  1.2808,  1.4613,  1.8071,  2.0773,  2.2794,  2.4201,  2.5510,  2.5395,  2.1493,  1.5705,  0.3991, -0.6072, -0.9600, -0.9600, -0.9208, -0.7700, -0.6630];
a12    = [ 0.0000,  0.0800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0181,  0.0309,  0.0409,  0.0491,  0.0619,  0.0719,  0.0800,  0.0800,  0.0800,  0.0800,  0.0800,  0.0800,  0.0800,  0.0800,  0.0800];
a13    = [-0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600];

a14    = [ 1.0800,  0.7000,  1.0800,  1.0800,  1.1331,  1.1708,  1.2000,  1.2000,  1.2000,  1.1683,  1.1274,  1.0956,  1.0697,  1.0288,  0.9971,  0.9395,  0.8985,  0.8409,  0.8000,  0.4793,  0.2518,  0.0754,  0.0000,  0.0000];
a15    = [-0.3500, -0.3900, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3191, -0.2629, -0.2230, -0.1668, -0.1270, -0.0708, -0.0309,  0.0000,  0.0000,  0.0000];
a16    = [ 0.9000,  0.6300,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.9000,  0.8423,  0.7458,  0.5704,  0.4460,  0.2707,  0.1463, -0.0291, -0.1535, -0.2500, -0.2500, -0.2500];
a18    = [-0.0067,  0.0000, -0.0067, -0.0067, -0.0067, -0.0067, -0.0076, -0.0093, -0.0093, -0.0093, -0.0083, -0.0069, -0.0057, -0.0039, -0.0025,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000];


s1e    = [  0.590,   0.590,   0.590,   0.590,   0.605,   0.615,   0.623,   0.630,   0.630,   0.630,   0.630,   0.630,   0.630,   0.630,   0.630,   0.630,   0.630,   0.615,   0.604,   0.589,   0.578,   0.570,   0.611,   0.640];
s2e    = [  0.470,   0.470,   0.470,   0.470,   0.478,   0.483,   0.488,   0.495,   0.501,   0.509,   0.514,   0.518,   0.522,   0.527,   0.532,   0.539,   0.545,   0.552,   0.558,   0.565,   0.570,   0.587,   0.618,   0.640];
s1m    = [  0.576,   0.576,   0.576,   0.576,   0.591,   0.602,   0.610,   0.617,   0.617,   0.616,   0.614,   0.612,   0.611,   0.608,   0.606,   0.602,   0.594,   0.566,   0.544,   0.527,   0.515,   0.510,   0.572,   0.612];
s2m    = [  0.453,   0.453,   0.453,   0.453,   0.461,   0.466,   0.471,   0.479,   0.485,   0.491,   0.495,   0.497,   0.499,   0.501,   0.504,   0.506,   0.503,   0.497,   0.491,   0.500,   0.505,   0.529,   0.579,   0.612];
s3     = [  0.470,   0.420,   0.420,   0.420,   0.462,   0.492,   0.515,   0.550,   0.550,   0.550,   0.520,   0.497,   0.479,   0.449,   0.426,   0.385,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350];
s4     = [  0.300,   0.300,   0.300,   0.300,   0.305,   0.309,   0.312,   0.317,   0.321,   0.326,   0.329,   0.332,   0.335,   0.338,   0.341,   0.346,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350,   0.350];
ro     = [  1.000,   0.740,   1.000,   1.000,   0.991,   0.982,   0.973,   0.952,   0.929,   0.896,   0.874,   0.856,   0.841,   0.818,   0.783,   0.680,   0.607,   0.504,   0.431,   0.328,   0.255,   0.200,   0.200,   0.200];

c1=6.75;
c4=4.5;
a3=0.265;
a4=-0.231;
a5=-0.398;
n=1.18;
c=1.88;
c2=50;

constants.period = period(index);
constants.lin    = lin(index);
constants.a1     = a1(index);
constants.a2     = a2(index);
constants.a3     = a3;
constants.a4     = a4;
constants.a5     = a5;
constants.a8     = a8(index);
constants.a10    = a10(index);
constants.a12    = a12(index);
constants.a13    = a13(index);
constants.a14    = a14(index);
constants.a15    = a15(index);
constants.a16    = a16(index);
constants.a18    = a18(index);
constants.b      = b(index);
constants.c      = c;
constants.c1     = c1;
constants.c2     = c2;
constants.c4     = c4;
constants.n      = n;

if FVS30==1;
    constants.s1    = s1e(index);
    constants.s2    = s2e(index);
else
    constants.s1    = s1m(index);
    constants.s2    = s2m(index);
end;
constants.s3     = s3(index);
constants.s4     = s4(index);
constants.ro     = ro(index);


T=period(index);
if index==2;
    constants.v1 = 862;
elseif T<=0.5;
    constants.v1 = 1500;
elseif T<=1;
    constants.v1 = exp(8.0-0.795*log(T/0.21));
elseif T<2;
    constants.v1 = exp(6.76-0.297*log(T));
else
    constants.v1 = 700;
end;
if T<0.35 || Vs30>1000;
    constants.e2=0;
elseif T<=2;
    constants.e2=-0.25*log(Vs30/1000)*log(T/0.35);
else
    constants.e2=-0.25*log(Vs30/1000)*log(2/0.35);
end;
if T<2;
    constants.a22=0;
else
    constants.a22=0.0625*(T-2);
end;

constants.Vs30s=min(Vs30,constants.v1);

function [tsigma, sigma, tau, sigmaB, tauB] = abrahamson_silva_sigma(M, pga_rock, Vs30, pga_sigmaB, pga_tauB, V)
% calculate the sigma

% use the published coefficients for the geometric mean
if M<5;
    sigma0 = V.s1;
elseif M<=7;
    sigma0 = V.s1 + (V.s2-V.s1)/2 * (M-5);
else
    sigma0 = V.s2;
end;
sigmaAMP=0.3;
sigmaB=sqrt(sigma0^2-sigmaAMP^2);
if M<5;
    tau0 = V.s3;
elseif M<=7;
    tau0 = V.s3 + (V.s4-V.s3)/2 * (M-5);
else
    tau0 = V.s4;
end;
tauB=tau0;

if Vs30>=V.lin;
    term1 = 0;
else
% from openSHA
    term1 = V.b * pga_rock * ( (-1/(pga_rock+V.c)) + (1/(pga_rock + V.c*((Vs30/V.lin)^V.n)))  );
end;
% from openSHA
sigma = sqrt(sigma0^2 + term1^2 * pga_sigmaB^2 + 2*term1 * sigmaB*pga_sigmaB*V.ro);
tau   = sqrt(tau0^2   + term1^2 * pga_tauB^2   + 2*term1 * tauB  *pga_tauB  *V.ro);

tsigma = sqrt(sigma^2 + tau^2);
