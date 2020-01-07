function [Y, sig_lnY, tau_lnY, phi_lnY] = gmpeV_bc_2016(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,T,W,Ztor,Z2p5,Zhyp)
% By N. Simon Kwong, nealsimonkwong@berkeley.edu
% 
% Bozorgnia and Campbell 2016 GMPM for spectral accelerations of vertical
% components of ground motion. Paper citation:
% 
% Bozorgnia, Y. and Campbell, K. W. (2016). Vertical Ground Motion Model 
% for PGA, PGV, and Linear Response Spectra Using the NGA-West2 Database.
% Earthquake Spectra, 32(2), 979-1004.
% 
% Required inputs:
% M       = Moment magnitude
% Rrup    = Closest distance to rupture plane (km)
% Rjb     = Closest distance to surface projection of rupture plane (km)
% Rx      = Horizontal distance from surface projection of top edge of 
%         rupture plane to site, measured perpendicular to average strike, 
%         and is negative for footwall but positive for hanging wall (km)
% FRV     = Indicator variable representing reverse and reverse-oblique
%         faulting (=1 when rake is within 30 to 150 deg)
% FNM     = Indicator variable representing normal and normal-oblique
%         faulting (=1 when rake is within -30 to -150 deg)
% dip     = Average dip angle of rupture plane (deg)
% Vs30    = Time-averaged shear wave velocity in top 30m of site (m/s)
% region  = 0 for Global
%         = 1 for California
%         = 2 for Japan and Italy
%         = 3 for eastern China
% Sj      = Indicator variable representing Japan's site effects
% T       = Vibration period of interest (sec); =0 for PGA and =-1 for PGV
% 
% Optional inputs (default values can be calculated):
% W       = Down-dip width of rupture plane (km)
% Ztor    = Depth to top of rupture plane (km)
% Z2p5    = Depth to 2.5km/s shear wave velocity, or sediment depth (km)
% Zhyp    = Hypocentral depth of earthquake measured from sea level (km)
% 
% Outputs:
% Y       = Median V component of: 5%-damped spectral acceleration (g) or PGA
%         (g) when T=0 or PGV (cm/s) when T=-1
% sig_lnY = Total aleatory standard deviation
% tau_lnY = Between-event standard deviation
% phi_lnY = Within-event standard deviation

%% If desired, can check inputs for validity here (e.g., Rrup>0)

%% Calculate default values for input vars (see Excel file)
if nargin<12
    Zbot = 15; % Depth to bottom of seismogenic crust (km)
    Zbor = 15; % Depth to bottom of rupture plane (km)
    Ztor = FRV*max([2.704-1.226*max([M-5.849 0]) 0])^2 + (1-FRV)*max([2.673-1.136*max([M-4.97 0]) 0])^2 ;
    if Ztor>20; display('Warning: Ztor exceeds 20km.'); end
    W = min([sqrt(10^((M-4.07)/0.98)) (Zbot-Ztor)/sind(dip)]);
    Z2p5 = 0.6068;
    temp1 = (M<6.75)*(-4.317+0.984*M) + (M>=6.75)*2.325 + (dip<40)*(0.0445*(dip-40));
    temp2 = log(0.9*(Zbor-Ztor));
    Zhyp = max([Ztor+exp(min([temp1 temp2])) 5]);
    if Zhyp>20; display('Warning: Zhyp exceeds 20km.'); end
end

%% Specify periods used in GMPM development
T_BC2016 = [0.010 0.020 0.030 0.050 0.075 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0 0 -1]; % in units of sec except =0 for PGA and =-1 for PGV
T_BC2016sub = T_BC2016(1:(end-2)); % Remove periods for PGA and PGV

%% Execute GMPM for input periods
nTinput = length(T);
Y = zeros(size(T));
sig_lnY = zeros(size(T));
tau_lnY = zeros(size(T));
phi_lnY = zeros(size(T));
for ii=1:nTinput
    Tcurr = T(ii);
    if ismember(Tcurr,T_BC2016) % Pre-defined period requested
        % Report prediction for pre-defined period
        ip = find(T_BC2016 == Tcurr);  
        [Y(ii), sig_lnY(ii), tau_lnY(ii), phi_lnY(ii)] = BC2016vert_sub(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,ip,W,Ztor,Z2p5,Zhyp);
    elseif Tcurr>=min(T_BC2016sub)&& Tcurr<=max(T_BC2016sub) % Assume Tcurr is between 0.01 and 10 sec
        % Determine neighboring periods        
        ip_lo = find(T_BC2016sub<=Tcurr,1,'last');
        ip_hi = find(T_BC2016sub>=Tcurr,1,'first');
        T_lo = T_BC2016sub(ip_lo);
        T_hi = T_BC2016sub(ip_hi);
        
        % Determine output for neighboring periods
        [Y_lo, sig_lnY_lo, tau_lnY_lo, phi_lnY_lo] = BC2016vert_sub(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,ip_lo,W,Ztor,Z2p5,Zhyp);
        [Y_hi, sig_lnY_hi, tau_lnY_hi, phi_lnY_hi] = BC2016vert_sub(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,ip_hi,W,Ztor,Z2p5,Zhyp);
        
        % Linearly interpolate
        Y(ii) = exp(interp1(log([T_lo T_hi]),log([Y_lo Y_hi]),log(Tcurr))); % Log-log scale
        sig_lnY(ii) = interp1(log([T_lo T_hi]),[sig_lnY_lo sig_lnY_hi],log(Tcurr)); % Semilog scale         
        tau_lnY(ii) = interp1(log([T_lo T_hi]),[tau_lnY_lo tau_lnY_hi],log(Tcurr)); % Semilog scale         
        phi_lnY(ii) = interp1(log([T_lo T_hi]),[phi_lnY_lo phi_lnY_hi],log(Tcurr)); % Semilog scale         
    else
        display(['Current input period is ' num2str(Tcurr,4)]);
        error('Invalid input period');
    end
end
end


function [Y, sig_lnY, tau_lnY, phi_lnY] = BC2016vert_sub(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,ip,W,Ztor,Z2p5,Zhyp)
%% Regression output
% Col headers: c0	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10	c11	c12	c13	c14	c15	c16	c17	c18	c19	c20
cVals = [...
-4.674	0.977	0.533	-1.485	-0.445	-2.665	0.214	7.136	0   	-0.229	0.759	-0.354	1.015	0.372	-0.1193	-0.094	0	0.1026	0.0452	0.00784	-0.0053
-4.548	0.976	0.549	-1.488	-0.453	-2.699	0.215	6.936	0   	-0.270	0.768	-0.344	0.950	0.400	-0.1454	-0.081	0	0.1059	0.0427	0.00786	-0.0052
-4.050	0.931	0.628	-1.494	-0.464	-2.772	0.216	7.235	0   	-0.315	0.766	-0.297	1.056	0.394	-0.1957	-0.091	0	0.1175	0.0410	0.00815	-0.0052
-3.435	0.887	0.674	-1.388	-0.552	-2.760	0.202	8.334	0   	-0.329	0.764	-0.363	1.316	0.422	-0.1870	-0.290	0	0.1238	0.0408	0.00783	-0.0062
-3.435	0.902	0.726	-1.469	-0.543	-2.575	0.177	8.761	0   	-0.290	0.795	-0.427	1.758	0.336	-0.0950	-0.261	0	0.1088	0.0516	0.00726	-0.0072
-3.930	0.993	0.698	-1.572	-0.470	-2.461	0.166	9.049	0   	-0.203	0.842	-0.429	1.411	0.314	-0.0999	-0.091	0	0.0918	0.0559	0.00644	-0.0072
-5.505	1.267	0.510	-1.669	-0.452	-2.349	0.164	8.633	0   	-0.203	0.736	-0.421	1.227	0.289	0.0017	-0.092	0	0.0720	0.0447	0.00745	-0.0066
-6.280	1.366	0.447	-1.750	-0.435	-2.335	0.175	8.742	0   	-0.203	0.801	-0.429	0.987	0.290	0.0402	-0.081	0	0.0602	0.0485	0.00789	-0.0056
-6.789	1.458	0.274	-1.711	-0.410	-2.332	0.183	8.400	0       -0.203	0.715	-0.438	0.577	0.303	0.0468	0.011	0	0.0500	0.0416	0.00629	-0.0049
-7.400	1.528	0.193	-1.770	-0.305	-2.297	0.190	7.643	0       -0.203	0.708	-0.421	0.279	0.336	0.0255	0.092	0	0.0382	0.0438	0.00524	-0.0046
-8.750	1.739	-0.020	-1.594	-0.446	-2.219	0.185	7.059	0       -0.203	0.683	-0.401	0.358	0.358	0.0606	0.122	0	0.0264	0.0307	0.00522	-0.0037
-9.740	1.872	-0.121	-1.577	-0.489	-2.205	0.191	6.375	0       -0.203	0.704	-0.417	0.229	0.432	0.0904	0.287	0	0.0163	0.0287	0.00539	-0.0031
-11.050	2.021	-0.042	-1.757	-0.530	-2.143	0.188	5.166	0.016	-0.203	0.602	-0.490	0.574	0.459	0.1776	0.292	0	-0.0016	0.0277	0.00501	-0.0021
-12.184	2.180	-0.069	-1.707	-0.624	-2.092	0.176	5.642	0.032	-0.115	0.394	-0.539	0.980	0.442	0.2389	0.316	0	-0.0072	0.0277	0.00506	-0.0012
-13.451	2.270	0.047	-1.621	-0.686	-1.913	0.144	5.963	0.128	-0.005	0.328	-0.611	0.819	0.520	0.2758	0.450	0	-0.0262	0.0293	0.00353	-0.0004
-13.700	2.271	0.149	-1.512	-0.840	-1.882	0.126	7.584	0.255	0.120	0.112	-0.630	0.044	0.566	0.3051	0.424	0	-0.0408	0.0221	0.00220     0
-13.900	2.150	0.368	-1.315	-0.890	-1.789	0.105	8.645	0.284	0.170	0.011	-0.562	-0.396	0.562	0.3482	0.300	0	-0.0512	0.0321	-0.00137	0
-14.594	2.132	0.726	-1.506	-0.885	-1.781	0.100	10.204	0.26112	0.17	0.000	-0.537	0.001	0.515	0.3527	0.257	0	-0.0567	0.0225	0.00053     0
-15.634	2.116	1.027	-1.721	-0.878	-1.690	0.098	8.386	0.28229	0.17747	0       -0.442	-0.592	0.511	0.3044	0.170	0	-0.0429	0.0237	0.00233     0
-17.129	2.223	0.169	-0.756	-1.077	-1.721	0.125	5.779	0.38692	0.38278	0       -0.343	-1.138	0.575	0.1679	0.219	0	-0.0308	0.0171	-0.00298	0
-17.657	2.132	0.367	-0.800	-1.282	-1.948	0.163	4.135	0.32216	0.33417	0       -0.199	-0.325	0.324	0.1686	0.127	0	0.0067	-0.0017	0.00092     0
-4.729	0.984	0.537	-1.499	-0.443	-2.666	0.214	7.166	0       -0.230	0.759	-0.356	1.019	0.373	-0.1172	-0.097	0	0.1020	0.0442	0.00784	-0.0053
-3.860	1.510	0.270	-1.299	-0.379	-2.383	0.196	6.274	0.111	-0.128	0.140	-0.395	0.338	0.407	-0.0016	0.382	0	0.0581	0.0294	0.00761	-0.0019];
nTfixed = size(cVals,1);

% Additional coefficients
deltaC_20_CA = zeros(nTfixed,1);
deltaC_20_JI = [-0.0018	-0.0018	-0.0020	-0.0026	-0.0021	-0.0018	-0.0018	-0.0022	-0.0025	-0.0027	-0.0024	-0.0025	-0.0025	-0.0023	-0.0013	-0.0004	0	0	0	0	0	-0.0018	0.0005]';
deltaC_20_CH = [0.0039	0.0036	0.0033	0.0039	0.0048	0.0050	0.0048	0.0041	0.0034	0.0031	0.0024	0.0021	0.0020	0.0012	0.0004	0	0	0	0	0	0	0.0039	0.0019]';
a2 = [0.168	0.166	0.167	0.173	0.198	0.174	0.198	0.204	0.185	0.164	0.160	0.184	0.216	0.596	0.596	0.596	0.596	0.596	0.596	0.596	0.596	0.167	0.596]';
 
% Col headers: h1	h2	h3	h4	h5	h6
hVals = [...
0.242	1.471	-0.714	1.000	-0.336	-0.270
0.244	1.467	-0.711	1.000	-0.339	-0.263
0.246	1.467	-0.713	1.000	-0.338	-0.259
0.251	1.449	-0.701	1.000	-0.338	-0.263
0.260	1.435	-0.695	1.000	-0.347	-0.219
0.259	1.449	-0.708	1.000	-0.391	-0.201
0.254	1.461	-0.715	1.000	-0.449	-0.099
0.237	1.484	-0.721	1.000	-0.393	-0.198
0.206	1.581	-0.787	1.000	-0.339	-0.210
0.210	1.586	-0.795	1.000	-0.447	-0.121
0.226	1.544	-0.770	1.000	-0.525	-0.086
0.217	1.554	-0.770	1.000	-0.407	-0.281
0.154	1.626	-0.780	1.000	-0.371	-0.285
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.117	1.616	-0.733	1.000	-0.128	-0.756
0.241	1.474	-0.715	1.000	-0.337	-0.270
0.117	1.616	-0.733	1.000	-0.128	-0.756];

% Col headers: k1	k2	k3
kVals = [...
865	0.000	0.000
865	0.000	0.000
908	0.000	0.000
1054	0.000	0.000
1086	0.000	0.000
1032	0.000	0.000
878	0.000	0.000
748	0.000	0.000
654	0.000	0.000
587	0.000	0.000
503	0.000	0.000
457	0.000	0.000
410	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
400	0.000	0.000
865	0.000	0.000
400	0.000	0.000];    
    
% Standard deviations (between or inter event)
tau_1 = [0.462	0.474	0.529	0.576	0.523	0.461	0.391	0.363	0.355	0.355	0.36	0.376	0.416	0.472	0.507	0.539	0.515	0.553	0.578	0.6	0.495	0.461	0.334]';
tau_2 = [0.345	0.375	0.416	0.468	0.427	0.39	0.343	0.308	0.288	0.265	0.28	0.284	0.322	0.311	0.329	0.345	0.335	0.331	0.294	0.379	0.442	0.347	0.24]';

% Standard deviations (within or intra event)
phi_1 = [0.695	0.7	0.722	0.751	0.74	0.723	0.731	0.701	0.687	0.668	0.628	0.606	0.568	0.536	0.511	0.507	0.474	0.466	0.43	0.386	0.395	0.694	0.608]';
phi_2 = [0.494	0.508	0.536	0.584	0.578	0.57	0.536	0.51	0.507	0.514	0.521	0.526	0.536	0.55	0.559	0.571	0.557	0.566	0.568	0.527	0.481	0.493	0.442]';

%% Magnitude term
if M<=4.5
    f_mag = cVals(ip,1) + cVals(ip,2)*M;
elseif M>4.5 && M<=5.5
    f_mag = cVals(ip,1) + cVals(ip,2)*M + cVals(ip,3)*(M-4.5);
elseif M>5.5 && M<=6.5
    f_mag = cVals(ip,1) + cVals(ip,2)*M + cVals(ip,3)*(M-4.5) + cVals(ip,4)*(M-5.5);    
else
    f_mag = cVals(ip,1) + cVals(ip,2)*M + cVals(ip,3)*(M-4.5) + cVals(ip,4)*(M-5.5) + cVals(ip,5)*(M-6.5);
end

%% Geometric attenuation term
f_dist = (cVals(ip,6) + cVals(ip,7)*M) * log(sqrt( Rrup^2 + cVals(ip,8)^2 ));

%% Style-of-faulting term
f_fltF = cVals(ip,9)*FRV + cVals(ip,10)*FNM;
if M<=4.5
    f_fltM = 0;
elseif M>4.5 && M<=5.5
    f_fltM = M-4.5;
else
    f_fltM = 1;
end
f_flt = f_fltF * f_fltM;

%% Hanging-wall term
% Preliminary calcs
R1 = W*cosd(dip);
R2 = 62*M - 350;
f1 = hVals(ip,1) + hVals(ip,2)*(Rx/R1) + hVals(ip,3)*(Rx/R1)^2;
f2 = hVals(ip,4) + hVals(ip,5)*((Rx-R1)/(R2-R1)) + hVals(ip,6)*((Rx-R1)/(R2-R1))^2;

% Rx intermediate term
if Rx<0
    f_hngRx = 0;
elseif Rx>=0 && Rx<R1
    f_hngRx = f1;
else
    f_hngRx = max([f2 0]);
end

% Rrup intermediate term
if Rrup==0
    f_hngRrup = 1;
else
    f_hngRrup = (Rrup-Rjb)/Rrup;
end

% M intermediate term
if M<=5.5
    f_hngM = 0;
elseif M>5.5 && M<=6.5
    f_hngM = (M-5.5)*(1 + a2(ip)*(M-6.5));
else
    f_hngM = 1 + a2(ip)*(M-6.5);
end

% Z intermediate term
if Ztor<=16.66
    f_hngZ = 1-0.06*Ztor;
else
    f_hngZ = 0;
end

% dip intermediate term
f_hngDip = (90-dip)/45;

% Combine terms
f_hng = cVals(ip,11) * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngDip;

%% Shallow site response term
f_siteG = cVals(ip,12)*log(Vs30/kVals(ip,1)); % Since k2=0 for all periods, no need to consider A_1100 in defining f_site
if Vs30<=200
    f_siteJ = cVals(ip,13)*(log(Vs30/kVals(ip,1)) - log(200/kVals(ip,1)));
else
    f_siteJ = cVals(ip,14)*log(Vs30/kVals(ip,1));
end
f_site = f_siteG + Sj*f_siteJ;

%% Basin response term
if Z2p5<=1
    f_sed = (cVals(ip,15) + cVals(ip,16)*Sj) * (Z2p5 - 1);
elseif Z2p5>1 && Z2p5<=3
    f_sed = 0;
else
    f_sed = cVals(ip,17)*kVals(ip,3)*exp(-0.75)*( 1 - exp( -0.25*(Z2p5-3) ) ); % Even though k3=0 for all periods, kept this for future
end

%% Hypocentral depth term
% Hypocentral depth effect
if Zhyp<=7
    f_hypH = 0;
elseif Zhyp>7 && Zhyp<=20
    f_hypH = Zhyp - 7;
else
    f_hypH = 13;
end
% Magnitude effect
if M<=5.5
    f_hypM = cVals(ip,18);
elseif M>5.5 && M<=6.5
    f_hypM = cVals(ip,18) + (cVals(ip,19)-cVals(ip,18))*(M-5.5);
else
    f_hypM = cVals(ip,19);
end
% Combine terms
f_hyp = f_hypH * f_hypM;

%% Fault dip term
if M<=4.5
    f_dip = cVals(ip,20)*dip;
elseif M>4.5 && M<=5.5
    f_dip = cVals(ip,20)*(5.5-M)*dip;
else
    f_dip = 0;
end

%% Anelastic attenuation term
% Identify regional effect
if region==2 % Japan or Italy
    deltaC_20 = deltaC_20_JI(ip);
elseif region==3 % Eastern China
    deltaC_20 = deltaC_20_CH(ip);
else
    deltaC_20 = deltaC_20_CA(ip); % Global or CA
end
% Determine regional anelastic attenuation
if Rrup>80
    f_atn = (cVals(ip,21) + deltaC_20)*(Rrup - 80);
else
    f_atn = 0;
end

%% Median prediction
lnY = f_mag + f_dist + f_flt + f_hng + f_site + f_sed + f_hyp + f_dip + f_atn;
Y = exp(lnY);

%% Aleatory variability
% Between/inter-event variability
if M<=4.5
    tau_lnY = tau_1(ip);
elseif M>=5.5
    tau_lnY = tau_2(ip);
else
    tau_lnY = tau_1(ip) + (tau_2(ip)-tau_1(ip))/(5.5-4.5)*(M-4.5);
end

% Within/intra-event variability
if M<=4.5
    phi_lnY = phi_1(ip);
elseif M>=5.5
    phi_lnY = phi_2(ip);
else
    phi_lnY = phi_1(ip) + (phi_2(ip)-phi_1(ip))/(5.5-4.5)*(M-4.5);
end

% Total variability
sig_lnY = sqrt( tau_lnY^2 + phi_lnY^2 );
end