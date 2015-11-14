% by Yoshifumi Yamamoto, 11/10/08
% Stanford University
% yama4423@stanford.edu
% fix bug 2009/05/05
%
% Campbell-Bozorgnia attenuation equation, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
% M             = Magnitude
% T             = Period (sec); Use -1 for PGV computation and -10 for PGD
%                 computation
%                 Use 1000 for output the array of Sa with period
% Rrup          = Closest distance coseismic rupture (km)
% Rjb           = Joyner-Boore distance (km)
% Ztor          = Depth to the top of coseismic rupture (km)
% delta         = average dip of the rupture place (degree)
% lambda        = rake angle (degree)
% Vs30          = shear wave velocity averaged over top 30 m (m/s)
% Zvs           = Depth to the 2.5 km/s shear-wave velocity horizon (km)
% arb           = 1 for arbitrary component sigma
%               = 0 for average component sigma
%
% Output Variables
% Sa            = Median spectral acceleration prediction
% sigma         = logarithmic standard deviation of spectral acceleration
%                 prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sa, sigma, period1] = CB_2008_nga (M, T, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb)

% Coefficients
period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10 0 -1 -10];

        % Style of faulting
if Ztor < 1;
    ffltz = Ztor;
else    
    ffltz = 1;
end;
if Rjb == 0;
    fhngr = 1;
else    
    if Ztor < 1;
        fhngr = (max(Rrup,sqrt(Rjb^2+1))-Rjb)/max(Rrup,sqrt(Rjb^2+1));
    else    
        fhngr = (Rrup - Rjb)/Rrup;
    end;
end;

Frv = (lambda > 30 & lambda < 150);
Fnm = (lambda > -150 & lambda < -30);
if M <= 6;
    fhngm = 0;
elseif M < 6.5;
    fhngm = 2 * (M - 6);
else    
    fhngm = 1;
end;

nT=length(T);
iflg=0;
if(nT==1);
    if(T==1000);
        iflg=1;
        %     Computation with original periods
        nperi=length(period)-3;
        Sa=zeros(1,nperi);
        sigma=zeros(1,nperi);
        period1=period(1:nperi);
        for i=1:1:nperi;
            [Sa(i) sigma(i)] = CB_2008_nga_sub (M, i, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, ffltz, fhngr, Frv, Fnm, fhngm);
        end;
    end;
end;
if(iflg==0);
    Sa=zeros(1,nT);
    sigma=zeros(1,nT);
    period1=T;
    for it=1:1:nT;
        Teach=T(it);
        % interpolate between periods if neccesary    
        if (length(find(period == Teach)) == 0)
            T_low = max(period(find(period<Teach)));
            T_hi = min(period(find(period>Teach)));

%             [sa_low sigma_low] = CB_2008_nga (M, T_low, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, ffltz, fhngr, Frv, Fnm, fhngm);
%             [sa_hi sigma_hi] = CB_2008_nga (M, T_hi, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, ffltz, fhngr, Frv, Fnm, fhngm);
            [sa_low sigma_low] = CB_2008_nga (M, T_low, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb);
            [sa_hi sigma_hi] = CB_2008_nga (M, T_hi, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb);

            x = [log(T_low) log(T_hi)];
            Y_sa = [log(sa_low) log(sa_hi)];
            Y_sigma = [sigma_low sigma_hi];
            Sa(it) = exp(interp1(x,Y_sa,log(Teach)));
            sigma(it) = interp1(x,Y_sigma,log(Teach));
        else
            i = find(period == Teach); 
            [Sa(it) sigma(it)] = CB_2008_nga_sub (M, i, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, ffltz, fhngr, Frv, Fnm, fhngm);
        end;
    end;
end;

function [Sa sigma] = CB_2008_nga_sub (M, ip, Rrup, Rjb, Ztor, delta, lambda, Vs30, Zvs, arb, ffltz, fhngr, Frv, Fnm, fhngm)

c0 = [-1.715	-1.68	-1.552	-1.209	-0.657	-0.314	-0.133	-0.486	-0.89	-1.171	-1.466	-2.569	-4.844	-6.406	-8.692	-9.701	-10.556	-11.212	-11.684	-12.505	-13.087	-1.715	0.954	-5.27];
c1 = [0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.656	0.972	1.196	1.513	1.6	1.6	1.6	1.6	1.6	1.6	0.5	0.696	1.6];
c2 = [-0.53	-0.53	-0.53	-0.53	-0.53	-0.53	-0.53	-0.446	-0.362	-0.294	-0.186	-0.304	-0.578	-0.772	-1.046	-0.978	-0.638	-0.316	-0.07	-0.07	-0.07	-0.53	-0.309	-0.07];
c3 = [-0.262	-0.262	-0.262	-0.267	-0.302	-0.324	-0.339	-0.398	-0.458	-0.511	-0.592	-0.536	-0.406	-0.314	-0.185	-0.236	-0.491	-0.77	-0.986	-0.656	-0.422	-0.262	-0.019	0];
c4 = [-2.118	-2.123	-2.145	-2.199	-2.277	-2.318	-2.309	-2.22	-2.146	-2.095	-2.066	-2.041	-2	-2	-2	-2	-2	-2	-2	-2	-2	-2.118	-2.016	-2];
c5 = [0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17];
c6 = [5.6	5.6	5.6	5.74	7.09	8.05	8.79	7.6	6.58	6.04	5.3	4.73	4	4	4	4	4	4	4	4	4	5.6	4	4];
c7 = [0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.255	0.161	0.094	0	0	0	0	0	0.28	0.245	0];
c8 = [-0.12	-0.12	-0.12	-0.12	-0.12	-0.099	-0.048	-0.012	0	0	0	0	0	0	0	0	0	0	0	0	0	-0.12	0	0];
c9 = [0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.371	0.154	0	0	0	0	0.49	0.358	0];
c10 = [1.058	1.102	1.174	1.272	1.438	1.604	1.928	2.194	2.351	2.46	2.587	2.544	2.133	1.571	0.406	-0.456	-0.82	-0.82	-0.82	-0.82	-0.82	1.058	1.694	-0.82];
c11 = [0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.077	0.15	0.253	0.3	0.3	0.3	0.3	0.3	0.3	0.04	0.092	0.3];
c12 = [0.61	0.61	0.61	0.61	0.61	0.61	0.61	0.61	0.70	0.75	0.85	0.883	1	1	1	1	1	1	1	1	1	0.61	1	1];
k1 = [865	865	908	1054	1086	1032	878	748	654	587	503	457	410	400	400	400	400	400	400	400	400	865	400	400];
k2 = [-1.186	-1.219	-1.273	-1.346	-1.471	-1.624	-1.931	-2.188	-2.381	-2.518	-2.657	-2.669	-2.401	-1.955	-1.025	-0.299	0	0	0	0	0	-1.186	-1.955	0];
k3 = [1.839	1.84	1.841	1.843	1.845	1.847	1.852	1.856	1.861	1.865	1.874	1.883	1.906	1.929	1.974	2.019	2.11	2.2	2.291	2.517	2.744	1.839	1.929	2.744];
sigmac = [0.166	0.166	0.165	0.162	0.158	0.17	0.18	0.186	0.191	0.198	0.206	0.208	0.221	0.225	0.222	0.226	0.229	0.237	0.237	0.271	0.29	0.166	0.19	0.29];
slny = [0.478   0.480   0.498   0.510   0.520   0.531   0.532   0.534   0.534   0.544   0.541   0.550   0.568   0.568   0.564   0.571   0.558   0.576   0.601   0.628   0.667   0.478   0.484   0.667];
roh  = [1.000   0.999   0.989   0.963   0.922   0.898   0.890   0.871   0.852   0.831   0.785   0.735   0.628   0.534   0.411   0.331   0.289   0.261   0.200   0.174   0.174   1.000   0.691   0.174];
tlny = [0.219   0.219   0.235   0.258   0.292   0.286   0.280   0.249   0.240   0.215   0.217   0.214   0.227   0.255   0.296   0.296   0.326   0.297   0.359   0.428   0.485   0.219   0.203   0.485];

c = 1.88;
n = 1.18;

% Magnitude dependence

if M <= 5.5;
    fmag = c0(ip) + c1(ip) * M;
elseif M<=6.5;
        fmag = c0(ip) + c1(ip) * M + c2(ip) * (M - 5.5); 
else    
        fmag = c0(ip) + c1(ip) * M + c2(ip) * (M - 5.5) + c3(ip) * (M - 6.5); 
end;

% Distance dependence

fdis = (c4(ip) + c5(ip) * M) * log(sqrt(Rrup^2 + c6(ip)^2));

% Style of faulting

fflt = c7(ip) * Frv * ffltz + c8(ip) * Fnm;

% Hanging-wall effects

fhngz = ((20 - Ztor)/20) * (Ztor >= 0 && Ztor < 20);

fhngdelta = (delta <= 70) + ((90 - delta)/20) * (delta > 70);

fhng = c9(ip) * fhngr * fhngm * fhngz * fhngdelta;

% Site conditions

if Vs30 < k1(ip);
    A1100 = CB_2008_nga (M, 0, Rrup,  Rjb,Ztor, delta, lambda, 1100, Zvs, arb);
    fsite = c10(ip) * log(Vs30/k1(ip)) + k2(ip) * (log(A1100 + c * (Vs30/k1(ip))^n) - log(A1100 + c));
elseif Vs30 < 1100;
    fsite = (c10(ip) + k2(ip) * n) * log(Vs30/k1(ip));
else    
    fsite = (c10(ip) + k2(ip) * n) * log(1100/k1(ip));
end;

% Sediment effects

if Zvs < 1;
    fsed = c11(ip) * (Zvs - 1);
elseif Zvs <= 3;
    fsed = 0;
else    
    fsed = c12(ip) * k3(ip) * exp(-0.75) * (1 - exp(-0.25 * (Zvs - 3)));
end;

% Median value
Sa = exp(fmag + fdis + fflt + fhng + fsite + fsed);

% Standard deviation computations
if (Vs30 < k1(ip));
    alpha = k2(ip) * A1100 * ((A1100 + c*(Vs30/k1(ip))^n)^(-1) - (A1100+c)^(-1));
else
    alpha = 0;
end;

slnaf=0.3;
slnpga=slny(length(slny)-2);
slnab=sqrt(slnpga^2    -slnaf^2);
slnyb=sqrt(slny(ip)^2  -slnaf^2);
sigmat = sqrt(slny(ip)^2 + alpha^2*slnab^2 + 2*alpha*roh(ip)*slnyb*slnab);
tau = tlny(ip);

sigma=sqrt(sigmat^2 + tau^2);
if arb == 1;
    sigma = sqrt(sigma^2 + sigmac(ip)^2);
end;