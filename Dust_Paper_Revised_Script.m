% Reference: "Dust detection via voltage power spectroscopy on a Cubesat in
% Earth's ionosphere" Maj R. and Cairns I.H., JGR: Space Physics (2018) 

% This script outlines the assumptions and working for the calculation of
% the dust impact spectra for a CubeSat in the Earth's ionosphere. This
% includes the predicted spectra from Space Debris, Interplanetary dust as
% well as the dust cloud created from an In-Orbit Satellite Collision.
% This also includes the commands necessary for producing the figures that 
% were presented in the revised submission to JGR: Space Physics (00/00/00)

% Firstly, we define some parameters and constants to be used throughout 
% the rest of the script:

% Physical constants
kB = 1.3806488e-23;     % Boltzmann constant
me = 9.109383e-31;      % Mass of electron
mi = 1.67262e-27;       % Mass of proton
e = 1.60217657e-19;     % Electron charge
eps0 = 8.85418782e-12;  % Vacuum Permittivity
% Antenna/spacecraft parameters
l = 0.3;                % Antenna Length
a = 2e-4;               % Antenna Radius
S = 6*0.01*(2/(3*pi));  % Exposed surface area to dust impact of 1U cubesat
S_a = 2*pi*a*l;         % Exposed surface area of the antennas to impact

%% Calculation of Spectra for Dust Impacts on CubeSat from Space Debris
tic
% Data for the dust distribution for space debris in the ionosphere are 
% taken from the NASA ORDEM 3.0 model and in this case assumed to be saved 
% as .xlsx files in the present working directory.

% Here we have used outputs for 300, 800 and 1500 km altitude circular
% polar orbits which we can import using:
Dust_dia =xlsread('ORDEM_Dust_Dist_300km_inc_90degrees.xlsx',1,'A1:A501');
Dust_flux300 = xlsread('ORDEM_Dust_Dist_300km_inc_90degrees.xlsx',...
    1,'B1:B501');
Dust_flux800 = xlsread('ORDEM_Dust_Dist_800km_inc_90degrees.xlsx',...
    1,'B1:B501');
Dust_flux1500 = xlsread('ORDEM_Dust_Dist_1500km_inc_90degrees.xlsx',...
    1,'B1:B501');
% The first variable 'Dust_dia' is the dust particle diameter while the
% others are the cumulative debris flux in (m^-2.yr^-1), which is the 
% debris flux for particles of length L or larger. 
% To convert to (m^-2.s^-1) we divide by 31536000:
DF_300 = Dust_flux300/31536000;
DF_800 = Dust_flux800/31536000;
DF_1500 = Dust_flux1500/31536000;
% Fit a power law to the ORDEM data. Only particles of size 2.5e-5 m or 
% less are considered:
options3 = fitoptions('power1');
[Fit300, ~, ~] = fit((Dust_dia(1:41)),(DF_300(1:41)),'power1',options3);
[Fit800, ~, ~] = fit((Dust_dia(1:41)),(DF_800(1:41)),'power1',options3);
[Fit1500, ~, ~]=fit((Dust_dia(1:41)),(DF_1500(1:41)),'power1',options3);

% Assign the coefficients from the power law fits to variables:
Coef300 = coeffvalues(Fit300);
Coef800 = coeffvalues(Fit800); 
Coef1500 = coeffvalues(Fit1500);
% Cumulative dust flux as simple power law wrt length based on ORDEM data
F300 = @(L) Coef300(1)*L.^(Coef300(2));
F800 = @(L) Coef800(1)*L.^(Coef800(2));
F1500 = @(L) Coef1500(1)*L.^(Coef1500(2));
% Determining the derivative (and then multiplying by -1) of the cumulative
% functions allows us to determine the non-cumulative functions, i.e.
% allowing us to choose what lengths to integrate over to determine the
% flux. Converting between length and mass means that the density of the
% particles would need to be known. We assume that the dust particles in
% the ionosphere are mostly due to space debris and therefore consist
% mainly of aluminium and titanium dioxide (paint), averaging these gives
% the density as rho = 3465 kg.m^{-3}. We also assume that they are
% spherical so that the mass of each particle is:
% m = V_p * rho (V_p the volume per particle, rho the density)
% m = (pi*rho*L^3)/6 (where we assume the radius of sphere is = L/2)
% So in this case:
rho = 3465;
m = @(L) (rho*pi/6)*(L.^3);
% Rearranging to make L the subject gives:
L = @(m) ((6*m)/(rho*pi))^(1/3);
% Therefore we can now convert the cumulative dust flux distributions to
% dust flux per kilogram distributions:
f300 = @(m) ((-Coef300(2))*Coef300(1))/3 * (6/(rho*pi))^(Coef300(2)/3)*...
    m.^(Coef300(2)/3 - 1);
f800 = @(m) ((-Coef800(2))*Coef800(1))/3 * (6/(rho*pi))^(Coef800(2)/3)*...
    m.^(Coef800(2)/3 - 1);
f1500 = @(m)((-Coef1500(2))*Coef1500(1))/3*(6/(rho*pi))^(Coef1500(2)/3)...
    *m.^(Coef1500(2)/3 - 1);
% Reference: "The New NASA Orbital Debris Engineering Model ORDEM 3.0"
% P.H. Krisko (2010)

% We now set a lower limit for the impact rate to calculate the 'average' 
% spectrum. This is done because much larger particles overwhelm the 
% contribution to the average spectrum even though they impact the 
% spacecraft very rarely. The minimum impact rate is set as:
N_I_min = 1.1e-3;
% which is equivalent to one impact every 15 minutes or at least a 
% 100*N_I_min percent chance of one particle impacting the spacecraft.

% To determine the maximum mass of dust particle for consideration in our 
% calculation of the dust spectrum at each height, we first calculate the 
% length of the particle that gives us that impact rate, which we can then 
% convert to a mass. The impact rate can be defined as:
%                       N_I = S * F(L) = S * A * L^-b 
% where S is the surface area, 'F(L)' is the distribution of the cumulative
% flux for particles of size 'L' or larger which has a simple power law
% relationship with 'A' constant, 'L' the diameter of the particle, and 'b'
% the coefficient where we assume b>0. 
% Calculating the coefficients based on the fits done for the ORDEM data:
b300 = -log(F300(exp(1))/F300(1));
b800 = -log(F800(exp(1))/F800(1));
b1500 = -log(F1500(exp(1))/F1500(1));
% then the length of the maximum particle size that will be considered
% in the calculation of the dust spectrum:
L300 = ((N_I_min)/(S*F300(1)))^(-1/b300);
L800 = ((N_I_min)/(S*F800(1)))^(-1/b800);
L1500 = ((N_I_min)/(S*F1500(1)))^(-1/b1500);
% finally converting the length to a mass:
m_cri300 = m(L300);
m_cri800 = m(L800);
m_cri1500= m(L1500);

% We estimate the relative speed of the spacecraft wrt the dust particle by
% considering the orbital speed of the spacecraft at each altitude:
v_orb = @(r) sqrt(6.67408e-11*5.972e24/(6371e3+r));
v_orb_300 = v_orb(300e3);
v_orb_800 = v_orb(800e3); 
v_orb_1500 = v_orb(1500e3);
% which we assume travels in a polar orbit, while the dust particles travel
% equatorially (at the same speed). Therefore the relative speed is:
v_kms_300 = v_orb_300*sqrt(2)/1000;
v_kms_800 = v_orb_800*sqrt(2)/1000;
v_kms_1500 = v_orb_1500*sqrt(2)/1000;
% We also need an expression for the charge released when the dust particle
% impacts the surface of the spacecraft. Here we use:
Q = @(m,v_kms_1) 0.704019*(v_kms_1^3.48)*(m.^1.02);
% which has been determined experimentally from impacts onto Al targets
% Reference: 'Meteoroid impacts on spacecraft: sporadics, streams, and the 
% 1999 Leonids', N. McBride and J.A.M. McDonnell (1999).
% The plasma cloud that is created by the dust impact expands out at speed:
v_ex = @(v,rho,rho_T) v/(1+sqrt(rho/rho_T));
% where rho_T is the density of the target material, which we consider here
% to be Aluminium:
rho_T_Al = 2700;
% and v is the relative speed between the particle and the spacecraft.
% Therefore at each altitude we have:
v_ex_300 = v_ex(v_orb_300*sqrt(2),rho,rho_T_Al);
v_ex_800 = v_ex(v_orb_800*sqrt(2),rho,rho_T_Al);
v_ex_1500 = v_ex(v_orb_1500*sqrt(2),rho,rho_T_Al);
% From which we can now calcualte the rise time at each altitude, using the
% average minimum distance from impact site to antenna, d = 0.0962 m
tau_r_300 = 0.0962/v_ex_300;
tau_r_800 = 0.0962/v_ex_800;
tau_r_1500 = 0.0962/v_ex_1500;

% Capacitance is different at each altitude due to the plasma sheath. To
% calculate it we use the values of electron density and temperature at
% each height in order to calculate the Debye length and therefore the
% capacitance using the low-frequency approximation:
LD = @(N,T) sqrt((kB*T*eps0)./(N*e^2)); % Debye length
C_A = @(Ld) (2*pi*eps0*l)./(log((Ld+a)/a)); % Capacitance Low-Freq. Approx.
% Loading electron temperatures and densities at each height considered:
load('Ne300D.mat');
load('Te300D.mat');
Xls_Ne_800 = xlsread('Ne_800km.xlsx');
Xls_Te_800 = xlsread('Te_800km.xlsx');
Xls_Ne_1500 = xlsread('Ne_1500km.xlsx');
Xls_Te_1500 = xlsread('Te_1500km.xlsx');
Ne800D = Xls_Ne_800(2:end,2:end);
Te800D = Xls_Te_800(2:end,2:end);
Ne1500D = Xls_Ne_1500(2:end,2:end);
Te1500D = Xls_Te_1500(2:end,2:end);
% Averages of electron density and temperature:
Ne300 = mean(mean(Ne300D));
Te300 = mean(mean(Te300D));
Ne800 = mean(mean(Ne800D));
Te800 = mean(mean(Te800D));
Ne1500 = mean(mean(Ne1500D));
Te1500 = mean(mean(Te1500D));
% Calculate Debye Lengths for each position around the globe
LD300 = LD(Ne300D,Te300D);
LD800 = LD(Ne800D,Te800D);
LD1500 = LD(Ne1500D,Te1500D);
% Calculate capacitance around globe
C_A_300 = C_A(LD300);
C_A_800 = C_A(LD800);
C_A_1500 = C_A(LD1500);
% Average capacitances:
C_300 = mean(mean(C_A_300));
C_800 = mean(mean(C_A_800));
C_1500 = mean(mean(C_A_1500));
% The gain in the capacitive coupling regime can be calculated using the
% capacitance of the sheath. Firstly, the capacitance of the load (or the
% capacitance of the receiver and the base of the enclosure) is approx.:
C_L = 10e-12;
% Formula for gain:
G = @(C_S) C_S./(C_S+C_L);
G_300 = G(C_300);
G_800 = G(C_800);
G_1500 = G(C_1500);
% And the resistance of the load is estimated as ~100Mohm:
R_L = 100e6;
% These estimates are both based on those from the CASSIOPE e-POP Radio 
% Receiver Instrument or RRI:
% Reference: "The e-POP Radio Receiver Instrument on CASSIOPE", James et.
% al. (2015)

% Decay time 
tau_d = R_L*C_L;
% alpha - Fraction of released charges that reach the antenna/s
alpha = 0.0055;
% Minimum particle mass considered
m_min = m(1e-8);
% Maximum frequency considered (Hz)
f_max = 1e8;
% Minimum frequency considered (Hz)
f_min = 1e4;
% Number of elements in frequency vector
n_f = 500;

% Calculate the dust spectrum for space debris at 300, 800 and 1500 km for
% a CubeSat in a circular polar orbit:
[V2d300, ~] = Dust_Gen_Spectrum(S,v_kms_300,C_300,tau_r_300,tau_d,alpha,...
    m_min,m_cri300,f_min,f_max,n_f,Q,f300);
[V2d800, ~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri800,f_min,f_max,n_f,Q,f800);
[V2d1500, fr] = Dust_Gen_Spectrum(S,v_kms_1500,C_1500,tau_r_1500,tau_d,alpha,...
    m_min,m_cri1500,f_min,f_max,n_f,Q,f1500);
disp('Dust Debris Spectra Calculation:')
toc

%% Calculation of Spectrum for Dust Impacts from Interplanetary Sources
tic
% For the interplanetary sources of dust, we use a different expression for
% the cumulative dust distribution, namely that from Grun et. al. (1985)
% and then calculate the flux per kilogram and critical particle mass based
% off of that.
% Reference: "Collisional Balance of the Meteoritic Complex" E.Grun et. al. 
% (1985)

% Firstly, the cumulative distribution of dust particles from Grun et. al.
% (1985) is [Equation (A3)]:
% Constants:
c_4 = 2.2e3;
c_5 = 15;
c_6 = 1.3e-9;
c_7 = 1e11;
c_8 = 1e27;
c_9 = 1.3e-16;
c_10 = 1e6;
% Function handle:
F_int = @(m) (c_4*m.^0.306 + c_5).^-4.38 + ...
    c_6*(m + c_7*m.^2 + c_8*m.^4).^-0.36 + c_9*(m + c_10*m.^2).^-0.85;


% Find where the impact rate goes down to N_I_min (if at all) in order to 
% calculate the critical mass, m_cri, for this distribution and spacecraft
% specifics. We choose a very large maximum mass for our search range such
% that the critical mass will be found, i.e. the mass of 5 micron particle 
m_vec = logspace(log10(m_min),log10(m(5e-6)),100001);
N_I = S*F_int(m_vec);
N_cri_idx = sum(N_I > N_I_min);
m_cri = m_vec(N_cri_idx);

% Now we find the dust flux per kilogram distribution by differentiating 
% and multiplying by (-1) as for the space debris case:
f_int = @(m) ((2.2e3*m.^(0.306) + 15).^-5.38)*2948.62 + (1.872e18*m.^3 +...
    93.6*m + 4.68e-10).*(m + 1e11*m.^2 + 1e27*m.^4).^-1.36 + ...
    (2.21e-10*m + 1.105e-16).*(m+1e6*m.^2).^-1.85;

% The speed of impact of the particles with the spacecraft will differ in
% the case of interplanetary dust compared to space debris. Using Grun
% et.al. (1985) we use the same estimate of speed:
v_int = 20;

% We also need to recalcualte the rise time due to the different density of
% interplanetary particles. In Grun et. al. (1985) it is argued that
% densities of rho_Int = 2500 kg.m^-3 are a good average estimate:
rho_int = 2500;
% Therefore we can calculate the rise time from this as:
v_ex_int = v_ex(v_int*1000,rho_int,rho_T_Al);
% From which we can now calcualte the rise time, using the average minimum
% distance from impact site to antenna, d = 0.0962 m
tau_r_int = 0.0962/v_ex_int;
% For later use we can also calculate the cumulative of dust particles from
% Grun et. al. (1985) as a function of L instead:
m_int = @(L) (rho_int*pi/6)*(L.^3);
F_int_L = @(L) F_int(m_int(L));

% And finally calculate the dust spectrum for interplanetary flux:
[V2dInt300, ~] = Dust_Gen_Spectrum(S,v_int,C_300,tau_r_int,tau_d,alpha,...
    m_min,m_cri,f_min,f_max,n_f,Q,f_int);
[V2dInt800, ~] = Dust_Gen_Spectrum(S,v_int,C_800,tau_r_int,tau_d,alpha,...
    m_min,m_cri,f_min,f_max,n_f,Q,f_int);
[V2dInt1500, ~] = Dust_Gen_Spectrum(S,v_int,C_1500,tau_r_int,tau_d,alpha,...
    m_min,m_cri,f_min,f_max,n_f,Q,f_int);

disp('Interplanetary Dust Spectrum Calculation:')
toc

%% Calculation of Spectra using Different Sampling (Averaging) Times
tic
% Although we have used N_I_min = 1.1e-3 as our minimum impact rate for the
% previous spectra, this need not be the case. We can also set rates that
% are larger or smaller, in effect changing the sampling (or averaging)
% time over which dust impacts are considered.

% In this section, we will calculate the dust spectra for a CubeSat at an 
% altitude of 800km using sampling times of:
% 1 minute
% 15 minutes
% 60 minutes
% 720 minutes (12 hours)
% This will reveal the effect that the sampling time has on the 'average'
% spectrum that is calculated. Therefore, firstly we define the impact
% rates to obtain the sampling times mentioned:
N_I_min_1min = 1/60;
N_I_min_60min = 1/3600;
N_I_min_720min = 1/43200;
% from which we can calculate the critical length of particle for each:
L800_1min = ((N_I_min_1min)/(S*F800(1)))^(-1/b800);
L800_60min = ((N_I_min_60min)/(S*F800(1)))^(-1/b800);
L800_720min = ((N_I_min_720min)/(S*F800(1)))^(-1/b800);
% and then convert the length to a critical mass:
m_cri800_1min = m(L800_1min);
m_cri800_60min = m(L800_60min);
m_cri800_720min= m(L800_720min);

% This now allows us to calculate the spectra 
[V2d800_1min,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri800_1min,f_min,f_max,n_f,Q,f800);
V2d800_15min = V2d800;
[V2d800_60min,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri800_60min,f_min,f_max,n_f,Q,f800);
[V2d800_720min,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri800_720min,f_min,f_max,n_f,Q,f800);

disp('Dust Spectra with different sampling times calculation:')
toc

%% Calculation of Single Particle Impact Spectrum
tic
% We can also calculate the spectrum from the impact of one particle onto
% the spacecraft. Here we consider a single 5 micron sized particle
% impacting the spacecraft and assume that alpha and all other variables
% remain the same (f = 0 as there is no flux distribution function
% necessary, mmax = mmin as there is just one particle and S is not used in
% the function either but it is set as it was before):
[V2d_one5um_300,~] = Dust_Gen_Spectrum(S,v_kms_300,C_300,tau_r_300,tau_d,alpha,...
    m(5e-6),m(5e-6),f_min,f_max,n_f,Q,0);
[V2d_one5um_800,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m(5e-6),m(5e-6),f_min,f_max,n_f,Q,0);

disp('Single particle spectrum calculation:')
toc

%% Calculation of Dust Spectra After Satellite Collision
tic
% After a collision two satellites, such as the 2009 Iridium-COSMOS
% collision, a large amount of new debris is created that can range from
% the very small to the very large. These particles will obviously lead to
% some dust impacts on other satellites that may go through the orbital
% path that the collided satellites followed. In the case here, we use a
% model from H.K. Springer et. al. (2010) to represent the dust particle
% distribution. This can be converted to a flux per kilogram distribution
% by converting to a mass, differentiating to obtain the number of
% particles per kilogram, multiplying by the relative speed of the
% particles wrt the spacecraft and then dividing by the volume that the
% dust particles take up. 
% Reference: "Satellite Collision Modeling with Physics-Based Hydrocodes: 
% Debris Generation Predictions of the Iridium-Cosmos Collision
% Event and Other Impact Events" H.K. Springer et. al. (2010)

% In Table 2 of H.K. Springer et. al. (2010) there is a simple "average 
% power-law fit" for the situation of maximum overlap of the satellies 
% during their collision. This is what we use for our calculation, i.e. 
% a_avg = 9.15, b_avg = -1.35 where N_cum = aL^b:
%                       N_cum(L) = 9.15*L^-1.35
% This is a cumulative distribution. Converting to a function of mass:
%                   N(m) = 9.15*(0.081991 * m^(1/3))^(-1.35)
%                        = 267.8163*m^(-0.45)
% Differentiating and multiplying by (-1) gives:
n_n = @(m) 120.5173*m.^(-1.45);
% which is the number of particles per unit mass (kg^-1).

% We can then use this to obtain the flux per kilogram distribution by
% noting:
%                        f(m) = (n(m) * v_ms_1) / V
% where v_ms_1 is the relative velocity of the dust to the CubeSat in
% m.s^(-1), and V is the volume that the dust particles from the satellite 
% collision take up in space in (m^3). 

% We are considering the CubeSat passing through the dust cloud from the
% collision at different times. We will assume that the collision takes 
% place at altitude 800 km. Therefore the volume of the dust cloud from
% after the collision will vary and necessitate separate flux distribution 
% functions. We can calculate the volume by considering the velocity at
% which the particles expand out from their point of collision.
% Using Table 3 from Springer et. al. (2010), the Cosmos satellite had mass 
% 898 kg and initial kinetic energy of 3.92e7 kJ. Under the maximum overlap 
% collision scenario, 58.82% of this mass and 11.82% of this energy is 
% present in particles of size less than 1cm after the collision. Therefore 
% if we assume that the kinetic energy is evenly distributed by mass we can
% use 
%                        KE = (1/2) * m* v_{exp}^2 
% to determine that the particles of the dust cloud expand out spherically 
% at a speed v_{exp} = 4.19 km.s^-1 relative to the Cosmos satellite after 
% the initial collision, with the center of mass following the orbital path
% of the Cosmos satellite at 7.46 km.s^-1 relative to the Earth:
vexp = sqrt((2*3.92e10*0.1182)/(898*0.5882));
% Therefore, the maximum relative speed between the CubeSat and the debris 
% cloud is v ~= 15 km.s^-1. Here we retain the average estimate for the
% relative speed between the spacecraft and dust particles at 10.9 km.s^-1:
v_ms_1 = 10900;

% To calculate the volume of the expanding dust cloud sphere we define:
V_sph = @(r) (4/3)*pi*r^3;
% with respect to the radius 'r'.
% Using v_{exp} = 4.19 km.s^-1 and considering a number of times after the
% collision, namely 10 sec, 30 sec, and 1 min after the collision, we
% expect that the maximum radial distances that the dust would have reached
% are 41.9 km, 126 km and 251 km, respectively. Calculating the volume of 
% these spheres:
V_sph41_9 = V_sph(41900);
V_sph126= V_sph(125700);
V_sph251= V_sph(251400);
% Now we can calculate the dust flux per kilogram distributions for these:
f_sat41_9 = @(m) (n_n(m)*v_ms_1)./V_sph41_9;
f_sat126 = @(m) (n_n(m)*v_ms_1)./V_sph126;
f_sat251 = @(m) (n_n(m)*v_ms_1)./V_sph251;
% Now calculating the critical mass based on
m_cri_sat41_9 = ((0.45*N_I_min*V_sph41_9)/(120.5173*S*v_ms_1))^(-1/0.45);
m_cri_sat126 = ((0.45*N_I_min*V_sph126)/(120.5173*S*v_ms_1))^(-1/0.45);
m_cri_sat251 = ((0.45*N_I_min*V_sph251)/(120.5173*S*v_ms_1))^(-1/0.45);

% Now we can calculate the spectra for each time period:
[V2d_sph_41_9km,~]= Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri_sat41_9,f_min,f_max,n_f,Q,f_sat41_9);
[V2d_sph_126km,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri_sat126,f_min,f_max,n_f,Q,f_sat126);
[V2d_sph_251km,~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_800,tau_d,alpha,...
    m_min,m_cri_sat251,f_min,f_max,n_f,Q,f_sat251);

disp('Dust Spectra after Satellite Collision calculation:')
toc

%% Justification of Constant Spherical Expansion of Dust After Collision
tic
% The assumption of uninhibited spherical expansion is a simplification as
% gravity and the drag force will mean that the particles will not expand
% out at a constant rate and with the same speed in all directions. However
% due to the short time scales used here, the assumption is justified. For
% particles moving in the direction of the Earth, the force of gravity will
% accelerate the particles moving them further away from the orbital path
% than would be the case for a constant velocity. For particles moving away
% from the Earth, gravity will decelerate and therefore the particles will
% not reach as far as they would have with a constant velocity.
% To calculate these distances we can use a vector 't' between 0.1 - 60 sec
% in timesteps of 0.1 sec:
t = linspace(0.1,60.1,601);
t_s = t(2)-t(1);
% and record the distance traversed at each time step:
s_record_down = zeros(size(t));
s_record_up = zeros(size(t));
% but first we need to establish the initial velocity:
vi_down = vexp; % m/s
vi_up = vexp;
% and acceleration at the given height:
H = 800e3;
accel = (9.806*6371000^2)./((6371000+H).^2);
% and accelerations that can be changed:
accel_down = (9.806*6371000^2)./((6371000+H).^2);
accel_up = (9.806*6371000^2)./((6371000+H).^2);
% which changes with time, therefore we run a for loop:
for i = 1:length(t)-1
    % Final velocity after timestep t_s
    vf_up = vi_up - accel_up*t_s;
    vf_down = vi_down + accel_down*t_s;
    % Distance traversed 
    s_i_up = ((vf_up+vi_up)/2)*t_s;
    s_i_down = ((vf_down+vi_down)/2)*t_s;
    % Record:
    s_record_up(i) = s_i_up;
    s_record_down(i) = s_i_down;
    % Update acceleration based on new height above the Earth:
    accel_up = (9.8*6371000^2)./((6371000+H+sum(s_record_up)).^2);
    accel_down = (9.8*6371000^2)./((6371000+H-sum(s_record_down)).^2);
    % and also the initial velocity for the next iteration:
    vi_up = vf_up;
    vi_down = vf_down;
end
% Therefore the distance traversed in 60 sec upwards and downwards is:
s_up = sum(s_record_up);
s_down = sum(s_record_down);
% which is very close to the case of a constant speed (about 6% diff.):
s_const = vexp*60;
diff_up = s_const - s_up;
diff_down = s_down - s_const;
frac_up = 1 - s_up/s_const;
frac_down = 1 - s_const/s_down;
% For the case of particles moving in opposition to the orbital motion,
% their speed relative to the Earth will only be 3.27 kms^-1 in the orbital
% direction (perpendicular to the surface of the Earth):
vopp = v_orb_800 - vexp;
% At orbital speed, i.e. 7.46 kms^-1, a particle will fall at the same rate
% as the curvature of the orbital path. At 800 km altitude, the 
% acceleration due to gravity is 7.74 ms^-2, therefore after 1 sec a 
% particle will fall 3.87 m:
s_k = @(vi,a,t) vi*t + 0.5*a*t^2;
fall_1sec = s_k(0,accel,1);
% and if traveling at orbital speed, will travel in an arc of 7.46e3 m 
% along the orbital path. The angle between the tangent and the arc is:
angle_arc = asin(fall_1sec/v_orb_800);
% based on the small angle approximation.
% For a particle traveling only at 3.27 kms^-1, the arc will be shorter, at
% 3.27e3 m, and therefore not remain at a 800 km altitude, but fall closer 
% to the Earth. For a particle in orbit, an arc of 3.27e3 m corresponds to
% a fall from the tangent of only 1.69 m:
fall_arc_3_27e3 = vopp*angle_arc;
% but a particle traveling at 3.27 kms^-1 will still fall 3.87 m relative 
% to the tangent and therefore will come closer to the Earth by 2.17 m:
fall_opp_1sec = fall_1sec - fall_arc_3_27e3;
% After 60 sec, this will amount to a fall of only 130 m:
fall_opp_60sec = fall_opp_1sec*60;
% which is only 0.02\% of the total height, 800 km.

disp('Approximation of Spherical Expansion calculation:')
toc

%% Calculation of Quasi-Thermal Noise and Shot Noise 
tic
% Lastly, we will commpare the spectrum expected from dust impacts against
% two effects that are known to create detectable voltage power spectra:
% quasi-thermal noise and shot noise.

% The negative potential for the spacecraft means a saturated photoemission
% current is expected. In Laakso et. al. (1995) the photoelectron current
% density is given as j_ph^0 = 4 nA cm^-2 (= 2.50e14 charges m^-2 s^-1)
% however other studies give a range of 6-8 nA cm^-2. Therefore here we
% have chosen j_ph^0 = 6.5 nA cm^-2 (= 4e14 charges m^-2 s^-1) as a middle
% ground between the two.
% Reference: "Plasma gradient effects on double-probe measurements in the 
% magnetosphere" H. Laakso, T. L. Aggson, and R. F. Pfaff Jr. (1995)
j = 4e14;

% ---------------------Load potentials from prior run ---------------------
load('Dust_Paper_Global_Potentials.mat')
% ---------------------------Calculate from scratch -----------------------

% We calculate the potential of both the spacecraft body and the antenna:

% Calculation of potential of the spacecraft at each height around globe
%[P300] = PotentialNew(Ne300D, Ne300D, Te300D, Te300D, 0.06, j, 300, 0);
%[P800] = PotentialNew(Ne800D, Ne800D, Te800D, Te800D, 0.06, j, 800, 0);
%[P1500]= PotentialNew(Ne1500D,Ne1500D,Te1500D,Te1500D,0.06, j, 1500,0);

% Calculation of potential of the antenna at each height around the globe
%[P_a300] = PotentialNew(Ne300D, Ne300D, Te300D, Te300D, S_a, j, 300, 1);
%[P_a800] = PotentialNew(Ne800D, Ne800D, Te800D, Te800D, S_a, j, 800, 1);
%[P_a1500]= PotentialNew(Ne1500D,Ne1500D,Te1500D,Te1500D,S_a, j, 1500,1);

%--------------------------------------------------------------------------

% The average spacecraft potential for a given height
AP300=mean(mean(P300));
AP800=mean(mean(P800));
AP1500=mean(mean(P1500));

% The average antenna potential for a given height
AP_a300=mean(mean(P_a300));
AP_a800=mean(mean(P_a800));
AP_a1500=mean(mean(P_a1500));

% Having the plasma conditions for the ionosphere at the three different
% heights we can now calculate the transition frequency.
% We can calculate the transition frequency from a resistively coupled to a
% capacitively coupled regime by considering the resistance of the sheath
% as well:
% Define correction factor for cylindrical wire surface:
A = 4/(pi*sqrt(2));
% Ambient ion current:
I_i = @(S,n_i,T_i,V) e*A*S*n_i.*sqrt((kB*T_i)/(2*pi*mi)).*...
(2*sqrt(-e*V./(kB*T_i))+exp(-e*V./(kB*T_i)).*erfc(-e*V./(kB*T_i)));
I_i300 = I_i(S_a,Ne300D,Te300D,P_a300);
I_i800 = I_i(S_a,Ne800D,Te800D,P_a800);
I_i1500 = I_i(S_a,Ne1500D,Te1500D,P_a1500);
% Photoemission current:
I_ph = @(S,j) (A*S/(2*sqrt(2)))*e*j;
I_pho = I_ph(S_a,j);
% Electron temperature in volts:
U_e = @(T_e) kB*T_e/e;
U_e300 = U_e(Te300D);
U_e800 = U_e(Te800D);
U_e1500 = U_e(Te1500D);
% Sheath resistance at each height and point around the globe:
R_S = @(U_e,I_ph,I_i) U_e./(I_ph+I_i);
R_S300 = R_S(U_e300,I_pho,I_i300);
R_S800 = R_S(U_e800,I_pho,I_i800);
R_S1500 = R_S(U_e1500,I_pho,I_i1500);
% Transition angular frequency at each point around globe:
tran_omega = @(R_S,C_S) 1./(R_S.*C_S);
tran_omega300 = tran_omega(R_S300,C_A_300);
tran_omega800 = tran_omega(R_S800,C_A_800);
tran_omega1500 = tran_omega(R_S1500,C_A_1500);
% Average transition frequency:
f_t_300 = mean(mean(tran_omega300))/(2*pi);
f_t_800 = mean(mean(tran_omega800))/(2*pi);
f_t_1500 = mean(mean(tran_omega1500))/(2*pi);

% ------------Load QTN and Shot Noise spectra from prior run---------------
load('Dust_Paper_QTN_Shot_Spectra.mat')
% ---------------Calculate QTN and Shot Noise from scratch-----------------

% Calculation of the voltage power spectra for QTN and shot noise at each
% altitude:
%[V2qtn300, V2s300, ~, ~] = plasma_QTN_shot_spectra(l,Ne300,Te300,a,...
%    0,f_min,f_max,n_f,1,AP_a300,1e5,1e-5,0.75);
%[V2qtn800, V2s800, ~, ~] = plasma_QTN_shot_spectra(l,Ne800,Te800,a,...
%    0,f_min,f_max,n_f,1,AP_a800,1e5,1e-5,0.75);
%[V2qtn1500,V2s1500,~,~] = plasma_QTN_shot_spectra(l,Ne1500,Te1500,a,...
%    0,f_min,f_max,n_f,1,AP_a1500,1e5,1e-5,0.75);

%--------------------------------------------------------------------------

% Although the above QTN spectra have been calculated over the 
% entire spectrum f = 1e4 Hz to 1e8 Hz, the assumptions in using the shot 
% noise expression are only true for f << f_p (much less than the plasma
% frequency). Therefore the code only caculates up to a fraction of the
% plasma frequency for the shot noise case. Here the cut-off is 0.75*f_p.

disp('QTN and Shot Noise spectra calculation:')
toc

%% Calculate/Display Numerical Results - Voltage Power and Voltage
tic
% Here we will calculate a number of variables.

% Power spectral density and voltage signal strength at 10 kHz (1e4):

% Function for calculating voltage signal at 10 kHz based on the
% definitions of frequency used in this script:
V_trap = @(a,b) sqrt(((a+b)/2)*(fr(39)-fr(1)));

% Space Debris
Vd300 = V_trap(V2d300(1),V2d300(39));
Vd800 = V_trap(V2d800(1),V2d800(39));
Vd1500 = V_trap(V2d1500(1),V2d1500(39));

% Interplanetary
VdInt300 = V_trap(V2dInt300(1),V2dInt300(39));
VdInt800 = V_trap(V2dInt800(1),V2dInt800(39));
VdInt1500 = V_trap(V2dInt1500(1),V2dInt1500(39));

% Satellite collision
Vd_41_9km = V_trap(V2d_sph_41_9km(1),V2d_sph_41_9km(39));
Vd_126km = V_trap(V2d_sph_126km(1),V2d_sph_126km(39));
Vd_251km = V_trap(V2d_sph_251km(1),V2d_sph_251km(39));

% Single particle impact
Vd_one5um_300 = V_trap(V2d_one5um_300(1),V2d_one5um_300(39));
Vd_one5um_800 = V_trap(V2d_one5um_800(1),V2d_one5um_800(39));

disp('Voltage Signal Calculation:')
toc

intro_V2_V = ['\nThe power spectral density and voltage signal strength'...
    ' at 10 kHz are \nshown in the table below for a number of dust '...
    'environments. The NASA\nORDEM 3.0 model is used to determine the '...
    'orbital debris environment\nat different altitudes. The '...
    'interplanetary dust is modeled\nfrom Grun et. al. (1985). The '...
    'collision scenario is from \nSpringer et. al. (2010), assumes a '...
    'COSMOS-Iridium type satellite\ncollision and that the dust cloud '...
    'is from a maximum overlap collision,\ni.e. the greatest '...
    'release of dust. \n \nVoltage Power is in units of (V^2.Hz^{-1}) '...
    '\nVoltage signal is in (V).\n\n'];

Dust_Environment = {...
    'Space Debris at 300 km';...
    'Space Debris at 800 km';...
    'Space Debris at 1500 km';...
    'Interplanetary Dust at 300 km';...
    'Interplanetary Dust at 800 km';...    
    'Interplanetary Dust at 1500 km';...
    '10 sec after Collision';...
    '30 sec after Collision';...
    '60 sec after Collision';...
    'One 5 micron particle (300km)';...
    'One 5 micron particle (800km)'};

Voltage_Power = [V2d300(1); V2d800(1); V2d1500(1); V2dInt300(1);...
    V2dInt800(1); V2dInt1500(1); V2d_sph_41_9km(1); V2d_sph_126km(1);...
    V2d_sph_251km(1); V2d_one5um_300(1);V2d_one5um_800(1)];

Voltage =[Vd300; Vd800; Vd1500; VdInt300; VdInt800; VdInt1500;Vd_41_9km;...
    Vd_126km; Vd_251km; Vd_one5um_300; Vd_one5um_800];

% Create table based on voltage power and voltage information above
T_V2_V = table(Dust_Environment,Voltage_Power,Voltage);

% Print information to console
fprintf(intro_V2_V)
disp(T_V2_V)

%% Calculate/Display Numerical Results - Minimum particle size and rate
% Size of particle necessary to produce spectrum 3x larger than shot noise
% and impact rate at which this can occur at given height is given below.
tic
% Set what factor the dust noise should be larger than the shot noise
B = 3;
% Set the frequency of interest
f = 1e4;
% Minimum particle and its rate at 300 km:
[min_mass300,min_len300,N_I_min300] = Dust_min_part(v_kms_300,C_300,tau_r_300,...
    tau_d,alpha,f,l,a,Ne300,Te300,AP_a300,S,F300,rho,Q,B);
% Same at 800 km
[min_mass800,min_len800,N_I_min800] = Dust_min_part(v_kms_800,C_800,tau_r_800,...
    tau_d,alpha,f,l,a,Ne800,Te800,AP_a800,S,F800,rho,Q,B);
% Same at 1500 km
[min_mass1500,min_len1500,N_I_min1500] = Dust_min_part(v_kms_1500,C_1500,...
    tau_r_1500,tau_d,alpha,f,l,a,Ne1500,Te1500,AP_a1500,S,F1500,rho,Q,B);

disp('Minimum Particle Size Calculation:')
toc

intro_min_dust = ['\nFor a single particle to be detected via voltage '...
    'power spectroscopy it \nmust rise above the shot noise for the '...
    'dipole antenna being considered. \nHere we show the minimum particle'...
    ' size that can create a dust voltage \npower spectrum 3 times larger'...
    ' than the shot noise for the given plasma \nconditions at 300 km, '...
    '800 km, and 1500km. \n \nHeight is above the surface of the Earth'...
    ' and is in (km) \nMass is mass of the particle and is in (kg) \n'...
    'Length is the diameter of the particle in (m) \nImpactRate is '...
    'for a CubeSat of 1U, i.e. a 10x10x10 cm cube, in (s^{-1}) \nTime '...
    'is the time taken for one dust impact to occur in (min)\n\n'];

Height = [300;800;1500];
Mass = [min_mass300; min_mass800; min_mass1500];
Length = [min_len300;min_len800; min_len1500];
ImpactRate = [N_I_min300;N_I_min800;N_I_min1500];

% Create a table based on the minimum particle impact info above
T_min_dust = table(Height,Mass,Length,ImpactRate);
T_min_dust.Time = (ImpactRate.^-1)/60;

% Print information to console
fprintf(intro_min_dust)
disp(T_min_dust)


%% Comparison of different charge release functions at 800 km
% As the surface of the hypothetical CubeSat may not be Al, we need to test
% whether the predictions made for Al hold up for other materials, and if
% there are changes, what are the new minimum sizes of particles necessary
% to create a detectable spectrum. We will use the parameters at 800 km
% altitude and compare the predicted dust spectra due to space debris,
% followed by comparing the minimum dust size and impact time at 800 km.

% Charge released upon impact functions from
% Reference: "Micrometeoroid impact charge yield for common spacecraft 
% materials" Collette et. al. (2014)
Q_MLI = @(m,v_kms_1) 1.7e-3*(v_kms_1^4.7)*m;
Q_Solar = @(m,v_kms_1) 4.7e-3*(v_kms_1^4.2)*m;
Q_BeCu = @(m,v_kms_1) 1.2e-2*(v_kms_1^3.8)*m;

% Due to the new target material, new estimates for the rise time must also
% be made as the speed at which the plasma expands out will change.
% MLI is multilayer thermal insulation and in Collette et. al. (2014) they
% use the MLI on the STEREO spacecraft in their experiments - indium tin
% oxide coated Teflon. Teflon's density is 2200 kg.m^-3 while that of ITO
% is 7120 kg.m^-3. Here we will assume 2% is ITO:
rho_MLI = 2220*0.98 + 7120*0.02;
% The STEREO spacecraft uses solar cells (or solar panels) made of GaAs 
% with a top surface made of cerium-doped glass 150 microns thick.
% Depending on the composition, cerium-doped glasses can have a range of
% densities from 2550 - 3729 kg.m^-3 (Reference: "Spectroscopic properties 
% of cerium-doped aluminosilicate glasses" Herrmann et. al. (2015)). Here
% we will use:
rho_Solar = 3000;
% Lastly, BeCu alloy is a copper alloy containing 0.5-3% beryllium. Copper 
% has a density of 8960 kg.m^-3 while beryllium's is 1850 kg.m^-3,
% therefore a 2% alloy would be:
rho_BeCu = 8960*0.98 + 1850*0.02;
% Now we can calculate the v_ex and therefore tau_r for each:
v_ex_MLI = v_ex(v_orb_800*sqrt(2),rho,rho_MLI);
v_ex_Solar = v_ex(v_orb_800*sqrt(2),rho,rho_Solar);
v_ex_BeCu = v_ex(v_orb_800*sqrt(2),rho,rho_BeCu);
% From which we can now calcualte the rise time at each altitude, using the
% average minimum distance from impact site to antenna, d = 0.0962 m
tau_r_MLI = 0.0962/v_ex_MLI;
tau_r_Solar = 0.0962/v_ex_Solar;
tau_r_BeCu = 0.0962/v_ex_BeCu;

% Voltage power spectra using different Q:
[V2d800_MLI, ~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_MLI,tau_d,alpha,...
    m_min,m_cri800,f_min,f_max,n_f,Q_MLI,f800);
[V2d800_Solar, ~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_Solar,tau_d,alpha,...
    m_min,m_cri800,f_min,f_max,n_f,Q_Solar,f800);
[V2d800_BeCu, ~] = Dust_Gen_Spectrum(S,v_kms_800,C_800,tau_r_BeCu,tau_d,alpha,...
    m_min,m_cri800,f_min,f_max,n_f,Q_BeCu,f800);

% Calculations for minimum mass/length and impact rate at 800 km with the
% new Q functions:
[min_mass800_MLI,min_len800_MLI,N_I_min800_MLI] = Dust_min_part(v_kms_800,...
    C_800,tau_r_MLI,tau_d,alpha,f,l,a,Ne800,Te800,AP_a800,S,F800,rho,Q_MLI,A);
[min_mass800_Solar,min_len800_Solar,N_I_min800_Solar] = Dust_min_part(...
    v_kms_800,C_800,tau_r_Solar,tau_d,alpha,f,l,a,Ne800,Te800,AP_a800,S,F800,rho...
    ,Q_Solar,A);
[min_mass800_BeCu,min_len800_BeCu,N_I_min800_BeCu] = Dust_min_part(...
    v_kms_800,C_800,tau_r_BeCu,tau_d,alpha,f,l,a,Ne800,Te800,AP_a800,S,F800,rho...
    ,Q_BeCu,A);

intro_min_dustQ = ['\nBelow we show the minimum particle'...
    ' size that can create a dust voltage \npower spectrum 3 times larger'...
    ' than the shot noise for different charge \nrelease functions at '...
    '800 km. These functions are derived for different \nmaterials '...
    'including Aluminium, MLI which is mulitlayer insulation, \nsolar '...
    'cell material, and BeCu, which is the alloy common for antennas.'...
    '\n \nQ_Function is the charge release material considered'...
    '\nMass_Q is mass of the particle and is in (kg) \n'...
    'Length_Q is the diameter of the particle in (m) \nImpactRate_Q is '...
    'for a CubeSat of 1U, i.e. a 10x10x10 cm cube, in (s^{-1}) \nTime_Q '...
    'is the time taken for one dust impact to occur in (min)\n\n'];

Q_Function = {'Al';'MLI (Insulation)';'Solar_Cell';'BeCu (Antenna)'};
Mass_Q = [min_mass800; min_mass800_MLI;min_mass800_Solar;min_mass800_BeCu];
Length_Q = [min_len800; min_len800_MLI; min_len800_Solar;min_len800_BeCu];
ImpactRate_Q =[N_I_min800;N_I_min800_MLI;N_I_min800_Solar;N_I_min800_BeCu];

% Create a table based on the minimum particle impact info above
T_min_dustQ = table(Q_Function,Mass_Q,Length_Q,ImpactRate_Q);
T_min_dustQ.Time_Q = (ImpactRate_Q.^-1)/60;

% Print information to console
fprintf(intro_min_dustQ)
disp(T_min_dustQ)


%% Plotting of Figures
% In this section we now plot the various figures used in the paper.
% The screen resolution used to create the plots was 1920x1080

% Figure 3 - Power Law Fits
B1 = Dust_dia(1:41);
B300 = DF_300(1:41);
B800 = DF_800(1:41);
B1500 = DF_1500(1:41);

Fig3 = figure('pos',[5 5  720  1000]);
Fig3.Name = 'Fig3-ORDEM_Fit_Plots';

subplot(3,2,1)
P3_1 = plot(Fit300,B1,B300,'predobs');
P3_1(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','log');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Particle Flux (m^{-2}s^{-1})');
title('300 km - ORDEM Output and Fit');
L3_1 = legend;
L3_1.Location = 'best';

subplot(3,2,2)
P3_2 = plot(Fit300,B1,B300,'stresiduals');
P3_2(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','linear');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Standard Deviations')
legend Location Best
title('300km - Std. Residuals')
L3_2 = legend;
L3_2.Location = 'best';

subplot(3,2,3)
P3_3 =plot(Fit800,B1,B800,'predobs');
P3_3(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','log');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Particle Flux (m^{-2}s^{-1})');
title('800 km - ORDEM Output and Fit');
L3_3 = legend;
L3_3.Location = 'best';

subplot(3,2,4)
P3_4 = plot(Fit800,B1,B800,'stresiduals');
P3_4(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','linear');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Standard Deviations')
title('800km - Std. Residuals')
L3_4 = legend;
L3_4.Location = 'best';

subplot(3,2,5)
P3_5 = plot(Fit1500,B1,B1500,'predobs');
P3_5(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','log');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Particle Flux (m^{-2}s^{-1})');
title('1500 km - ORDEM Output and Fit');
L3_5 = legend;
L3_5.Location = 'best';

subplot(3,2,6)
P3_6 = plot(Fit1500,B1,B1500,'stresiduals');
P3_6(1).DisplayName = 'ORDEM Output';
set(gca,'XScale','log','YScale','linear');
xlim([1e-5 2.5e-5]);
xlabel('Particle Diameter (m)')
ylabel('Standard Deviations')
title('1500km - Std. Residuals')
L3_6 = legend;
L3_6.Location = 'best';
%saveas(Fig3,Fig3.Name,'epsc')

% Figure 4 - Dust Spectra with Different Sampling Times
Fig4 = figure('pos',[10 10 675 950]);
Fig4.Name = 'Fig4-Dust_Spec_Samp_Times';
P4 = loglog(...
    fr,V2d800_1min,...
    fr,V2d800_15min,...
    fr,V2d800_60min,...
    fr,V2d800_720min,...
    'LineWidth',2);
ax4 = gca;
ax4.FontSize = 15;
ax4.FontWeight = 'bold';
ax4.YLim = [1e-30 1e-12];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold'); 
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
L4=legend(...
    'Dust Debris Spectrum - 1 min Sampling',...
    'Dust Debris Spectrum - 15 min Sampling',...
    'Dust Debris Spectrum - 60 min Sampling',...
    'Dust Debris Spectrum - 12 hr Sampling',...
    'Location','best');
%saveas(Fig4,Fig4.Name,'epsc')

% Figure 5 - Dust Spectra from Space Debris and Interplanetary Sources
Fig5 = figure('pos',[10 10 675 950]);
Fig5.Name = 'Fig5-Dust_Spec_Debris_IntPlan';
P5 = loglog(...
    fr,V2d300,...
    fr,V2d800,...
    fr,V2d1500,...
    fr,V2dInt300,...
    fr,V2dInt800,...
    fr,V2dInt1500,...
    'LineWidth',2);
for i = 4:6
    P5(i).LineStyle = '--';
end
ax5 = gca;
ax5.FontSize = 15;
ax5.FontWeight = 'bold';
ax5.YLim = [1e-30 1e-17];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold'); 
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
L5=legend(...
    'Dust Debris Spectrum (300km)',...
    'Dust Debris Spectrum (800km)',...
    'Dust Debris Spectrum (1500km)',...
    'Interplanetary Dust Spectrum (300km)',...
    'Interplanetary Dust Spectrum (800km)',...
    'Interplanetary Dust Spectrum (1500km)',...
    'Location','northeast');
%saveas(Fig5,Fig5.Name,'epsc')

% Figure 6 - Dust Spectra of CubeSat passing through Satellite Collision
Fig6 = figure('pos',[10 10 675 950]); 
Fig6.Name = 'Fig6-Dust_Spec_Sat_Coll';
P6 = loglog(...
    fr,V2d_sph_41_9km,...
    fr,V2d_sph_126km,...
    fr,V2d_sph_251km,...
    fr,V2dInt800,...
    'LineWidth',2);
P6(4).LineStyle = '--';
ax6 = gca;
ax6.FontSize = 15;
ax6.FontWeight = 'bold';
ax6.YLim = [1e-34 1e-20];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold'); 
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
L6=legend(...
    'Dust Spectrum 10 sec after Collision',...
    'Dust Spectrum 30 sec after Collision',...
    'Dust Spectrum 1 min after Collision',...
    'Interplanetary Dust Spectrum (800km)',...
    'Location','best');
%saveas(Fig6,Fig6.Name,'epsc')

% Figure 7 - QTN, Shot Noise and Dust Spectra (incl. Single Particle Hit)
Fig7 = figure('pos',[10 10 675 950]); 
Fig7.Name = 'Fig7-QTN_Shot_Dust_Spec';
P7 = loglog(...
    fr,(V2qtn300+V2s300),...
    fr,(V2qtn800+V2s800),...
    fr,(V2qtn1500+V2s1500),...
    fr,V2d300,...
    fr,V2d_sph_41_9km,...
    fr,V2d_one5um_300,...    
    fr,V2d_one5um_800,...
    'LineWidth',2);
P7(4).LineStyle = ':';
P7(5).LineStyle = ':';
P7(6).LineStyle = '-.';
P7(7).LineStyle = '-.';
ax7 = gca;
ax7.FontSize = 15;
ax7.FontWeight = 'bold';
ax7.YLim = [1e-23 1e-9];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold'); 
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
L7=legend(...
    'QTN + Shot Noise (300km)',...
    'QTN + Shot Noise (800km)',...
    'QTN + Shot Noise (1500km)',...
    'Dust Debris Spectrum - 300 km',...
    'Collision Dust Noise - 10 sec',...
    'Dust Noise - one 5\mum particle (300km)',...
    'Dust Noise - one 5\mum particle (800km)',...
    'Location','best');
box on
%saveas(Fig7,Fig7.Name,'epsc')

% Figure 8 - Comparison of different charge release functions - Q
Fig8 = figure('pos',[10 10 675 950]);
Fig8.Name = 'Fig8-Dust_Spec_Diff_Q';
SP8_1 = subplot(10,1,1:5);
P8_1 = loglog(...
fr,V2d800,...
fr,V2d800_MLI,...
fr,V2d800_Solar,...
fr,V2d800_BeCu,...
'LineWidth',2);
ax8_1 = gca;
ax8_1.FontSize = 15;
ax8_1.FontWeight = 'bold';
ax8_1.YLim = [1e-26 1e-18];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold');
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
'FontWeight','bold');
L8_1=legend(...
'Dust Debris Spectrum - Aluminum',...
'Dust Debris Spectrum - MLI',...
'Dust Debris Spectrum - Solar Cell',...
'Dust Debris Spectrum - BeCu Alloy',...
'Location','best');
box on

SP8_2 = subplot(10,1,7:10);
P8_2 = loglog(...
fr,V2d800,...
fr,V2d800_MLI,...
fr,V2d800_Solar,...
fr,V2d800_BeCu,...
'LineWidth',2);
ax8_2 = gca;
ax8_2.FontSize = 15;
ax8_2.FontWeight = 'bold';
ax8_2.XLim = [1e4 1.2e4];
ax8_2.YLim = [5e-22 1e-21];
xlabel('Frequency (Hertz)','FontSize',16,'FontWeight','bold');
ylabel('Power Spectral Density (V^2 Hz^{-1})','FontSize',16,...
'FontWeight','bold');
box on
%saveas(Fig8,Fig8.Name,'epsc')
