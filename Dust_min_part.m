function [min_mass,min_len,N_I] = Dust_min_part(v_kms_1,C,tau_r,tau_d,alpha,f,l,a,...
    N,T,P,S,F,rho,Q,A)
% Calculates the minimum size particle hit necessary to produce a dust
% voltage power value 'A' times above the shot noise at a given frequency
% for a dipole antenna. The mass and length of the particle are quoted as
% final results. This function also calculates the impact rate for such a
% particular to impact an object (such as a spacecraft) given the
% cumulative flux distribution (as a function of length L).
%
%% Inputs:
% v_kms_1   - speed of dust particles relative to the CubeSat   (km.s^-1)
% C         - capacitance of antenna                            (F)
% tau_r     - rise time constant of signal                      (s)
% tau_d     - decay time constant of signal                     (s)
% alpha     - fraction of charge contributing to signal         
% f         - frequency of interest  (single value)             (Hz)
% l         - length of one arm of the dipole antenna           (m)
% a         - radius of antenna wire                            (m)
% N         - electron density of plasma                        (m^-3)
% T         - Temperature of plasma                             (K)
% P         - potential of antenna                              (V)
% S         - surface area of object impacted by dust           (m)
% F         - cumulative dust flux distribution function @(L)   (m^-2.s^-1)
% rho       - density of dust to convert length to mass         (kg.m^-3)
% Q         - charge released on impact function   @(v_kms_1,m) (C)
% A         - desired factor of dust spectrum above the shot noise
%
%% Output:
% min_mass  - the minimum mass of the dust particle             (kg)
% min_len   - the minimum length (diameter) of the particle     (m)
% N_I       - given F and S, the impact rate on the object      (s^-1)

%% Calculations
% Calculates the voltage power spectrum V2S for shot noise in a
% thermal plasma with temperature T and density N, on an antenna of
% length l, thickness a, with wire dipole configuration. Using V2S, we see
% what mass particle will produce a dust signal at the given frequency that
% is A times above the shot noise value V2S

% Frequency converted into omega/w
    omega = 2*pi*f;

% Physical constants
    kB = 1.3806488e-23;     % Boltzmann constant
    me = 9.109383e-31;      % Mass of electron    
    e = 1.60217657e-19;     % Electron charge
    eps0 = 8.85418782e-12;  % Vacuum Permittivity 

% Plasma parameters
    LD = sqrt((kB*T*eps0)/(N*e^2));
    vT = sqrt(2*kB*T/me);   % Thermal velocity of electrons

% Surface area, impact rate, conversion between z and omega
    S_a=2*pi*a*l;
    N_e = ((4*pi)^(-0.5))*N*vT*S_a;
    z = @(x) omega./(x.*vT);     
    
% Functions in integrand
    F1 = @(x) ((sinint(x))-(sinint(2*x)./2)-((2*((sin(x/2)).^4))./x))./x;
    W = @(x) faddeeva1(x);
    epsL = @(k) (1 + ((1 + (1i)*sqrt(pi)*z(k).*W(z(k)))./((k.^2)*(LD^2))));

% Integral calculation
    Int = @(k) (((F1(k*l)).*((besselj(0,k*a)).^2))./(epsL(k)));
    kmax = 1e5;
    kmin = 1e-5;
    Integral1 = quadgk(Int,kmin,kmax,'MaxIntervalCount',1000000); 

% Impedance and calculation of shot noise
    Z = (((4i)/((pi^2)*eps0*omega))*Integral1);
    V2S = (2*(e^2)*N_e)*((abs(Z))^2)*exp(e*P/(kB*T));
    
% Dust voltage power spectrum 
% Firstly, the cumulative dust flux distribution is given as a function of
% the length of the particle. We assume that this length is the diameter of
% the dust and that it follows a simple power law form with coefficients of
% B and c:
%                     F(L) = B * L^c.
% We also assume the dust particles are spherical in shape so that we
% convert between mass and length via the function:
    L = @(m) ((6*m)/(rho*pi)).^(1/3);
% with rho the density of the dust.

% The factor in dust spectrum not dependent on m:
    V2d_C = (alpha^2)/(2*C^2*tau_r^2*omega^4)*...
            ((tau_r^2*omega^2)/(tau_r^2*omega^2+1))*...
            ((tau_d^2*omega^2)/(tau_d^2*omega^2+1)); 
% Dust spectrum as a function of m:
    V2d = @(m) V2d_C*Q(m,v_kms_1).^2;
    
% Solve for 'x' such that V2d(x) = A*V2S to find the minimum mass:
    syms x positive
    eqn = V2d(x) == A*V2S;
    min_mass = double(solve(eqn,x,'Real',true));
% Convert to length:
    min_len = L(min_mass);
% Find impact rate
    N_I = S*F(min_len);
end
    
    
    