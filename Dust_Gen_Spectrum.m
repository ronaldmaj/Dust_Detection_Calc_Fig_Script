function [V2d, fr] = Dust_Gen_Spectrum(S,v_kms_1,C,tau_r,tau_d,alpha,...
    mmin,mmax,fmin,fmax,n,Q,f)
    % This function takes in variables and functions in order to predict &
    % output the voltage power spectrum due to the impact of dust particles
    % onto the surface of a CubeSat (or other spacecraft) or its antenna/s. 
    
    % This function is general and can be applied to a number of different 
    % scenarios, specifically, the function for the amount of charge
    % released when the dust particle impacts the spacecraft 'Q' and the 
    % distribution of flux per kilogram of the dust particles in the
    % environment of interest 'f' are left to the user to choose.
    
    % Reference: "Voyager 2 at Uranus: Grain impacts in the ring plane"
    % N. Meyer-Vernet, M.G. Aubier and B.M. Pedersen (1986)
    
    %% INPUTS:
    % S         - surface area of spacecraft                    (m^2)
    % v_kms_1   - speed of dust particles wrt to the CubeSat    (km.s^{-1})
    % C         - capacitance of antenna                        (F)
    % tau_r     - rise time constant of signal                  (s)
    % tau_d     - decay time constant of signal                 (s)
    % alpha     - fraction of charge contributing to signal         
    % mmin      - minimum mass of particles considered          (kg)
    % mmax      - maximum mass of particles considered          (kg)
    % fmin      - minimum of frequency range                    (Hz)
    % fmax      - maximum of frequency range                    (Hz)
    % n         - number of data points
    % N_I_min   - minimum impact rate cut-off                   (s^-1)
    % Q         - charge released on impact(function handle of mass, speed)
    % f         - flux of dust particles (funtion handle)  (m^2.s^-1.kg^-1)
    %             In case of single dust particle set f = 0
    
    %% OUTPUTS:
    % V2d       - dust spectrum vector                          (V2.Hz{-1})
    % fr        - frequency vector separated logarithmically    (Hz)
    
    %% Input Functions
    % Firstly, it is assumed that the charge function Q is a function of 
    % the mass of the particle and speed, having the general form:
    %                   Q(m,v_kms_1) = A * m^a * v_kms_1^b
    % where 'A' is a constant, 'a','b' are the coefficents for the 
    % variables 'm', the mass in kg, and 'v_kms_1' the speed in km.s^-1
    %
    % Secondly, the flux per kilogram function 'f' is also assumed to be a
    % simple power law and a function of mass 'm' only. This represents the
    % flux of the dust particles on a surface in the ionosphere per 
    % kilogram mass of the particle itself. Therefore we expect the form to
    % be:
    %                           f(m) = A * m^c 
    % where A is a constant and we assume c<0 as the coefficient of 'm' the 
    % mass in kg.
    % Therefore we expect, for example, that the flux of smaller
    % particles to be larger than the flux of larger particles. To find the
    % flux on its own, integrate over the mass range of interest to find
    % the number of particles per metre squared that impact the surface.
    % This assumption of a simple power law for the flux distribution is
    % important for calculating the critical mass
    
    %% Spectrum Calculation
    % To determine the voltage power spectrum from the dust impacts, we use
    % the following general expression:
    % 
    % V2d = 2 * N * |V(omega)|^2
    % 
    % where V(omega) is the Fourier transform of the voltage signal V(t)
    %
    % V(t) = H(t) * (1-exp(-t/tau_r)) * exp(-t/tau_d)
    %
    % where H(t) is the Heaviside step function and tau_r and tau_d are the
    % rise time and decay time constant respectively, and omega is defined:

    fr = logspace(log10(fmin),log10(fmax),n);
    omega = 2*pi*fr;
    
    % This function can either find the spectrum from a single dust
    % particle impact or the average from a distribution of dust particles
    % of varying sizes. Therefore setting mmin = mmax will initiate the 
    % single particle case otherwise the function will proceed to integrate
    % over the mass range chosen. Therefore for a single impact:
    if mmax == mmin
        V2d = (alpha^2*Q(mmax,v_kms_1)^2)./(2*C^2*tau_r^2.*omega.^4).*...
            ((tau_r^2*omega.^2)./(tau_r^2*omega.^2+1)).*...
            ((tau_d^2*omega.^2)./(tau_d^2*omega.^2+1)); 
    else
    % Else we are not dealing with impacts from identical particles but 
    % rather impacts from dust particles of a range of different sizes and
    % impact rates. Therefore we effectively need to find the integral:
    % 
    % V2d = 2 * integral(N(m)*|V(m,omega)|^2, mmin,mmax)
    %
    % = integrate((alpha^2*f(m)*S*Q(m)^2)./(2*C^2*tau_r^2*omega.^4).*((...
    % tau_r^2*omega.^2)./(tau_r^2*omega.^2 + 1)).*((tau_d^2*omega.^2)./(...
    % tau_d^2*omega.^2 + 1)) ,mmin,mmax)
    %
    % = (alpha^2*S)./(2*C^2*tau_r^2*omega.^4).*((tau_r^2*omega.^2)./(...
    % tau_r^2*omega.^2+1)).*((tau_d^2*omega.^2)./(tau_d^2*omega.^2+1)).*...
    % integrate(f(m)*Q(m)^2, mmin, mmax)
    
    % The factor not dependent on m:
        V2d_C = (alpha^2)./(2*C^2*tau_r^2*omega.^4).*...
            ((tau_r^2*omega.^2)./(tau_r^2*omega.^2+1)).*...
            ((tau_d^2*omega.^2)./(tau_d^2*omega.^2+1)); 
    
    % Integral function handle:
        Int = @(m) (f(m)*S.*Q(m,v_kms_1).^2);
    
        V2d_int = integral( @(m) Int(m) ,mmin ,mmax);
    
        V2d = V2d_C*V2d_int;
        
    end
    
end
    
    