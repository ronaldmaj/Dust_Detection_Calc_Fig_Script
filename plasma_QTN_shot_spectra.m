function [V2QTN, V2S, Z_a, w_vec] = plasma_QTN_shot_spectra(l,N,T,a,...
    sph,fmin,fmax,n,dipole,P,k_max,k_min,C)
    % Code for calculating the voltage power spectrum similar to that done 
    % in the Meyer-Vernet/Perche 1989 Paper "Tool Kit for Antennae and 
    % Thermal Noise Near the Plasma Frequency". Produces a plot of V^2 at 
    % the end, and outputs this voltage spectrum and the frequency range.
    %
    % Inputs: 
    % l     - length of antenna                 (m)
    % N     - density of electrons              (m^-3)
    % T     - temperature of electrons          (K)
    % a     - radius of antenna                 (m)
    % sph   - spherical antenna?                (1 = yes, 0 = cylindrical)
    % fmin  - lower bound for frequency vector  (power of 10)
    % fmax  - upper bound for frequency vector  (power of 10)
    % n     - no. of elements in freq. vector
    % dipole- considering dipole antenna?       (1 = yes, 0 = mono)
    % P     - potential of the antenna          (V)
    % k_max - maximum k value in integral       (m^-1)
    % k_min - minimum k value in integral       (m^-1)
    % C     - fraction of plasma freq. after which shot noise is set to 0
    % 
    % Physical constants
    kB = 1.3806488e-23;             % Boltzmann constant
    me = 9.109383e-31;              % Mass of electron
    e = 1.60217657e-19;             % Electron charge
    eps0 = 8.85418782e-12;          % Vacuum Permittivity 

    % Plasma parameters
    vT = sqrt(2*kB*T/me);           % Thermal velocity of electrons
    LD = sqrt((kB*T*eps0)/(N*e^2)); % Debye length
    wp = vT/(sqrt(2)*LD);           % Plasma frequency

    % Antenna Response Function in integrand
    if dipole == 1
        F1 =@(x)((sinint(x))-(sinint(2*x)./2)-((2*((sin(x/2)).^4))./x))./x;
    else
        F1 =@(x)0.125*(x.^2 + 2*cos(x) - 2)./(x.^2);
    end

    % Faddeeva function:
    W = @(x) faddeeva(x);

    % Frequency vector
    f = logspace(log10(fmin),log10(fmax),n);
    w_vec = f*(2*pi); % Determine omega from f

    % Surface area based on whether spherical or cylindrical antenna
    if sph == 0
        S=2*pi*a*l;
    elseif sph == 1
        S=4*pi*a^2;
    else
        S = 0.01*6;
    end

    % Impact rate of electrons onto antenna
    N_e = ((4*pi)^(-0.5))*N*vT*S;

    % Vector for recording resistance
    V2QTN = zeros(size(w_vec));
    V2S = zeros(size(w_vec));
    Z_a = zeros(size(w_vec));

    % Progress bar initialization:
    textprogressbar('Progress: ');
    % For loops going through all frequency values in vector
    for i = 1:length(w_vec)

        % Display which frequency data point 
        %disp([int2str(i) ' out of ' int2str(length(w_vec)) ' complete.'])
        textprogressbar((i/length(w_vec))*100);
        % Permittivity
        w = w_vec(i);           % Set specific omega value
        z = @(x) w./(x.*vT); 
        DL = @(k) (1 + ((1+(1i)*sqrt(pi)*z(k).*W(z(k)))./((k.^2)*(LD^2))));

        % Integral calculation
        
        % Close to the plasma frequency the integral has a problem with
        % properly calculating the result. Therefore, to get good a good
        % result for the QTN, omega (w) values close to the plasma
        % frequency will have a different integral carried out compared to
        % the spectrum below that value. The shot noise is also only valid
        % for w << wp, so beyond the wp, the shot noise does not need to be
        % calculated.
        
        if w >= wp
            % Integral
            Int = @(k) (F1(k*l).*(besselj(0,k*a)).^2).*abs(imag(1./DL(k)));
            kmin = k_min;
            kmax = k_max;
            Integrand = integral(Int,kmin,kmax); 
            
            % Impedance
            Z = (4/(pi^2*eps0*w))*Integrand;
            Z_a(i) = Z;
            % Voltage power for QTN at frequency w
            V2QTN(i) = abs(4*kB*T*real(Z));
            % Shot noise expression is not valid at for w >~ wp
            V2S(i) = 0;
            
        else
            % Integral
            Int = @(k) (F1(k*l).*(besselj(0,k*a)).^2).*(1./(DL(k)));
            kmin = k_min;
            kmax = k_max;
            Integrand = integral(Int,kmin,kmax); 
            % Impedance
            Z = (4i/(pi^2*eps0*w))*Integrand;
            Z_a(i) = Z;

            % Voltage power for QTN at frequency w
            V2QTN(i) = abs(4*kB*T*real(Z));
            % For the shot noise save only shot noise values where w < C*wp 
            if w < C*wp
                V2S(i) = (2*(e^2)*N_e)*((abs(Z))^2)*exp(e*P/(kB*T));
            else
                V2S(i) = 0;
            end
            
        end
        clear Z
        
    end
    % End progress bar
    textprogressbar(' Complete');
    
end