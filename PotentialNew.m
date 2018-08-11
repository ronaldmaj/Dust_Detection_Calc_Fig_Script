function [V] = PotentialNew(n_e, n_i, T_e, T_i, S, j_ph, height, cyl)

    % Function for finding the potential of a spacecraft due to electron and
    % ion thermal motion, photoemission and from direct impact of ions with
    % surface due to relative velocity of plasma with spacecraft

    % n_e: electron density (m^-3)
    % n_i: ion density (m^-3)
    % T_e: electron temperature (K)
    % T_i: ion temperature (K)
    % j_ph: photon impact rate (s^-1)
    % height: altitude above the Earth (m)
    % cyl: 1 = cylinder, 0 = cube 
    
    kB = 1.3806488e-23;     % Boltzmann constant
    me = 9.109383e-31;      % Mass of electron
    e = 1.60217657e-19;     % Electron charge
    % eps0 = 8.85418782e-12;  % Vacuum Permittivity
    mi = 1.67262e-27;       % Mass of proton
    Deg2Rad = @(x) (x/180)*pi;
    
    % Orbital velocity of spacecraft at height r:
    v_orb = @(r) sqrt(6.67408e-11*5.972e24/(6371e3+r));

    % Corrective factor for surface impacts
    if cyl == 1 % Cylinder
        A = 4/(pi*sqrt(2));
    elseif cyl == 0 % Cube
        A = 2/(3*pi);
    else
        error('Choose cyl either 0 or 1, i.e. 0 = cube, 1 = cylinder')
    end
    
    % Currents:
    I_e = @(S,n_e,T_e,V) e*S*n_e.*sqrt((kB*T_e)./(2*pi*me)).*exp((e*V)./...
        (kB*T_e));
    I_ph = @(S,j) (A*S/(2*sqrt(2)))*e*j;
    I_i = @(S,n_i,T_i,V) e*A*S*n_i.*sqrt((kB*T_i)/(2*pi*mi)).*...
        (2*sqrt(-e*V./(kB*T_i))+exp(-e*V./(kB*T_i)).*erfc(-e*V./(kB*T_i)));
    I_r = @(n_i, S, vorb,lat) e*n_i*A*S*sqrt(vorb^2 + (465.1*cos(lat))^2);

    % Vector of voltages used to find current flow potential
    V_vec = linspace(-2.5,0,1000);

    % Set up vector for final potential/s:
    [a,b] = size(n_e);
    V = zeros(a,b);

    for i=1:a
        for j=1:b
            
            if (j >10) || (j<27)
                % Calculate net current over all V_vec for (i,j)th
                % density/temp (excluding photoemission)
                R = I_e(S,n_e(i,j),T_e(i,j),V_vec) - I_i(S,n_i(i,j),...
                    T_i(i,j),V_vec) - I_r(n_i(i,j),S,v_orb(height),...
                    Deg2Rad((-100 + i*10)));
            else
                % Calculate net current over all V_vec for (i,j)th
                % density/temp (including photoemission)
                R = I_e(S,n_e(i,j),T_e(i,j),V_vec) - I_i(S,n_i(i,j),...
                    T_i(i,j),V_vec) - I_ph(S,j_ph) - I_r(n_i(i,j),...
                    S,v_orb(height),Deg2Rad((-100 + i*10)));                
            end
            
            % Create 2x1000 matrices to represent two lines; one is R, the
            % other is x=0;
            L1 = [V_vec; R];
            L2 = [V_vec; zeros(1,length(V_vec))];
            
            % Find intersection point in terms of (x,y) coordinate 
            P = InterX(L1,L2);
            
            % Save x-coordinate as the voltage
            V(i,j) = P(1);
            
            clear R
            
        end
    end



    