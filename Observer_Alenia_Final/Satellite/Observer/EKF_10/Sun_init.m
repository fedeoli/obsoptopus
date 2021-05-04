%%%%%%%%%%%%%%%%%%%%%%% SUN SENSOR INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine Sun Vec in Inertial Frame - see https://en.m.wikipedia.org/wiki/Position_of_the_Sun
    global Agent ObserverTest
    
    %%%% ECLIPTIC COORDINATES %%%%
    % # days since greenwich noon 01/01/2000
    T = 7551; % 3 settembre 2020
    CONST.ua = 1.496e+08;

    % mean anomaly (ref 4.1)
    Msun=357.528+0.9856003*T;  
    
    %ecliptic longitude (ref 4.2)
    Vsun=280.461+0.9856474*T+1.915*sind(Msun)+0.020*sind(2*Msun); 
    
    % ecliptic latitude
    Lsun = 0; 
    
    %distance sun-Earth in ua
    R = 1.000014-0.01671*cosd(Msun)-0.00014*cosd(2*Msun); 
    
    %%%% EQUATORIAL COORDINATES %%%%
    % ref 4.3
    obliquity=23.4393-0.0000004*T;
    
    %%% Conversion in equatorial coordinates %%%
    right_asc = atan2(cos(obliquity)*sind(Vsun), cosd(Vsun));
    decl = asin(sind(obliquity)*sind(Vsun));
    
    %%% Conversion in Rectangular equatorial coordinates %%%
    % ref 4.4
    Xsun=R*cosd(Vsun);
    Ysun=R*cosd(obliquity)*sind(Vsun);
    Zsun=R*sind(obliquity)*sind(Vsun);
    Si=[Xsun Ysun Zsun];
    Si=Si./norm(Si);
   
    
    params.Si = Si;
    

    % symbolic definition - sun position in ECI
    % Volt = sym('V', [1 3]);
    % Ncell = sym('N', [1 3]);
    Psun = sym('P', [1 3]);
    params.Psun = Psun;

    % symbolic - sun position in Body (ref 4.5)
    hss = dcm*transpose(Psun);
%     hss = transpose(Psun);
    
    %%%%%%%%%%%%%%%%%%%%%%% ALBEDO EFFECT INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %data from "Coarse Sun Sensing for Attitude Determination of a CubeSat"
    % John Gaebler Summa Cum Laude Thesis Bachelor of Science in Aerospace Engineering
    
    % Area of the cell - albedo model (m^2)
    A=30.16/100^2;
    
    % reflectivity of the cell - max value
    eff=0.251; 
    
    % W/m^2 Sun incidence energy avg - EM0 in Bhanderi
    K=1367; 
    
    % Max theoretical reflected intensity - sun normal to sensor (0 deg)
    params.I0 = A*K*cos(0)*eff; % ref 4.6

    % dcm matrix from ECI to ECI Fixed
    dcmecef = dcmeci2ecef('IAU-2000/2006',[2019 12 15 10 20 36]);
    params.dcmecef = dcmecef;
    
    % to ECEF in Km;
    for k = 1:ObserverTest.Nagents
        Agent(k).sat_ecef = dcmecef*satellites_iner_ECI(1+6*(k-1):3+6*(k-1))*1E03;
    end
    params.sun_ecef = params.dcmecef*params.Si'*CONST.ua;