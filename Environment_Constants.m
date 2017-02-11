%% Phisical constants-------------------------------------------------------

c=299792458;                            % Speed of light            [m/s]
AU=149597870.691000015;                 % Astronomical Unit         [Km]
GM_moon=(GM/1E9)*0.012300036905526;     % Moon Gravity constant     [Km3/s2]
GM_sun=1.327124400410210*10^11;         % Sun Gravity constant      [km3/s2]
SRS=1.97412574336e-8*AU;                % Schwarzschild Sun Radius  [Km]

%% EGM2008 Parameters (woring in TT)---------------------------------------

R_E=6378136.3;                          % Mean Equator radius       [m]
GM=398600.4415*1E9;                     % Earth gravity constant    [m3/s2]
G=6.67428*10^-11;                       % Constant of Gravitation   [m3/(kg*s2)]

%% WGS84 Parameters--------------------------------------------------------

a_wgs=6378137;                          % Mean Equatorial radius    [m]
w_E=7.292115*10^-5;                     % Earth Angular Velocity    [rad/s]
f_E=1/298.25642;                        % Flattering
e_E=sqrt(f_E*(2-f_E));                  % Eccentricity
g_eq=9.7803278;                         % Gravity at the Equator    [m/s2]
A=8.0091029*10^37;
B=8.0092559*10^37;
C=8.0354872*10^37;
Mass=5.9733328*10^24;                   % Earth Mass                [Kg]

%% Time Constants----------------------------------------------------------

DJ2000=2451545.0;                       % J2000 Reference TT        [days]