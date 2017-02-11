
global R_E GM G f_E e_E w_E c g_eq

R_E=6378136.3;                          % Mean Equator radius       [m]
GM=398600.4415*1E9;                     % Earth gravity constant    [m3/s2]
G=6.67428*10^-11;                       % Constant of Gravitation   [m3/(kg*s2)]
f_E=1/298.25642;                        % Flattering
e_E=sqrt(f_E*(2-f_E));                  % Eccentricity
w_E=7.292115*10^-5;                     % Earth Angular Velocity    [rad/s]
c=299792458;                            % Speed of light [m/s]
g_eq=9.7803278;     

format long
addpath('c:\mice\src\mice\');       
addpath('c:\mice\lib\');
cspice_furnsh('de425.bsp')
order=input('Choose the Model degree (until 100): ');

EGM=dlmread('EGM2008_to2190_TideFree.txt','',[0,0,5147,3]);
EOP=importdata('EOP C04 IAU2000A.txt','\n');
EOP=EOP(15:end);
DAT=importdata('tai-utc.txt','\n');

x0=[7100.000125179327,1.584308459001052e-05,999.9991112263747,0,7.4,1];
x0=[42000,0,1300,1,1.3,1];
time=linspace(jday(2003,1,1,11,59,28),jday(2003,1,2,11,59,28));
options=odeset('AbsTol',1E-12,'RelTol',1E-12);
A_m=0;
[~,Stato_sat]=ode113(@(t,x)Earth_dynamic(t,x,time(1),order,EGM,EOP,DAT,A_m),(time*86400-time(1)*86400),x0',options);