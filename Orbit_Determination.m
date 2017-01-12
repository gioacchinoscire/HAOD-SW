function [Solution]=Orbit_Determination(Measures,Obs_data)
%----------------------------Orbit_Determination---------------------------%
% This function Perform the Orbit Determination
%   INPUT:
%           Obs_data:           Observer postion in georaphical coordinates
%           Measure:            File Name of the Measure input
%   OUTPUT:
%           Solution:           Structure containing the Solution data
%   NOTE:   This Function Need an Active Internet connection
%-------------------------------------------------------------------------%

global R_E f_E e_E w_E GM G c g_eq

format long
addpath('c:\mice\src\mice\');       
addpath('c:\mice\lib\');
cspice_furnsh('de432s.bsp')
order=input('Choose the Model degree (until 100): ');

% Measure Import...........................................................

Meas_data=importdata(Measures,' ');

Meas=Meas_data.data;        % Optical measure matrix [RA (deg), DEC(deg)]

N=length(Meas);             % Measure Number

Time_str=Meas_data.textdata;% Measure Time vector YYYY/MM/DD h:m:s

Time_UTC=zeros(N,6);        % Measure Times UTC

Time_JD=zeros(N,1);         % Measure Times JD

for i=1:N
    
    data_str=Time_str{i,1};                          % Day
    day_str=Time_str{i,2};                           % Day instant(h:m:s)
    Time_UTC(i,1)=round(str2double(data_str(1:4)));  % year
    Time_UTC(i,2)=round(str2double(data_str(6:7)));  % month
    Time_UTC(i,3)=round(str2double(data_str(9:10))); % day
    Time_UTC(i,4)=round(str2double(day_str(1:2)));   % hours
    Time_UTC(i,5)=round(str2double(day_str(4:5)));   % minutes
    Time_UTC(i,6)=str2double(day_str(7:end));        % seconds
       
    Time_JD(i)=jday(Time_UTC(i,1),Time_UTC(i,2),Time_UTC(i,3),...
        Time_UTC(i,4),Time_UTC(i,5),Time_UTC(i,6))+3/86400;  % JD calculation   
end

% WGS84 Parameters.........................................................

R_E=6378136.6;                          % Mean Equator radius       [m]
GM=398600.4415*1E9;                     % Earth gravity constant    [m3/s2]
G=6.67428*10^-11;                       % Constant of Gravitation   [m3/(kg*s2)]
f_E=1/298.25642;                        % Flattering
e_E=sqrt(f_E*(2-f_E));                  % Eccentricity
w_E=7.292115*10^-5;                     % Earth Angular Velocity    [rad/s]
c=299792458;                            % Speed of light [m/s]
g_eq=9.7803278;                                              % Gravity at the Equator    [m/s2]

% Earth Orientation Parameters.............................................

EGM=dlmread('EGM2008_to2190_TideFree.txt','',[0,0,5147,3]);
EOP=importdata('EOP C04 IAU2000A.txt','\n');
EOP=EOP(15:end);
DAT=importdata('tai-utc.txt','\n');

disp('-----------------------Orbit determination Data-----------------------');
fprintf('Observer latitude: %d deg\n', Obs_data(1));
fprintf('Observer longitude: %d deg\n', Obs_data(2));
fprintf('Observer height from MSL: %d m\n', Obs_data(3));
fprintf('Start time:  %i/%i/%i %i:%i:%5.3f\n',Time_UTC(1,1),Time_UTC(1,2),...
        Time_UTC(1,3),Time_UTC(1,4),Time_UTC(1,5),Time_UTC(1,6));
fprintf('Stop time:  %i/%i/%i %i:%i:%5.3f\n',Time_UTC(end,1),Time_UTC(end,2),...
        Time_UTC(end,3),Time_UTC(end,4),Time_UTC(end,5),Time_UTC(end,6));


% Geoid Height.............................................................

lat=Obs_data(1);                                   % Latitude        [deg]
long=Obs_data(2);                                  % Longitude       [deg]
alt=Obs_data(3);                                   % Height from MSL [m]
c_lat=cosd(lat);
s_lat=sind(lat);
r_d=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*c_lat^2)))*c_lat;
r_k=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*c_lat^2)))*s_lat;
r=sqrt(r_d^2+r_k^2);                               % Ellipsoid radius [m]

% Observer position over time..............................................

obs_pos=zeros(N,3);
recef=(r+alt)*1E-3.*[cosd(lat)*cosd(long);cosd(lat)*sind(long);sind(lat)];
timezone=0;

for i=1:length(Time_JD)
    [year,mon,day,hr,min,sec] = invjday (Time_JD(i));
    EOP_vector=find_EOP(Time_JD(i),EOP,DAT);
    xp=EOP_vector(1);
    yp=EOP_vector(2);
    dut1=EOP_vector(3);
    lod=EOP_vector(4);
    dX=EOP_vector(5);
    dY=EOP_vector(6);
    dat=EOP_vector(7);
    [~, ~, jdut1, ~, ~, ~, ttt]= convtime (year,mon,day,hr,min,sec,timezone,dut1,dat);
    pos=itrs2gcrs(recef,zeros(3,1),zeros(3,1),ttt,jdut1,lod,xp,yp,dX,dY); 
    obs_pos(i,:)=pos';
    [dRA_arcsec,dDEC_arcsec]=Aberration_Astrometry(Meas(i,1),Meas(i,2)',lat,long,alt,Time_JD(i),EOP,DAT);
    Meas(i,1)=Meas(i,1)+dRA_arcsec'/3600;
    Meas(i,2)=Meas(i,2)+dDEC_arcsec'/3600;
end

% Aberration Correction....................................................

clear i

% Initial Guess (SGP4).....................................................

TLE=importdata('TLE_EGEOS.txt',';');
lg1=(TLE{1,1}); % riga 1
lg2=(TLE{2,1}); % riga 2
Guess=(SGP4_propagator(lg1,lg2,Time_JD(1)))';

disp('Initial Guess:');

fprintf('\nX: %d Km\nY: %d Km\nZ: %d Km\nVx: %d Km/s\nVy: %d Km/s\nVz: %d Km/s\n',...
    Guess(1),Guess(2),Guess(3),Guess(4),Guess(5),Guess(6));
disp('----------------------------------------------------------------------');

% Orbit Determination (Powell's dogleg)....................................

kmax=70;                                           % Max Iteration Number
A_m=1;
Sol_struct=LSQ_PW_Optimization(Guess,kmax,obs_pos,Meas,Time_JD,order,EGM,EOP,DAT,A_m);

disp('----------------------------Iteration stop----------------------------');
disp('Solution State:');
fprintf('X: %d Km\nY: %d Km\nZ: %d Km\nVx: %d Km/s\nVy: %d Km/s\nVz: %d Km/s\n',...
    Sol_struct.sol(1),Sol_struct.sol(2),Sol_struct.sol(3),...
    Sol_struct.sol(4),Sol_struct.sol(5),Sol_struct.sol(6));
disp('----------------------------------------------------------------------');

Solution=Sol_struct;

end