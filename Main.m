clc
close all


%% Open Matlab pool to parallelize some operation..........................

% if isempty(gcp('nocreate'))==1;
%     parpool;
% end

%% Observer site geografic coordinates-------------------------------------

prompt={'Enter Geocentric Latitude [deg]:','Enter Longitude [deg]:',...
            'Enter Height from MSL [m]:'};
dlg_title='Observer Site Geographical data';
% def={'41.958055555555561',' 12.505277777777778','0'};
def={num2str(rad2deg(dms2rad(40,38,54.76580))),num2str(rad2deg(dms2rad(16,42,16.13192))),'529.627'};
def={num2str(rad2deg(dms2rad(50,38,18))),num2str(rad2deg(dms2rad(13,50,48.3))),'275'};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
obs_site=inputdlg(prompt,dlg_title,1,def,options);
 
Obs_data=zeros(3,1);

Obs_data(1)=str2double(obs_site{1});
Obs_data(2)=str2double(obs_site{2});
Obs_data(3)=str2double(obs_site{3});

%% Orbit determination Process.............................................

TLE_file='TLE_NAVSTAR.txt';
Measure_file='TDM_Measures_NAVSTAR 43 (USA 132)_old.txt';

Solution=Orbit_Determination(Measure_file,Obs_data,TLE_file);

%% Data Post-Processing....................................................

cost_fun=Solution.costfun;
Initial_State=Solution.sol;
C_r=Solution.Cr;

lat=str2double(obs_site{1});
long=str2double(obs_site{2});
alt=str2double(obs_site{3});

% Measure Import...........................................................

Meas_data=importdata(Measure_file,' ');

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
        Time_UTC(i,4),Time_UTC(i,5),Time_UTC(i,6));  % JD calculation   
end

global R_E f_E e_E w_E GM G c g_eq

EGM=dlmread('EGM2008_to2190_TideFree.txt','',[0,0,5147,3]);
EOP=importdata('EOP C04 IAU2000A.txt','\n');
EOP=EOP(15:end);
DAT=importdata('tai-utc.txt','\n');


R_E=6378136.3;                          % Mean Equator radius       [m]
GM=398600.4415*1E9;                     % Earth gravity constant    [m3/s2]
G=6.67428*10^-11;                       % Constant of Gravitation   [m3/(kg*s2)]
f_E=1/298.25642;                        % Flattering
e_E=sqrt(f_E*(2-f_E));                  % Eccentricity
w_E=7.292115*10^-5;                     % Earth Angular Velocity    [rad/s]
c=299792458;                            % Speed of light [m/s]
g_eq=9.7803278;                                              % Gravity at the Equator    [m/s2]

obs_pos=zeros(N,3);
timezone=0;
[x,y,z]=lla2ecef_deg(lat,long,alt);
recef=[x;y;z]./1E3;
dut1_vec=zeros(1,length(Time_JD));
dat_vec=zeros(1,length(Time_JD));
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
    dut1_vec(i)=dut1;
    dat_vec(i)=dat;
    [~, ~, jdut1, ~, ~, ~, ttt]= convtime (year,mon,day,hr,min,sec,timezone,dut1,dat);
    pos=itrs2gcrs(recef,zeros(3,1),zeros(3,1),ttt,jdut1,lod,xp,yp,dX,dY); 
    obs_pos(i,:)=pos';
end


% Aberration Correction....................................................

[dRA,dDEC]=Aberration_Astrometry(Time_JD,Meas(:,1),Meas(:,2),dut1_vec,dat_vec);
Meas(:,1)=Meas(:,1)+dRA;
Meas(:,2)=Meas(:,2)+dDEC;
options=odeset('AbsTol',1E-12,'RelTol',1E-12);

order=30;
[~,Stato_sat]=ode113(@(t,x)Earth_dynamic(t,x,Time_JD(1),order,EGM,EOP,DAT,0),(Time_JD*86400-Time_JD(1)*86400),Initial_State,options);

% Simulated measurements...................................................

Cel_Coord_model=Telescope_measures(Stato_sat,obs_pos,1);

% Covariance Estimation....................................................
    
H=observ_grad(Initial_State,Time_JD,obs_pos,order,EGM,EOP,DAT,Meas,0);
R_meas=(1/3600)^2*eye(numel(Meas));
Covariance=inv(H'*(inv(R_meas))*H);

[t_tran,Stato_for_tran]=ode113(@(t,x)Earth_dynamic(t,x,Time_JD(1),order,EGM,EOP,DAT,0),0:300:(Time_JD(end)*86400-Time_JD(1)*86400),Initial_State,options);

options=odeset('AbsTol',1E-6,'RelTol',1E-6);
[~,Transition_Mat]=ode113(@(t,x)StateTransition_dynamic(t,x,order,EGM,EOP,DAT,0,Time_JD(1),[t_tran,Stato_for_tran]),Time_JD*86400-Time_JD(1)*86400,reshape(eye(6),[36,1]),options);


Covariance_time=zeros(6,6,length(Time_JD));
Covariance_RTN=zeros(3,3,length(Time_JD));
Covariance_time(:,:,1)=Covariance(:,:,end);
first=squeeze(Covariance_time(:,:,1));
for i=1:length(Time_JD)
    last_trans=reshape(Transition_Mat(i,:),[6,6]);
    cov_i=(last_trans')*first*last_trans;
    Covariance_time(:,:,i)=cov_i;
    R_vec=Stato_sat(i,1:3)'/norm(Stato_sat(i,1:3));
    N_vec=cross(R_vec,Stato_sat(i,4:6)')/norm(cross(R_vec,Stato_sat(i,4:6)));
    T_vec=cross(N_vec,R_vec);
    Rot_matrix=[R_vec,N_vec,T_vec];
    Covariance_RTN(:,:,i)=Rot_matrix'*cov_i(1:3,1:3)*Rot_matrix;
end

%% Plot Results------------------------------------------------------------

% Cost function............................................................

figure('Name','Cost Function')
semilogy(cost_fun)
xlabel('Iteration number')
ylabel('Cost function [deg^2]')

% Residuals plot...........................................................

figure('Name','Residuals')
plot(3600*(Cel_Coord_model(:,1)-Meas(:,1)),3600*(Cel_Coord_model(:,2)-Meas(:,2)),'*');
hold on
plot(3600*mean(Cel_Coord_model(:,1)-Meas(:,1)),3600*mean(Cel_Coord_model(:,2)-Meas(:,2)),'sr','MarkerFaceColor','r');
RA_std=std(Cel_Coord_model(:,1)-Meas(:,1)-mean(Cel_Coord_model(:,1)-Meas(:,1)))*3600;
DEC_std=std(Cel_Coord_model(:,2)-Meas(:,2)-mean(Cel_Coord_model(:,2)-Meas(:,2)))*3600;
angle=linspace(0,2*pi);
plot(RA_std*cos(angle)+3600*mean(Cel_Coord_model(:,1)-Meas(:,1)),DEC_std*sin(angle)+3600*mean(Cel_Coord_model(:,2)-Meas(:,2)),'k');
plot(3*RA_std*cos(angle)+3600*mean(Cel_Coord_model(:,1)-Meas(:,1)),3*DEC_std*sin(angle)+3600*mean(Cel_Coord_model(:,2)-Meas(:,2)),'g');
grid on
xlabel('Right Ascension Residual \DeltaRA [arc sec]');
ylabel('Declination Residual \DeltaDEC [arc sec]');
axis equal

% Covariance plot..........................................................

figure
subplot(1,2,1)
Pos_1=plot_gaussian_ellipsoid(zeros(3,1), Covariance_time(1:3,1:3,1)*R_E^2, 1, 30);
set(Pos_1,'EdgeColor','b')
view(35,20)
grid on
title('Position Covariance ellipsoids')
hold on
Pos_3=plot_gaussian_ellipsoid(zeros(3,1), Covariance_time(1:3,1:3,1)*R_E^2, 3, 30);
view(35,20)
grid on
set(Pos_3,'EdgeColor','g')
subplot(1,2,2)
title('Velocity Covariance ellipsoids')
Vel_1=plot_gaussian_ellipsoid(zeros(3,1), Covariance_time(4:6,4:6,1)*sqrt(GM/R_E)^2, 1, 30);
view(35,20)
grid on
set(Vel_1,'EdgeColor','b')
hold on
Vel_3=plot_gaussian_ellipsoid(zeros(3,1), Covariance_time(4:6,4:6,1)*sqrt(GM/R_E)^2, 3, 30);
view(35,20)
grid on
set(Vel_3,'EdgeColor','g')