clc
close all
global R_E e_E g_eq

%% Open Matlab pool to parallelize some operation..........................

% if isempty(gcp('nocreate'))==1;
%     parpool;
% end

%% Observer site geografic coordinates--------------------------------------

prompt={'Enter Geocentric Latitude [deg]:','Enter Longitude [deg]:',...
            'Enter Height from MSL [m]:'};
dlg_title='Observer Site Geographical data';
% def={'41.958055555555561',' 12.505277777777778','0'};
def={num2str(rad2deg(dms2rad(40,38,54.76580))),num2str(rad2deg(dms2rad(16,42,16.13192))),'529.627'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
obs_site=inputdlg(prompt,dlg_title,1,def,options);
 
Obs_data=zeros(3,1);

Obs_data(1)=rad2deg(dms2rad(40,38,54.76580));
Obs_data(2)=rad2deg(dms2rad(16,42,16.13192));
Obs_data(3)=str2double(obs_site{3});

%% Orbit determination Process.............................................

Solution=Orbit_Determination('GEOSC_Measures.txt',Obs_data);

%% Data Post-Processing....................................................


cost_fun=Solution.costfun;
Covariance=Solution.cov;
Initial_State=Solution.sol;

% % Tracklet identification--------------------------------------------------
% 
% format long
% B_s=1E-3;
% BC=1/(12.741621*B_s);
% prop_time=Time_JD*86400;
% options=odeset('AbsTol',1E-10,'RelTol',1E-10);
% 
% [~,Stato_sat]=ode113(@(t,x)Earth_dynamic(t,x,prop_time(1)/86400,order,Coeff,BC),(prop_time-prop_time(1)),Initial_State,options);
% 
% % observer position over time.---------------------------------------------
% 
% lat=Obs_data(1); % Latitudine [deg]
% lon=Obs_data(2); % Longitude [deg]
% alt=Obs_data(3); % Quota [m]
% 
% % Geoid Height
% 
% N_geo=WGS84_Height(lat,lon,70);
% % position vetor norm
% 
% r_d=((R_E*sqrt(1-e_E^2)/sqrt(1-e_E^2*cosd(lat)^2)))*cosd(lat);
% r_k=((R_E*sqrt(1-e_E^2)/sqrt(1-e_E^2*cosd(lat)^2)))*sind(lat);
% r_site=sqrt(r_d^2+r_k^2);
% 
% % Vettori Posizione nel riferimento ECEF ad ogni istante di misura
% 
% obs_pos=zeros(length(prop_time),3);
% 
% obs_vec=(r_site+(alt+N_geo)*1E-3)*[cosd(lat).*cosd(lon),cosd(lat).*sind(lon),...
%     sind(lat)]';
% 
% % vettore posizione nel riferimento ECI ad ogni istante di misura
% 
% for i=1:length(prop_time)
%     [year,mon,day,hr,min,sec] = invjday (prop_time(i)/86400);
%     [~, ~, jdut1, ~, ~, ~, ttt]= convtime ( year, mon, day, hr, min,...
%         sec, timezone(i), dut1(i), dat(i) );   % Time conversion
%     obs_pos(i,:)=ecef2eci(obs_vec/1E3,zeros(3,1),zeros(3,1),ttt,...
%         jdut1,lod(i),xp(i),yp(i),ast,dPsi(i),dEps(i));
% end
% 
% [Cel_Coord] =Telescope_measures(Stato_sat,obs_pos,prop_time);
% 
% % Transformation of the Covariance matrix in the RTN ref.frame-------------
% 
% R=Initial_State(1:3)'./norm(Initial_State); %Radial direction
% N=cross(R,Initial_State(4:6)')./norm(cross(R,Initial_State(4:6))); % cross-track direction
% T=cross(N,R); % along-track direction
% 
% gamma=[T;N;R]; % Transfrormation Matrix
% % gamma=eye(3);
% Covariance_RTN=zeros(3,3,length(cost_fun));
% clear i;
% 
% for i=1:length(cost_fun)+1
%     Covariance_RTN(:,:,i)=gamma*squeeze(Covariance(1:3,1:3,i))*gamma'/9;
% end
% % % Ellipsoid covariance plot------------------------------------------------
% 
% figure('Name','Position Covariace Ellipsoids');hold on
% 
% ellipsoid_color={'b','r','g'};
% 
% for j=1:3   
%     h=plot_gaussian_ellipsoid(zeros(3,1),(R_E/1E3)^2*squeeze(Covariance_RTN(1:3,1:3,end)),j,20);
%     set(h,'EdgeColor',ellipsoid_color{j});
%     set(h,'FaceColor','none');
% end
% 
% % % RA DEC Covariance--------------------------------------------------------
% % 
% % H=[-sind(Cel_Coord(1,1))/(norm(Stato_sat(1,1:2))),cosd(Cel_Coord(1,1))/(norm(Stato_sat(1,1:2))),0;-cosd(Cel_Coord(1,1))*sind(Cel_Coord(1,2))/(norm(Stato_sat(1,1:3))),-sind(Cel_Coord(1,1))*sind(Cel_Coord(1,2))/(norm(Stato_sat(1,1:3))),cosd(Cel_Coord(1,2))/(norm(Stato_sat(1,1:3)))]; 
% % Covariance_RaDEC=zeros(2,2,length(cost_fun));
% % clear i;
% % 
% % for i=1:length(cost_fun)
% %     Covariance_RaDEC(:,:,i)=(3600*180/pi)^2*H*squeeze(Covariance(1:3,1:3,i))*H';
% % end
% % 
% % 
% % figure('Name','Position Covariace Ellipsoids');hold on
% % 
% % ellipsoid_color={'b','r','g'};
% % 
% % for j=1:3   
% %     plot_gaussian_ellipsoid(zeros(2,1),R_E^2*squeeze(Covariance_RaDEC(1:2,1:2,end)),j,100);
% % end
% % % % Cost function plot-------------------------------------------------------
% % % 
% % figure('Name', 'Cost Function');
% % semilogy(cost_fun);
% % xlabel('Iteration number');
% % ylabel('Cost Function [deg.^2]'); grid on;
% % 
% % %% Read second Nigth data..................................................
% % Measure_files={'Measure_object_x.txt','Measure_object_y.txt','Measure_object_z.txt'};
% % clear i j k;
% % figure('Name', 'second night residual');hold on
% % 
% % for i=1:1
% %     
% %     % Import delle misure dal file di input
% %     
% %     Meas_data2=importdata(Measure_files{i},' ');
% %     % Matrice contenente le misure ottiche [RA (deg), DEC(deg)]
% %     
% %     Meas2=Meas_data2.data;
% %     
% %     % Vettore di celle contenente gli istanti temporali di misura espressi nel
% %     % formatio g/m/y h:m:s
% %     
% %     Time_str_2=Meas_data2.textdata;
% %     
% %     % Inizializzazione dei vettori per la scansione degli istanti temporali di
% %     % misura
% %     
% %     N_2=length(Meas2); %Numero di misure
% %     
% %     Time_UTC_2=zeros(N_2,6); % Vettore degli istanti temporali in UTC
% %     Time_JD_2=zeros(N_2,1);  % Vettore degli istanti temporali in julian date
% %     
% %     for j=1:N_2
% %         
% %         data_str=Time_str_2{j,1};      % Giorno
% %         day_str=Time_str_2{j,2};       % Istante del giorno
% %         
% %         Time_UTC_2(j,1)=round(str2double(data_str(1:4)));  % year
% %         Time_UTC_2(j,2)=round(str2double(data_str(6:7)));  % month
% %         Time_UTC_2(j,3)=round(str2double(data_str(9:10))); % day
% %         Time_UTC_2(j,4)=round(str2double(day_str(1:2)));   % hours
% %         Time_UTC_2(j,5)=round(str2double(day_str(4:5)));   % minutes
% %         Time_UTC_2(j,6)=str2double(day_str(7:end));        % seconds
% %         
% %         % Calcolo delle date juliane per ogni istante di misura
% %         
% %         Time_JD_2(j)=jday(Time_UTC_2(j,1),Time_UTC_2(j,2),Time_UTC_2(j,3),...
% %             Time_UTC_2(j,4),Time_UTC_2(j,5),Time_UTC_2(j,6));
% %     end
% %     Tot_Time_JD=[Time_JD;Time_JD_2];
% % 
% %     % Vettori Posizione nel riferimento ECEF ad ogni istante di misura
% %     
% %     N_time=length(Tot_Time_JD);
% %     obs_pos=zeros(N_time,3);
% %     
% %     obs_vec=r_site*[cosd(lat).*cosd(lon),cosd(lat).*sind(lon),...
% %         sind(lat)]';
% %     
% %     % vettore posizione nel riferimento ECI ad ogni istante di misura
% %     
% %     for k=1:length(Tot_Time_JD)
% %         [year,mon,day,hr,min,sec] = invjday (Tot_Time_JD(k));
% %         [~, ~, jdut1, ~, ~, ~, ttt]= convtime ( year, mon, day, hr, min,...
% %             sec, timezone, dut1, dat );   % Time conversion
% %         obs_pos(k,:)=ecef2eci(obs_vec,zeros(3,1),zeros(3,1),ttt,...
% %             jdut1,lod,xp,yp,ast,dPsi,dEps);
% %     end
% %     clear k
% %      for k=1:length(Meas2(:,1))
% %         [dRA_arcsec,dDEC_arcsec]=aberration_corr(deg2rad(Meas2(k,1)),deg2rad(Meas2(k,2)),obs_pos(k,:),Time_UTC_2(k,:),[dPsi dEps dut1]);
% % 
% %         Meas2(k,1)=Meas2(k,1)+dRA_arcsec/3600;
% %         Meas2(k,2)=Meas2(k,2)+dDEC_arcsec/3600;
% %      end
% % %     
% % 
% %     [Res_vec,Diff_vec]=Residual_vec(Initial_State',Tot_Time_JD,obs_pos,order,Coeff,[Meas;Meas2]);
% %     rgb=zeros(3,1);
% %     rgb(i)=1;
% %     
% %     subplot(3,1,1)
% %     
% %     hold on;
% %     plot(Time_JD_2*24-Time_JD(1)*24,(Diff_vec(end-N_2+1:end,1))*3600,'s-','Color',rgb);
% %     xlabel('Time: Hours from previous night');
% %     ylabel('RA residuals ['''']');box on;
% %     
% %     subplot(3,1,2)
% %     
% %     hold on;
% %     plot(Time_JD_2*24-Time_JD(1)*24,(Diff_vec(end-N_2+1:end,2))*3600,'s-','Color',rgb);
% %     xlabel('Time: Hours from previous night');
% %     ylabel('DEC residuals ['''']');box on;
% %     
% %     subplot(3,1,3)
% %     
% %     hold on;
% %     plot(Time_JD_2*24-Time_JD(1)*24,(Res_vec(end-N_2+1:end))*3600,'s-','Color',rgb);
% %     xlabel('Time: Hours from previous night');
% %     ylabel('Residuals norm ['''']');box on;
% % 
% %     rgb(i)=0;
% %     disp(mean((Res_vec(end-N_2+1:end))*3600));
% %     
% % end
% % 
% % 
