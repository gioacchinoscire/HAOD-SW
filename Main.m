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
C_r=Solution.Cr;