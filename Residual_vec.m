function [Res_vec,Diff_vec,Cel_Coord_model]=Residual_vec(x0,time,obs_pos,order,EGM,EOP,DAT,Meas,A_m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the residual between the estimated dynamic 
% and the measures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% X0: satellite orbital dynamic at start ime
% obs_pos: observer inertial position
% time: measures time
% OUTPUT
% Res_vec: vector of the celestial coordinates residual at each measure
% time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Target dynamics..........................................................

options=odeset('AbsTol',1E-12,'RelTol',1E-12);

[~,Stato_sat]=ode113(@(t,x)Earth_dynamic(t,x,time(1),order,EGM,EOP,DAT,A_m),(time*86400-time(1)*86400),x0,options);

% Simulated measurements...................................................

Cel_Coord_model=Telescope_measures(Stato_sat,obs_pos,1);

% Fitness Function evaluation..............................................

Diff_vec=Meas-Cel_Coord_model; % (Nobs,2)

Diff_vec=rad2deg(atan2(sind(Diff_vec),cosd(Diff_vec)));

vers_meas=[cosd(Meas(:,2)).*cosd(Meas(:,1)),cosd(Meas(:,2)).*sind(Meas(:,1)),sind(Meas(:,2))];
vers_calc=[cosd(Cel_Coord_model(:,2)).*cosd(Cel_Coord_model(:,1)),cosd(Cel_Coord_model(:,2)).*sind(Cel_Coord_model(:,1)),sind(Cel_Coord_model(:,2))];

Res_vec=acosd(dot(vers_meas',vers_calc')');
Res_vec(Res_vec<0)=Res_vec(Res_vec<0)+360;

% Res_vec=sqrt(Diff_vec(:,1).^2+Diff_vec(:,2).^2);

end

