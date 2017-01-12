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

Res_vec=sqrt(Diff_vec(:,1).^2+Diff_vec(:,2).^2);

end

