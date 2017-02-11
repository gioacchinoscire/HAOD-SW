function [Cel_Coord] =Telescope_measures(Stato_sat,obs_pos,Vel_corr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulate the measures that we can obtain from an optical
% sensor without errors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Stato_sat: satellite orbital dynamic at the visibility instant
% obs_pos: observer inertial position
% OUTPUT
% Cel_Coord: Celestial coordinates obtained from the observation. These
% measures are reported in a matrix with N_obs row and 2 columw. The first
% column contains all the RA values, while the second all the DEC values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c

rho=Stato_sat(:,1:3)-obs_pos;  % Relative position in J2000

normRHO=sqrt(rho(:,1).^2+rho(:,2).^2+rho(:,3).^2);           % Relative range 

Cel_Coord=[atan2(rho(:,2),rho(:,1)),asin(rho(:,3)./normRHO)];                 % RA & DEC

Cel_Coord(Cel_Coord(:,1)<0,1)=Cel_Coord(Cel_Coord(:,1)<0,1)+2*pi;

Cel_Coord=rad2deg(Cel_Coord);

% satellite velocity effect................................................

if(Vel_corr==0)
    unc_vers=[cosd(Cel_Coord(:,2)).*cosd(Cel_Coord(:,1)),cosd(Cel_Coord(:,2)).*sind(Cel_Coord(:,1)),sind(Cel_Coord(:,2))];
    
    corr_vers=unc_vers+Stato_sat(:,4:6)./(c/1E3);

    corr_vers=corr_vers./(sqrt(dot(corr_vers',corr_vers'))'*ones(1,3));

    corr_RA=rad2deg(atan2(corr_vers(:,2),corr_vers(:,1)));
    corr_DEC=asind(corr_vers(:,3));
    
    Cel_Coord=[corr_RA,corr_DEC];
    
    Cel_Coord(Cel_Coord(:,1)<0,1)=Cel_Coord(Cel_Coord(:,1)<0,1)+360;
end

end

