function [dRA_arcsec,dDEC_arcsec]=Aberration_Astrometry(RA,DEC,lat,long,H_alt,JD,EOP,DAT)
%---------------------------Aberration_Astrometry--------------------------%
% This function is used with the NASA Matlab Toolbox for the SPK format
% Ephemeris
%   INPUT:
%           RA:             Right Ascension         [deg][1xN]            
%           DEC:            Declination             [deg][1xN]
%           lat:            latitude                [deg]            
%           long:           longitude               [deg]
%           H_alt:          Geoig Height+Altitude from geoid
%           JD:             Julian Data             [deg][1xN]
%   OUTPUT:
%           dRA_arcsec:     Delta Right Ascension   [arcsec][1xN]
%           dDEC_arcsec:    Delta Declination       [arcsec][1xN]
%   NOTE:   use the function  cspice_furnsh(SPKFILENAME) before This Func
%-------------------------------------------------------------------------%

global R_E e_E c

size_RA=size(RA);
size_DEC=size(DEC);
size_JD=size(JD);

if((size_RA(1)~=1 && size_DEC(1)~=1 && size_JD(1)~=1) && (size_RA(2)==1 || size_DEC(2)==1 || size_JD(2)==1 || size_RA(2)~=size_DEC(2)|| size_RA(2)~=size_JD(2)))
   error('Input Vectors shall be row vectors 1xN of equal lenghts') 
end

r_d=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*cosd(lat)^2)))*cosd(lat);
r_k=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*cosd(lat)^2)))*sind(lat);
r=sqrt(r_d^2+r_k^2);                                            % Ellipsoid radius [m]

recef=(r+H_alt)*1E-3.*[cosd(lat)*cosd(long);cosd(lat)*sind(long);sind(lat)];
timezone=0;
[year,mon,day,hr,min,sec] = invjday (JD);
[EOP_vector]=find_EOP(JD,EOP,DAT);
xp=EOP_vector(1);
yp=EOP_vector(2);
dut1=EOP_vector(3);
lod=EOP_vector(4);
dX=EOP_vector(5);
dY=EOP_vector(6);
dat=EOP_vector(7);
[~, ~, jdut1, ~, ~, ~, ttt]= convtime (year,mon,day,hr,min,sec,timezone,dut1,dat);
[~,veci]=itrs2gcrs(recef,zeros(3,1),zeros(3,1),ttt,jdut1,lod,xp,yp,dX,dY);

vers_stelle_top_mis=[cosd(RA).*cosd(DEC);cosd(DEC).*sind(RA);sind(DEC)];

Earth_from_Sun=Planet_Ephemeris(JD,'EARTH','SUN',dut1,dat);

b=(Earth_from_Sun(:,4:6)'+veci)/(c/1E3);

u=vers_stelle_top_mis;

v_new=u+b;
norm_tot=ones(3,1)*sqrt(sum(v_new.^2));
v_new=(v_new./norm_tot);

RA_new=atan2(v_new(2,:),v_new(1,:));
RA_new=RA_new+pi*(1-sign(RA_new));
RA_new=rad2deg(RA_new);
DEC_new=asind(v_new(3,:));

dRA_arcsec=(RA_new-RA)*3600;
dDEC_arcsec=(DEC_new-DEC)*3600;

end
