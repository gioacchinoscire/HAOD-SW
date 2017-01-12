function [da_SRP] =ERP_effect(r,v,Sun,JD,A_m,Nring,eps_m,ttt,jdut1,lod,xp,yp,ast,dPsi,dEps)
% Earth radiation pressure contribution to satellite accelration
global R_E c

r_norm=norm(r);
r_sun=norm(Sun);
R=r/r_norm;
N=cross(R,v./norm(v))./norm(cross(R,v./norm(v)));       % cross-track direction
T=cross(N,R);                                           % along-track direction
A_m=A_m/1E6;            % A/m ration in Km^2/kg

% Total number of ring elements............................................

Ner=6;                    % element number in each ring
Ntot=1;
for n=1:Nring
   Ntot=Ntot+n*Ner; 
end

csi=zeros(Nring,1);
gamma=zeros(Nring,1);
beta=zeros(Nring+1,1);
beta_c=zeros(Nring+1,1);

% Calculate pSR & MB/c.....................................................

year=invjday(JD);
JDap=jday(year,07,4,0,0,0);         % 4th july

Daph=(JD-JDap)*360/365.25;          % Daphelion in angle [deg]
SF=(1358/(1.004+0.0334*cosd(Daph)));    % W/m2-->[kg/s3]
pSR=SF/(c/1E3);                     % W/kgm2-->[kg /kms2]
MBc=pSR/4;                          % W/kgm2-->[1/s3]

% Limits of earth surface visible from sat.................................

limit_cs=asind((R_E/1E3)/r_norm);
beta_M=90-limit_cs;
beta(Nring+1)=beta_M;
% Limits of central cap....................................................

csi(1)=acosd((Ntot-1+cosd(limit_cs))/Ntot);     %[deg]
gamma(1)=asind(r_norm*sind(csi(1))/(R_E/1E3));  %[deg]
beta(1)=gamma(1)-csi(1);                        %[deg]

% Remaining ring boundaries & centers......................................

k=1;
for i=2:Nring
    k=k+6*(i-1);
    csi(i)=acosd(k*cosd(csi(1))-k+1);               %[deg]
    gamma(i)=asind(r_norm*sind(csi(i))/(R_E/1E3));  %[deg]
    beta(i)=gamma(i)-csi(i);                        %[deg]
end

for i=2:Nring+1
    beta_c(i)=(beta(i)+beta(i-1))/2;                %[deg]
end


% Radiation pressure due to central cap....................................

A_=2*(1-csi(1));
CTheta_sun=dot(r,Sun)/(r_norm*r_sun);
if(CTheta_sun<0)
    CTheta_sun=0;
end
tau=shadow_func(Sun,r);
r_ecef=eci2ecef(r,zeros(3,1),zeros(3,1),ttt,jdut1,lod,xp,yp,ast,dPsi,dEps); % ECI to ECEF conversion
lat_elem=asind(r_ecef(3)./norm(r_ecef));        %[deg]
[al,em]=EarthAlbedo_Emissivity(lat_elem,JD);
FluxOP=A_*(1+eps_m)*A_m*pSR;          % Optical flux    Km2/s3
FluxIR=A_*(1+eps_m)*A_m*MBc;          % IR flux         Km2/s3   

da_cap=(tau*al*FluxOP*CTheta_sun+em*FluxIR)*r/r_norm;

% Contribution of each ring................................................

da_tot=da_cap;
for i=2:Nring+1
   Nsi=Ner*(i-1);
  
   for l=1:Nsi
       lamd=l*360/Nsi;
       Relem=(R_E/1E3)*(cosd(beta_c(i-1))*R+sind(beta_c(i-1))*(cosd(lamd)*T+sind(lamd)*N));
       tau=shadow_func(Sun,Relem);
       Relem_ecef=eci2ecef(Relem,zeros(3,1),zeros(3,1),ttt,jdut1,lod,xp,yp,ast,dPsi,dEps); % ECI to ECEF conversion
       lat_elem=asind(Relem_ecef(3)./norm(Relem_ecef));
       Relm_v=Relem/(R_E/1E3);
       r_j=r-Relem;
       [al,em]=EarthAlbedo_Emissivity(lat_elem,JD);
       CTheta_sunj=dot(Relm_v,Sun)/(r_sun);
       if(CTheta_sunj<0)
           CTheta_sunj=0;
       end
       da_elem=(tau*al*FluxOP*CTheta_sunj+em*FluxIR)*r_j/(norm(r_j));
       da_tot=da_tot+da_elem;
   end

end    

da_SRP=da_tot;

end