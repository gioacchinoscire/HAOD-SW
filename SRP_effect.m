function [da_SRP] =SRP_effect(r,Sun,JD,C_r)
% solar radiation pressure contribution to satellite accelaration
global c
year=invjday(JD);
JDap=jday(year,07,4,0,0,0);         % 4th july

Daph=(JD-JDap)*360/365.25;          % Daphelion in angle
SF=1358/(1.004+0.0334*cosd(Daph));
pSR=SF/(c/1E3);
A_m=0.001;
Sun_sat=Sun-r;
v=shadow_func(Sun,r);
da_SRP=-v*pSR*C_r*abs(A_m)*Sun_sat/(norm(Sun_sat)^2);

end

