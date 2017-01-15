function [al,em]=EarthAlbedo_Emissivity(lat,JD)

% Earth Albedo and Emissivity model as in Knocke et.al 

global w_E

JD0=jday(1981,12,22,0,0,0);
P1=legendre(1,sind(lat));
P2=legendre(2,sind(lat));
a0=0.34;
a2=0.29;
c0=0;
c1=0.1;
c2=0;
a1=c0+c1*cos(w_E*(JD-JD0))+c2*sin(w_E*(JD-JD0));
al=a0+a1*P1+a2*P2;

e0=0.68;
e2=-0.18;
k0=0;
k1=-0.07;
k2=0;
e1=k0+k1*cos(w_E*(JD-JD0))+k2*sin(w_E*(JD-JD0));
em=e0+e1*P1+e2*P2;

end