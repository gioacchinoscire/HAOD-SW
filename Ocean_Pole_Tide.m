function [Delta_Cnm,Delta_Snm]=Ocean_Pole_Tide(n,m,m1,m2,ABnmCoeff,hlkCoeff)
global GM R_E G g_eq w_E

% m1 and m2 in arcsec

[row,~]=find(ABnmCoeff(:,1)==n & ABnmCoeff(:,2)==m);
m1=deg2rad(m1/3600);
m2=deg2rad(m2/3600);
rhow=1025;
kn=hlkCoeff(n,3);
R=(w_E^2*R_E^4/GM)*(4*pi*G*rhow/g_eq)*(1+kn)/(2*n+1);
Anm_R=ABnmCoeff(row,3);
Bnm_R=ABnmCoeff(row,4);
Anm_I=ABnmCoeff(row,5);
Bnm_I=ABnmCoeff(row,6);
gamma_r=0.6870;
gamma_i=0.0036;
Delta_Cnm=R*(Anm_R*(m1*gamma_r+m2*gamma_i)+Anm_I*(m2*gamma_r-m1*gamma_i));
Delta_Snm=R*(Bnm_R*(m1*gamma_r+m2*gamma_i)+Bnm_I*(m2*gamma_r-m1*gamma_i));

end