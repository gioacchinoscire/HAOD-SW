function [Delta_Cnm,Delta_Snm]=Ocean_Pole_Tide(n_max,m_max,m1,m2)
global GM R_E G g_eq w_E

% m1 and m2 in arcsec

OPT_coefficient=OPT_coeff();
[row,~]=find(OPT_coefficient(:,1)==n_max);
row=row(end);
n=OPT_coefficient(1:row,1);
OPT_coefficient=OPT_coefficient(1:row,3:6);

hlkCoeff=load_def_coeff();
hlkCoeff=hlkCoeff(1:n_max-1,3);

kn=ones(3,1)*hlkCoeff(1);
for i=2:n_max-1
    kn=[kn;ones(i+2,1)*hlkCoeff(i)];
end

m1=deg2rad(m1/3600);
m2=deg2rad(m2/3600);
rhow=1025;

R=(w_E^2*R_E^4/GM)*(4*pi*G*rhow/g_eq)*(1+kn)./(2*n+1);
Anm_R=OPT_coefficient(1:row,1);
Bnm_R=OPT_coefficient(1:row,2);
Anm_I=OPT_coefficient(1:row,3);
Bnm_I=OPT_coefficient(1:row,4);
gamma_r=0.6870;
gamma_i=0.0036;
Delta_Cnm=R.*(Anm_R.*(m1*gamma_r+m2*gamma_i)+Anm_I.*(m2*gamma_r-m1*gamma_i));
Delta_Snm=R.*(Bnm_R.*(m1*gamma_r+m2*gamma_i)+Bnm_I.*(m2*gamma_r-m1*gamma_i));

end