function [g]=harmonics_model(r_ecef,order,EGM,Sun,Moon,GM_sun,GM_moon,JD,xp,yp,ttt,jdut1)

global GM R_E

JD2010=jday(2010,01,01,00,00,00);   % Julian date at 2010/01/01 00:00:00 (limit for xp(t) and yp(t) IERS model)

%% Satellite data..........................................................

r_m=norm(r_ecef);       % vector position norm [m]
r_I=r_ecef(1);          % x ecef position component [m]
r_J=r_ecef(2);          % y ecef position component [m]
r_K=r_ecef(3);          % z ecef position component [m]
long=atan2(r_J,r_I);    % satellite longitude [rad]
lat=asin(r_K/r_m);      % satellite latitude [rad]

%% Sun & Moon data.........................................................

r_sun=norm(Sun);
lat_sun=asin(Sun(3)/r_sun);         % Sun latitude [rad]
long_sun=atan2(Sun(2),Sun(1));      % Sun longitude [rad]
r_moon=norm(Moon);
lat_moon=asin(Moon(3)/r_moon);      % Moon latitude [rad]
long_moon=atan2(Moon(2),Moon(1));   % Moon longitude [rad]

%% Fundamental values......................................................

n=EGM(:,1);
[row,col]=find(n==order);
row=row(end);
n=n(1:row);
m=EGM(1:row,2);
Cnm=EGM(1:row,3);
Snm=EGM(1:row,4);

if JD<JD2010
    coeff_x=[55.974,1.8243,0.18413,0.007024];
    coeff_y=[346.346,1.7896,-0.10729,-0.000908];
else
    coeff_x=[23.513,7.6141,0,0];
    coeff_y=[358.891,-0.6287,0,0];
end

for i=1:4
        xp_t=(ttt*100)^(i-1)*coeff_x(i)^(i-1);   % milli arcsec [mas]
        yp_t=(ttt*100)^(i-1)*coeff_y(i)^(i-1);   % milli arcsec [mas]
end

xp_t=deg2rad(xp_t/(1000*3600));
yp_t=deg2rad(yp_t/(1000*3600));
m1=rad2deg(xp-xp_t)*3600;
m2=-rad2deg(yp-yp_t)*3600;

Cnm(1)=-0.48416948*10^-3+11.6*10^-12*ttt*100;
Cnm(4)=0.9571612*10^-6+4.9*10^-12*ttt*100;
Cnm(8)=0.5399659*10^-6+4.7*10^-12*ttt*100;
Cnm(2)=sqrt(3)*xp_t*Cnm(1)-xp_t*Cnm(3)+yp_t*Snm(3);
Snm(2)=-sqrt(3)*yp_t*Cnm(1)-yp_t*Cnm(3)-xp_t*Snm(3);

low=1;
Leg=zeros(row,1);
Leg_plus1=zeros(row,1);
for i=2:order
    up=low+i;
    leg_fun=legendre(i,sin(lat));
    Leg(low:up)=leg_fun;
    Leg_plus1(low:up)=[leg_fun(2:end);0];
    low=up+1;
end

k=2*ones(size(n));
k(m==0)=1;

prod_nm=sqrt((factorial(n-m).*(2*n+1).*k)./factorial(n+m));

Leg=(-1).^(m).*Leg.*prod_nm;
n_mplus1=n-(m+1);
[row_min,col_min]=find(n_mplus1<0);
n_mplus1(row_min)=0;
prod_nm_plus1=sqrt((factorial(n_mplus1).*(2*n+1)*2)./factorial(n+m+1));
prod_nm_plus1(row_min)=0;
Leg_plus1=(-1).^(m+1).*Leg_plus1.*prod_nm_plus1;

%% Solid Earth Tides.......................................................
% Permanent Tide...........................................................

n_perm=n(1:7);
m_perm=m(1:7);
knm=Nominal_Love_numbers(n_perm,m_perm);
Leg_sun=[legendre(2,sin(lat_sun));legendre(3,sin(lat_sun))];
Leg_moon=[legendre(2,sin(lat_moon));legendre(3,sin(lat_moon))];

Delta_CS=2*knm./(2.*n_perm+1).*((GM_sun/GM)*(R_E/r_sun).^(n_perm+1).*Leg_sun.*exp(-1j*long_sun)+(GM_moon/GM)*(R_E/r_moon).^(n_perm+1).*Leg_moon.*exp(-1j*long_moon));
Delta_CS(1)=Delta_CS(1)-(4.4228*10^-8)*(-0.31460)*knm(1);
Delta_Cnm_PST=real(Delta_CS);
Delta_Snm_PST=imag(Delta_CS);

n_perm=n(1:3);
m_perm=m(1:3);
[~,knm_p]=Nominal_Love_numbers(n_perm,m_perm);

Delta_CS=(knm_p./5).*((GM_sun/GM)*(R_E/r_sun).^3.*Leg_sun(1:3).*exp(-1j*long_sun)+(GM_moon/GM)*(R_E/r_moon).^3.*Leg_moon(1:3).*exp(-1j*long_moon));
Delta_Cnm_PST=[Delta_Cnm_PST;real(Delta_CS)];
Delta_Snm_PST=[Delta_Snm_PST;imag(Delta_CS)];

Cnm(1:10)=Cnm(1:10)+Delta_Cnm_PST;
Snm(1:10)=Snm(1:10)+Delta_Snm_PST;

% Frequency correction.....................................................

[GMST]=Greenswich_MST(ttt,jdut1);

Delta_Freq_0=Solid_Tide_Frequency_correction(2,0,GMST,ttt);
Delta_Freq_1=Solid_Tide_Frequency_correction(2,1,GMST,ttt);
Delta_Freq_2=Solid_Tide_Frequency_correction(2,2,GMST,ttt);
Delta_Freq=[Delta_Freq_0;Delta_Freq_1;Delta_Freq_2];

Delta_Freq_Cnm=real(Delta_Freq);
Delta_Freq_Snm=imag(Delta_Freq);

Cnm(1:3)=Cnm(1:3)+Delta_Freq_Cnm;
Snm(1:3)=Snm(1:3)+Delta_Freq_Snm;

%% Solid Earth Pole Tide...................................................

Delta_Cnm_EPT=-1.333*10^-9*(m1+0.0115*m2);
Delta_Snm_EPT=-1.333*10^-9*(m2-0.0115*m1);

Cnm(2)=Cnm(2)+Delta_Cnm_EPT;
Snm(2)=Snm(2)+Delta_Snm_EPT;

%% Ocean Earth Pole Tide...................................................

[Delta_Cnm_OPT,Delta_Snm_OPT]=Ocean_Pole_Tide(order,order,m1,m2);

Cnm=Cnm+Delta_Cnm_OPT;
Snm=Snm+Delta_Snm_OPT;

%% Gravity field in ITRS Frame.............................................

dU_dr=-(GM/r_m^2)*(sum(((R_E/r_m).^n).*(n+1).*Leg.*(Cnm.*cos(m.*long)+Snm.*sin(m.*long)))+1);
dU_dlong=(GM/r_m)*sum((R_E/r_m).^n.*Leg.*m.*(Snm.*cos(m.*long)-Cnm.*sin(m.*long)));
dU_dlat=(GM/r_m)*sum((R_E/r_m).^n.*((Cnm.*cos(m.*long)+Snm.*sin(m.*long))).*(Leg_plus1-m.*tan(lat).*Leg));
            
g=zeros(3,1);

g(1)=((1/r_m)*dU_dr-(r_K/(r_m^2*sqrt(r_I^2+r_J^2)))*dU_dlat)*r_I-(1/(r_I^2+r_J^2)*dU_dlong)*r_J;
g(2)=((1/r_m)*dU_dr-(r_K/(r_m^2*sqrt(r_I^2+r_J^2)))*dU_dlat)*r_J+(1/(r_I^2+r_J^2)*dU_dlong)*r_I;
g(3)=(1/r_m)*dU_dr*r_K+(sqrt(r_I^2+r_J^2)/r_m^2)*dU_dlat;

end

