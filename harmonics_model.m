function [g]=harmonics_model(r_ecef,order,EGM,Sun,Moon,GM_sun,GM_moon,JD,xp,yp,ttt,jdut1)

global GM R_E
% ABnm=Global_Coeff.OceanPoleTide;
% hlk=Global_Coeff.LoveShidaNum;
% K20_corr=Global_Coeff.K20corr;
% K21_corr=Global_Coeff.K21corr;
% OT_Coeff=Global_Coeff.OceanTideCoeff;
% Coeff=Global_Coeff.EGM;
% Global_k=zeros(length([K20_corr(:,1);K21_corr(:,1)]),17);
% Global_k(:,1)=[K20_corr(:,1);K21_corr(:,2)];
% Global_k(:,2)=[K20_corr(:,2);K21_corr(:,1)];
% Global_k(:,3:end)=[K20_corr(:,3:end);K21_corr(:,3:end)];
% Doodson=OT_Coeff.textdata(:,1);
% OT_data=zeros(length(OT_Coeff.data(:,1)),7);
% for i=1:length(OT_Coeff.data(:,1))
%     OT_data(i,1)=str2double(Doodson{i});
% end
% OT_data(:,2:7)=OT_Coeff.data;

JD2010=jday(2010,01,01,00,00,00);   % Julian date at 2010/01/01 00:00:00 (limit for xp(t) and yp(t) IERS model)

% Satellite data...........................................................

r_m=norm(r_ecef);       % vector position norm [m]
r_I=r_ecef(1);          % x ecef position component [m]
r_J=r_ecef(2);          % y ecef position component [m]
r_K=r_ecef(3);          % z ecef position component [m]
long=atan2(r_J,r_I);    % satellite longitude [rad]
lat=asin(r_K/r_m);      % satellite latitude [rad]

% r_sun=norm(Sun);
% lat_sun=asin(Sun(3)/r_sun);      % Sun latitude [rad]
% long_sun=atan2(Sun(2),Sun(1));
% r_moon=norm(Moon);
% lat_moon=asin(Moon(3)/r_moon);      % Sun latitude [rad]
% long_moon=atan2(Moon(2),Moon(1));

% Gravity field in ITRF Frame..............................................

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


dU_dr=-(GM/r_m^2)*(sum(((R_E/r_m).^n).*(n+1).*Leg.*(Cnm.*cos(m.*long)+Snm.*sin(m.*long)))+1);
dU_dlong=(GM/r_m)*sum((R_E/r_m).^n.*Leg.*m.*(Snm.*cos(m.*long)-Cnm.*sin(m.*long)));
dU_dlat=(GM/r_m)*sum((R_E/r_m).^n.*((Cnm.*cos(m.*long)+Snm.*sin(m.*long))).*(Leg_plus1-m.*tan(lat).*Leg));
            
            
% dU_dr=1;                % Geopotential radius-Partial Derivative 
% dU_dlat=0;              % Geopotential latitude-Partial Derivative 
% dU_dlong=0;             % Geopotential longitude-Partial Derivative 
% GMTrad = gstime(jdut1);
% for n=2:order
% 
%     dU_dr_p=0;
%     dU_dlat_p=0;
%     dU_dlong_p=0;
%     
%     Leg=legendre(n,sin(lat));
%     if n<=3
%         Leg_sun=legendre(n,sin(lat_sun));
%         Leg_moon=legendre(n,sin(lat_moon));
%     end
%     for m=0:n
%         if m==0
%             k=1;
%         else
%             k=2;
%         end
%         prod_nm=sqrt((factorial(n-m)*(2*n+1)*k)/factorial(n+m));
% 
%         if(n>=2)
%             [row,~]=find(Coeff(:,1)==n & Coeff(:,2)==m);
%             Cnm=Coeff(row,3);
%             Snm=Coeff(row,4);
%         else
%             Cnm=0;
%             Snm=0;
%         end
%         
%         % Values different from EGM 2008 Model.............................
%         if JD<JD2010
%                 coeff_x=[55.974,1.8243,0.18413,0.007024];
%                 coeff_y=[346.346,1.7896,-0.10729,-0.000908];
%         else
%             coeff_x=[23.513,7.6141,0,0];
%             coeff_y=[358.891,-0.6287,0,0];
%         end
%         for i=1:4
%                 xp_t=((JD-JD2000)/365.25)^(i-1)*coeff_x(i)^(i-1);   % milli arcsec [mas]
%                 yp_t=((JD-JD2000)/365.25)^(i-1)*coeff_y(i)^(i-1);   % milli arcsec [mas]
%         end
% 
%         xp_t=deg2rad(xp_t/(1000*3600));
%         yp_t=deg2rad(yp_t/(1000*3600));
%         
%         if m==0 && n==2
%             Cnm=-0.48416948*10^-3+11.6*10^-12*(JD-JD2000)/365.25;
%         elseif m==0 && n==3
%             Cnm=0.9571612*10^-6+4.9*10^-12*(JD-JD2000)/365.25;
%         elseif m==0 && n==4 
%             Cnm=0.5399659*10^-6+4.7*10^-12*(JD-JD2000)/365.25;
%         elseif m==1 && n==2
%             row20=find(Coeff(:,1)==2 & Coeff(:,2)==0);
%             row22=find(Coeff(:,1)==2 & Coeff(:,2)==2);
%             
%             Cnm=sqrt(3)*xp_t*Coeff(row20,3)-xp_t*Coeff(row22,3)+yp_t*Coeff(row22,4);
%             Snm=-sqrt(3)*yp_t*Coeff(row20,3)-yp_t*Coeff(row22,3)-xp_t*Coeff(row22,4);
%         end
%         
%         % Solid Earth Tides................................................
%         
%         % Frequency invariant value
% %         if m==0 && n==2
% %             knm=Nominal_Love_numbers(n,m);
% %             Delta_CS=(knm/(2*n+1))*(GM_sun/GM)*(R_E/r_sun)^(n+1)*Leg_sun(m+1)*exp(-1j*long_sun)+(GM_moon/GM)*(R_E/r_moon)^(n+1)*Leg_moon(m+1)*exp(-1j*long_moon);
% %             Delta_CS_perm=(4.4228*10^-8)*(-0.31460)*knm;
% %             Delta_CS=Delta_CS-Delta_CS_perm;
% %             Delta_Cnm=real(Delta_CS);
% %             Delta_Snm=imag(Delta_CS);
% %         elseif n==4 && m<=2
% %             [~,knm_p]=Nominal_Love_numbers(2,m);
% %             Delta_CS=(knm_p/(2*n+1))*(GM_sun/GM)*(R_E/r_sun)^(n+1)*Leg_sun(m+1)*exp(-1j*long_moon)+(GM_moon/GM)*(R_E/r_moon)^(n+1)*Leg_moon(m+1)*exp(-1j*long_moon);
% %             Delta_Cnm=real(Delta_CS);
% %             Delta_Snm=imag(Delta_CS);
% %         elseif n==1
% %             Delta_Cnm=0;
% %             Delta_Snm=0;
% %         else
% %             if n<=3 && m<=3
% %                 knm=Nominal_Love_numbers(n,m);
% %                 Delta_CS=(knm/(2*n+1))*(GM_sun/GM)*(R_E/r_sun)^(n+1)*Leg_sun(m+1)*exp(-1j*long_moon)+(GM_moon/GM)*(R_E/r_moon)^(n+1)*Leg_moon(m+1)*exp(-1j*long_moon);
% %                 Delta_Cnm=real(Delta_CS);
% %                 Delta_Snm=imag(Delta_CS);
% %             else
% %                 Delta_Cnm=0;
% %                 Delta_Snm=0;
% %             end
% %         end
%         
%         %Frequency dependent correction
%         
% %         if n==2
% %             
% %             if m==0
% %                 Coeff_ST=K20_corr;
% %             else
% %                 Coeff_ST=K21_corr;
% %             end
% %             Delta_Freq=Solid_Tide_Frequency_correction(n,m,GMTrad,Coeff_ST,ttt);
% %             Delta_Freq_Cnm=real(Delta_Freq);
% %             Delta_Freq_Snm=imag(Delta_Freq);
% %         else
% %             Delta_Freq_Cnm=0;
% %             Delta_Freq_Snm=0;
% %         end
%         
%         % Ocean Tides......................................................
% 
% %         [Delta_Cnm_OT_Freq,Delta_Snm_OT_Freq]=Ocean_Tide_Frequency_correction(n,m,GMTrad,OT_data,Global_k,ttt);
%         
%         % Solid Earth Pole Tide............................................
%         
% %         if n==2 && m==1
% %             m1=rad2deg(xp-xp_t)*3600;
% %             m2=-rad2deg(yp-yp_t)*3600;
% %             Delta_EPT_Cnm=-1.333*10^-9*(m1+0.0115*m2);
% %             Delta_EPT_Snm=-1.333*10^-9*(m2-0.0115*m1);
% %         else
% %             Delta_EPT_Cnm=0;
% %             Delta_EPT_Snm=0;
% %         end
%         
%         % Ocean Pole Tide..................................................
%        
% %         m1=rad2deg(xp-xp_t)*3600;
% %         m2=-rad2deg(yp-yp_t)*3600;
% %         [Delta_OPT_Cnm,Delta_OPT_Snm]=Ocean_Pole_Tide(n,m,m1,m2,ABnm,hlk);
%         
%         %..................................................................
%         
% %         Cnm=Cnm+Delta_Cnm+Delta_Freq_Cnm+Delta_EPT_Cnm+Delta_OPT_Cnm+Delta_Cnm_OT_Freq;
% %         Snm=Snm+Delta_Snm+Delta_Freq_Snm+Delta_EPT_Snm+Delta_OPT_Snm+Delta_Snm_OT_Freq;
%         
%         Leg(m+1)=(-1)^(m)*Leg(m+1)*prod_nm;
%         dU_dr_p=dU_dr_p+Leg(m+1)*(Cnm*cos(m*long)+Snm*sin(m*long));
%         dU_dlong_p=dU_dlong_p+Leg(m+1)*m*(Snm*cos(m*long)-Cnm*sin(m*long));
%         
%         if m<n
%             prod_nm_plus1=sqrt((factorial(n-(m+1))*(2*n+1)*2)/factorial(n+m+1));
%             
%             Leg(m+2)=(-1)^(m+1)*Leg(m+2)*prod_nm_plus1;
%             
%             dU_dlat_p=dU_dlat_p+((Cnm*cos(m*long)+Snm*sin(m*long)))*...
%                 (Leg(m+2)-m*tan(lat)*Leg(m+1));
%         end
%         
%         
%     end
%     
%     dU_dr=dU_dr+((R_E/r_m)^n)*(n+1)*dU_dr_p;
%     dU_dlong=dU_dlong+(R_E/r_m)^n*dU_dlong_p;
%     dU_dlat=dU_dlat+(R_E/r_m)^n*dU_dlat_p;
% 
% end
% 
% dU_dr=-(GM/r_m^2)*dU_dr;
% dU_dlong=(GM/r_m)*dU_dlong;
% dU_dlat=(GM/r_m)*dU_dlat;

g=zeros(3,1);

g(1)=((1/r_m)*dU_dr-(r_K/(r_m^2*sqrt(r_I^2+r_J^2)))*dU_dlat)*r_I-(1/(r_I^2+r_J^2)*dU_dlong)*r_J;
g(2)=((1/r_m)*dU_dr-(r_K/(r_m^2*sqrt(r_I^2+r_J^2)))*dU_dlat)*r_J+(1/(r_I^2+r_J^2)*dU_dlong)*r_I;
g(3)=(1/r_m)*dU_dr*r_K+(sqrt(r_I^2+r_J^2)/r_m^2)*dU_dlat;
end

