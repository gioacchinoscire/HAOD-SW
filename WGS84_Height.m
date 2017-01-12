function [H] =WGS84_Height(lat,long,order,EGM2008)
%-------------------------------WGS84_Height------------------------------%
%   INPUT:
%           lat:             latitude               [deg]            
%           long:            longitude              [deg]
%           order:           Series Max order       [70 max]
%   OUTPUT:
%           H:               Height from WGS84 Ellipsoid              [m]
%-------------------------------------------------------------------------%
global R_E f_E e_E GM g_eq

if(order>150)
    order=150;
end

% WGS84 Parameters---------------------------------------------------------

r_d=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*cosd(lat)^2)))*cosd(lat);
r_k=(R_E*sqrt(1-e_E^2)/(sqrt(1-e_E^2*cosd(lat)^2)))*sind(lat);
lat_geod=rad2deg(atan(r_k/(r_d*(1-e_E^2))));                    % Geodetic Latitude         [deg]
r=sqrt(r_d^2+r_k^2);                                            % Ellipsoid radius          [m]
 
b_E=f_E*R_E;                                                    % Minor semiaxis            [m]
g_pole=9.8321863685;                                            % Gravity at the pole       [m/s2]
kg=b_E*g_pole/(R_E*g_eq);
g_th=g_eq*((1+kg*sind(lat_geod)^2)/(sqrt(1-e_E^2*sind(lat_geod)^2)));

%--------------------------------------------------------------------------
H=0;
for n=2:order
    Leg=legendre(n,sind(lat));
    for m=0:n
        if m==0
            k=1;
        else
            k=2;
        end
        
        prod_nm=sqrt((factorial(n-m)*(2*n+1)*k)/factorial(n+m));
        [row,~]=find(EGM2008(:,1)==n & EGM2008(:,2)==m);
        
        Cnm=EGM2008(row,3);
        Snm=EGM2008(row,4);
        
       
        % Coefficient Correction...........................................
        
        if(mod(n,2)==0 && m==0)
            delta_C2n0=((-1)^(n/2))*3*e_E^(n)*(1-(n/2)-5^1.5*(n/2)*EGM2008(1,3)/e_E^2)/((n+1)*(n+3)*sqrt(2*n+1));
            Cnm=Cnm-delta_C2n0;
        end
        
        %..................................................................
        
        Leg(m+1)=(-1)^(m)*Leg(m+1)*prod_nm;
        H=H+(R_E/r)^n*Leg(m+1)*(Cnm*cosd(m*long)+Snm*sind(m*long));
    end
end
H=H*GM/(r*g_th)-0.41;

end

