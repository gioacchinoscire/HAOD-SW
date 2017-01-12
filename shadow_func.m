function [shadow_func]=shadow_func(Sun,r_sat)
% Geometric shadow function to consider penumbra regions
global R_E
alpha_umbra=0.26411888; % [deg]
alpha_penumbra=0.26900424; % [deg]

sat_hor=dot(r_sat,-Sun)/norm(Sun);
sat_ver=sqrt(norm(r_sat)^2-sat_hor^2);
c1=sat_hor+(R_E/1E3)/sind(alpha_penumbra);
c2=sat_hor-(R_E/1E3)/sind(alpha_umbra);
l1=c1*tand(alpha_penumbra);
l2=c2*tand(alpha_umbra);

if(sat_ver>l2 && sat_ver<l1)
    shadow_func=(sat_ver-l2)/(l1-l2);
elseif(sat_ver<l2)
    shadow_func=0;
elseif(sat_ver>l1)
    shadow_func=1;
end

    
end