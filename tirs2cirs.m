function [R] = tirs2cirs(jdut1)

Tu=jdut1-2451545;
era=2*pi*(0.7790572732640+1.00273781191135448*Tu);
cos_era=cos(-era);
sin_era=sin(-era);

R=[cos_era sin_era 0;-sin_era cos_era 0;0 0 1];

end

