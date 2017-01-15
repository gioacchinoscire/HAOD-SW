function [R] = tirs2cirs(jdut1)

era=ERA(jdut1);
cos_era=cos(-era);
sin_era=sin(-era);

R=[cos_era sin_era 0;-sin_era cos_era 0;0 0 1];

end

