function [W] = itrs2tirs(ttt,xp,yp)

s_=-((47/1E6)/3600)*ttt;
sin_s=sind(-s_);
cos_s=cosd(-s_);
sin_xp=sin(xp);
cos_xp=cos(xp);
sin_yp=sin(yp);
cos_yp=cos(yp);

R1_yp=[1 0 0;0 cos_yp sin_yp;0 -sin_yp cos_yp];
R2_xp=[cos_xp 0 -sin_xp;0 1 0;sin_xp 0 cos_xp];
R3_s=[cos_s sin_s 0;-sin_s cos_s 0;0 0 1];

W=R3_s*R2_xp*R1_yp;

end

