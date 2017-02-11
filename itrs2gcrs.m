function [r_gcrs,v_gcrs,a_gcrs] = itrs2gcrs(r_itrs,v_itrs,a_itrs,ttt,jdut1,lod,xp,yp,dX,dY)

w_E=7.292115*10^-5*(1-lod/86400)*[0;0;1];                     % Earth Angular Velocity    [rad/s]

[W] = itrs2tirs(ttt,xp,yp);
[R] = tirs2cirs(jdut1);
[Q] = cirs2gcrs(ttt,dX,dY);

r_gcrs=Q*R*W*r_itrs;
r_tirs=W*r_itrs;

v_tirs=W*v_itrs;
v_gcrs=Q*R*(W*v_itrs+cross(w_E,r_tirs));
a_gcrs=Q*R*(W*a_itrs+cross(w_E,cross(w_E,r_tirs))+2*cross(w_E,v_tirs));

end

