function [r_itrs,v_itrs,a_itrs] = gcrs2itrs(r_gcrs,v_gcrs,a_gcrs,ttt,jdut1,lod,xp,yp,dX,dY)

w_E=7.292115*10^-5*(1-lod/86400)*[0;0;1];                     % Earth Angular Velocity    [rad/s]
[W] = itrs2tirs(ttt,xp,yp);
[R] = tirs2cirs(jdut1);
[Q] = cirs2gcrs(ttt,dX,dY);

r_itrs=(Q*R*W)'*r_gcrs;
r_tirs=(Q*R)'*r_gcrs;
v_tirs=(Q*R)'*v_gcrs;

v_itrs=W'*((Q*R)'*v_gcrs-cross(w_E,r_tirs));
a_itrs=W'*((Q*R)'*a_gcrs-cross(w_E,cross(w_E,r_tirs))-2*cross(w_E,v_tirs));

end
