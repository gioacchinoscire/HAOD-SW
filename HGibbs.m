function [v2] =HGibbs(r,JD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function use the Herrick-Gibbs Method for the velocity estimation 
% with three position vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r: 3x3 matrix of position vector[r(t1)|r(t2)|r(t3)]
% JD: 3x1 vector of julian dates [t1;t2;t3]
% OUTPUT
% v2: velocity estimation at t2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=398600.4415; % Costante gravitazionale della terra [km3/s2]

r1=r(:,1);
r2=r(:,2);
r3=r(:,3);

t31=(JD(3)*86400-JD(1)*86400);
t32=(JD(3)*86400-JD(2)*86400);
t21=(JD(2)*86400-JD(1)*86400);

v2=-t32*(1/(t21*t31)+mu/(12*norm(r1)^3))*r1+...
    (t32-t21)*(1/(t21*t32)+mu/(12*norm(r2)^3))*r2+...
    t21*(1/(t32*t31)+mu/(12*norm(r3)^3))*r3;
end

