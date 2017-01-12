function [v2] =Gibbs(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function use the Gibbs Method for the velocity estimation 
% with three position vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% r: 3x3 matrix of position vector[r(t1)|r(t2)|r(t3)]
% OUTPUT
% v2: velocity estimation at t2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=398600.4415; % Costante gravitazionale della terra [km3/s2]

r1=r(:,1);
r2=r(:,2);
r3=r(:,3);

Z12=cross(r1,r2);
Z23=cross(r2,r3);
Z31=cross(r3,r1);
N=norm(r1)*Z23+norm(r2)*Z31+norm(r3)*Z12;
D=Z12+Z23+Z31;
S=(norm(r2)-norm(r3))*r1+...
    (norm(r3)-norm(r1))*r2+...
    (norm(r1)-norm(r2))*r3;

B=cross(D,r2);
Lg=sqrt(mu/(norm(N)*norm(D)));


v2=(Lg/norm(r2))*B+Lg*S;


end

