function [State]=Gauss_OD(kmax,Cel_Coord,time_JD,obs_pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function perform the Gauss orbit determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% kmax: max iteration value
% Cel_Coord: matrix of celestial coordinates [RA|DEC]
% time_JD: Julian Date time at the measure instants
% obs_pos: observer inertial position
% OUTPUT
% State: Dynamic state at time_JD(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=398600.4415;

% Direction of observation in ECI Reference frame..........................

RA_gauss=(Cel_Coord(1:3,1));
DEC_gauss=(Cel_Coord(1:3,2));

% Matrix of the three direction of measurement.............................

L=[cosd(RA_gauss).*cosd(DEC_gauss),...
    sind(RA_gauss).*cosd(DEC_gauss),sind(DEC_gauss)]';

% time_JD parameters.......................................................

t1=time_JD(1)*86400;
t2=time_JD(2)*86400;
t3=time_JD(3)*86400;

tau_1=time_JD(1)*86400-time_JD(2)*86400;
tau_3=time_JD(3)*86400-time_JD(2)*86400;

a1=tau_3/(tau_3-tau_1);
b1=a1*((tau_3-tau_1)^2-tau_3^2)/6;
a3=-tau_1/(tau_3-tau_1);
b3=-a3*((tau_3-tau_1)^2-tau_1^2)/6;

% observer positions at the measurements time_JD...........................

r_site=obs_pos;

M=L\r_site;
d1=M(2,1)*a1-M(2,2)+M(2,3)*a3;
d2=M(2,1)*b1+M(2,3)*b3;

C=dot(L(:,2),r_site(:,2));

% range estimation.........................................................

r2_range=roots([1 0 -(d1^2+2*C*d1+norm(r_site(:,2))^2) 0 0 ...
    -2*mu*(C*d2+d1*d2) 0 0 -mu^2*d2^2]);

for i=1:8
    if isequal(imag(r2_range(i)),0)==0 || real(r2_range(i))<0;
        r2_range(i)=0;
    end
end

r2_range=sort(nonzeros(r2_range),'descend');

if length(r2_range)>1
    r2_range=r2_range(2);
end

u=mu/r2_range^3;

c1=a1+b1*u;
c2=-1;
c3=a3+b3*u;

c=[c1;c2;c3];

rho=(-M*c)./c;

% Positions Estimation.....................................................

r_est=[rho(1)*L(:,1)+r_site(:,1),rho(2)*L(:,2)+r_site(:,2),...
    rho(3)*L(:,3)+r_site(:,3)];

% Velocity estimation(Gibbs Method)........................................

r2=r_est(:,2);

d_theta1=acosd(dot(r_est(:,1)./norm(r_est(:,1)),r_est(:,2)./norm(r_est(:,2))));
d_theta3=acosd(dot(r_est(:,3)./norm(r_est(:,3)),r_est(:,2)./norm(r_est(:,2))));

if(mean([d_theta1,d_theta3])<3)
    v2=HGibbs(r_est,time_JD);
else
    v2=Gibbs(r_est);
end

State=[r2;v2];

err=100;

k=1;

% Iterative correction.....................................................

while err>1e-8 && k<kmax
    
    [Orb_par,E_A2]=State2Orb(State');
    
    a2=Orb_par(1);
    ec2=Orb_par(2);
    M2=Orb_par(6);
    
    
    % Eccentric anomalies..................................................
    
    if ec2<1
        n=sqrt(mu/a2^3);
        M1=M2+n*(t1-t2);
        M3=M2+n*(t3-t2);
        E_A1=fzero(@(x) x-ec2.*sin(x)-M1,0);
        E_A3=fzero(@(x) x-ec2.*sin(x)-M3,0);
    else
        n=sqrt(mu/-a2^3);
        M1=M2+n*(t1-t2);
        M3=M2+n*(t3-t2);
        
        E_A1=fzero(@(x) -x+ec2.*sinh(x)-M1,0);
        E_A3=fzero(@(x) -x+ec2.*sinh(x)-M3,0);
    end
    
    % f and g series applicaition..........................................
    
    f1=1-(1-cos(E_A1-E_A2))/(1-ec2*cos(E_A2));
    f3=1-(1-cos(E_A3-E_A2))/(1-ec2*cos(E_A2));
    g1=(t1-t2)-(1/n)*((E_A1-E_A2)-sin(E_A1-E_A2));
    g3=(t3-t2)-(1/n)*((E_A3-E_A2)-sin(E_A3-E_A2));
    
    c1=g3/(f1*g3-f3*g1);
    c3=-g1/(f1*g3-f3*g1);
    
    c=[c1;c2;c3];
    
    rho_n=(-M*c)./c;
    
    err=norm(rho_n-rho,Inf);
    
    rho=rho_n;
    
    % Positions Estimation.................................................
    
    r_est=[rho(1)*L(:,1)+r_site(:,1),rho(2)*L(:,2)+r_site(:,2),...
        rho(3)*L(:,3)+r_site(:,3)];

    % Velocity estimation(Gibbs Method)....................................
    
    d_theta1=acosd(dot(r_est(:,1)./norm(r_est(:,1)),r_est(:,2)./norm(r_est(:,2))));
    d_theta3=acosd(dot(r_est(:,3)./norm(r_est(:,3)),r_est(:,2)./norm(r_est(:,2))));
    
    if(mean([d_theta1,d_theta3])<3)
        v2=HGibbs(r_est,time_JD);
    else
        v2=Gibbs(r_est);
    end

    %  updating.......................................................
    
    k=k+1;
    State=[r2;v2];
    
end
end

