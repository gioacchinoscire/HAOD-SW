function [J]=Ad_Jacobian(X0_ad,obs_pos,Cel_Coord,time,order,EGM,EOP,DAT,C_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the Jacobian of the residual vector with the
% central differences method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% X0_ad: non-dimensional dynamic state
% obs_pos: intertial observer position
% Cel_Coord: matrix of celestial coordinates [RA|DEC]
% time: Julian Date time at the measure instants
% C_lm,S_lm: Harmonics coefficient matrices
% OUTPUT
% J: Jacobian Nx6 (N:number of measures);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global R_E GM

adim=[(R_E/1E3).*ones(3,1);sqrt((GM/1E9)/(R_E/1E3))*ones(3,1)];


epsJ=1E-8;

v=zeros(6,1);

J=zeros(length(Cel_Coord(:,1)),7);

% f0=Residual_vec(X0_ad.*adim,time,obs_pos,order,Coeff,Cel_Coord,C_r);


for i=1:6
    v(i)=1;
    
    Xperturb_pos=X0_ad+epsJ*v;

    fpos=Residual_vec(Xperturb_pos.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r);
    
    Xperturb_neg=X0_ad-epsJ*v;
    
    fneg=Residual_vec(Xperturb_neg.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r);
    
%     J(:,i)=(Residual_vec(Xperturb_pos.*adim,time,obs_pos,order,...
%         Coeff,Cel_Coord,C_r)-f0)/(epsJ);
    
    J(:,i)=(fpos-fneg)/(2*epsJ);
    v(i)=0;
    
end

C_r_pos=C_r+epsJ;
fpos=Residual_vec(X0_ad.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r_pos);
C_r_neg=C_r-epsJ;
fneg=Residual_vec(X0_ad.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r_neg);

J(:,7)=(fpos-fneg)/(2*epsJ);
% J(:,7)=(Residual_vec(X0_ad.*adim,time,obs_pos,order,...
%         Coeff,Cel_Coord,C_r+epsJ)-f0)/(epsJ);

end

