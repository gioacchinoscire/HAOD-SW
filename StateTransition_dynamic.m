function [dxdt]=StateTransition_dynamic(time,phi,order,Coeff,prop_time,Initial_State)
% global R_E GM
% 
% adim=[(R_E/1E3).*ones(3,1);sqrt((GM/1E9)/(R_E/1E3))*ones(3,1)];
adim=ones(6,1);
options=odeset('AbsTol',1E-3,'RelTol',1E-3);
Initial_State=Initial_State';

if(time~=prop_time(1)*86400)
    [~,Stato_sat]=ode113(@(t,x)Earth_dynamic(t,x,prop_time(1)/86400,order,Coeff),[prop_time(1),time/86400],Initial_State',options);
else
    Stato_sat=Initial_State;
end

A=zeros(6);
v=zeros(size(adim));
epsA=1E-6;
g=Earth_dynamic(time,Stato_sat(end,:)',prop_time(1),order,Coeff);
for i=1:6
    v(i)=epsA;
    gh_p=Earth_dynamic(time,Stato_sat(end,:)'+v,prop_time(1),order,Coeff);
%     gh_m=Earth_dynamic(time,Stato_sat(end,:)'-v,prop_time(1),order,Coeff);
%     gh_2p=Earth_dynamic(time,Stato_sat(end,:)'+2*v,prop_time(1),order,Coeff);
%     gh_2m=Earth_dynamic(time,Stato_sat(end,:)'-2*v,prop_time(1),order,Coeff);
%     A(:,i)=((gh_2m-8*gh_m+8*gh_p-gh_2p)/(12*epsA));
     A(:,i)=(gh_p-g)/(epsA);
    v(i)=0;
end

phi_mat=vec2mat(phi,6)';
phi_dot=A*phi_mat;

dxdt=phi_dot(:,1);
for i=2:6
    dxdt=[dxdt;phi_dot(:,i)];
end

end