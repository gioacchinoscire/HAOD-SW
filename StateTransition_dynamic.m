function [phi_dot]=StateTransition_dynamic(time,phi,order,EGM,EOP,DAT,A_m,JD0,state_interp)
% time=vettore dei tempi
% JD0

% Target dynamics..........................................................

curr_state=zeros(6,1);

[row_time,~]=find((state_interp(:,1)-time)>=-210 & (state_interp(:,1)-time)<=210);
delta=60;
while(isempty(row_time)==1)
    [row_time,~]=find((state_interp(:,1)-time)>=-(210+delta) & (state_interp(:,1)-time)<=(210+delta));
end
row_low=row_time-5;
row_high=row_time+5;

row_low=row_low(end);
row_high=row_high(end);

row_low=max(row_low,1);
row_high=max(min(row_high,length(state_interp(:,1))),11);


for i=2:7  
    curr_state(i-1)=lagrange_interpolation(state_interp(row_low:row_high,1)/86400,state_interp(row_low:row_high,i),time/86400);
end

A=zeros(6);
A(1:3,4:6)=eye(3);
v=zeros(6,1);
epsA=1E-5;

for i=1:6
    v(i)=epsA;
    gh_p=Earth_dynamic(time,curr_state+v,JD0,order,EGM,EOP,DAT,A_m);
    gh_m=Earth_dynamic(time,curr_state-v,JD0,order,EGM,EOP,DAT,A_m);
    A(4:6,i)=(gh_p(4:6)-gh_m(4:6))/(2*epsA);
    v(i)=0;
end

phi_mat=reshape(phi,[6,6]);
phi_dot=A*phi_mat;
phi_dot=reshape(phi_dot,[36,1]);


end
