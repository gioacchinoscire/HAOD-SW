function [H]=observ_grad(x0,time,obs_pos,order,EGM,EOP,DAT,Meas,A_m)
global R_E GM

adim=[(R_E/1E3).*ones(3,1);sqrt((GM/1E9)/(R_E/1E3))*ones(3,1)];

[~,~,Cel_Coord]=Residual_vec(x0,time,obs_pos,order,EGM,EOP,DAT,Meas,A_m);
v=zeros(6,1);
H=zeros(numel(Cel_Coord),6);
epsG=1E-8;
meas_vec=zeros(numel(Cel_Coord),1);
meas_vec(1:2:end)=Cel_Coord(:,1);
meas_vec(2:2:end)=Cel_Coord(:,2);
meas_vec_pert=zeros(numel(Cel_Coord),1);
for i=1:6
    v(i)=epsG;
    v=v.*adim;
    [~,~,Cel_Coord_pert]=Residual_vec(x0+v,time,obs_pos,order,EGM,EOP,DAT,Meas,A_m);
    meas_vec_pert(1:2:end)=Cel_Coord_pert(:,1);
    meas_vec_pert(2:2:end)=Cel_Coord_pert(:,2);
    H(:,i)=(meas_vec_pert-meas_vec)/epsG;
    v=zeros(6,1);
end

end