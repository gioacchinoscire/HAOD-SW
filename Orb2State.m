function [State]=Orb2State(Orb_el)

mu=398600.4415;
dim=size(Orb_el);

a=Orb_el(1:6:end);
ec=Orb_el(2:6:end);
in=Orb_el(3:6:end);
RAAN=Orb_el(4:6:end);
Ap=Orb_el(5:6:end);
M0=Orb_el(6:6:end);

ind_c_e=(ec==0 & in==0);
ind_c_i=(ec==0)&(in~=0);
ind_e_e=(ec~=0)&(in==0);
ind_e=ec<=1;
ind_h=ec>1;
dim_v=size(ec);

M0(ind_c_e)=RAAN(ind_c_e)+Ap(ind_c_e)+M0(ind_c_e);
RAAN(ind_c_e)=0;
Ap(ind_c_e)=0;

M0(ind_c_i)=Ap(ind_c_i)+M0(ind_c_i);
Ap(ind_c_i)=0;

Ap(ind_e_e)=RAAN(ind_e_e)+Ap(ind_e_e);
RAAN(ind_e_e)=0;

p=abs(a).*abs((1-ec.^2));

EA=Kepler(ec,M0);

C_EA=zeros(dim_v);
S_EA=zeros(dim_v);
C_EA(ind_e)=cos(EA(ind_e));
S_EA(ind_e)=sin(EA(ind_e));
C_EA(ind_h)=cosh(EA(ind_h));
S_EA(ind_h)=sinh(EA(ind_h));
S_True_Anom=zeros(dim_v);
C_True_Anom=zeros(dim_v);
S_True_Anom(ind_e)=(sqrt(1-ec(ind_e).^2).*S_EA(ind_e))./(1-ec(ind_e).*C_EA(ind_e));
C_True_Anom(ind_e)=(C_EA(ind_e)-ec(ind_e))./(1-ec(ind_e).*C_EA(ind_e));
S_True_Anom(ind_h)=(sqrt(ec(ind_h).^2-1).*S_EA(ind_h))./(ec(ind_h).*C_EA(ind_h)-1);
C_True_Anom(ind_h)=(ec(ind_h)-C_EA(ind_h))./(ec(ind_h).*C_EA(ind_h)-1);

% Dynamic State in the Orbital Reference frame.............................

State_orb=zeros(dim);
State_orb(1:6:end)=p./(1+ec.*C_True_Anom).*C_True_Anom;
State_orb(2:6:end)=p./(1+ec.*C_True_Anom).*S_True_Anom;
State_orb(4:6:end)=-sqrt(mu./p).*S_True_Anom;
State_orb(5:6:end)=sqrt(mu./p).*(ec+C_True_Anom);
X_p=State_orb(1:6:end);
Y_p=State_orb(2:6:end);
X_v=State_orb(4:6:end);
Y_v=State_orb(5:6:end);

c_RA=cos(RAAN);
s_RA=sin(RAAN);
c_Ap=cos(Ap);
s_Ap=sin(Ap);
c_in=cos(in);
s_in=sin(in);

State=zeros(dim);
State(1:6:end)=(c_RA.*c_Ap-s_RA.*s_Ap.*c_in).*X_p-...
    (c_RA.*s_Ap+s_RA.*c_Ap.*c_in).*Y_p;
State(2:6:end)=(s_RA.*c_Ap+c_RA.*s_Ap.*c_in).*X_p-...
    (s_RA.*s_Ap-c_RA.*c_Ap.*c_in).*Y_p;
State(3:6:end)=s_Ap.*s_in.*X_p+c_Ap.*s_in.*Y_p;
State(4:6:end)=(c_RA.*c_Ap-s_RA.*s_Ap.*c_in).*X_v-...
    (c_RA.*s_Ap+s_RA.*c_Ap.*c_in).*Y_v;
State(5:6:end)=(s_RA.*c_Ap+c_RA.*s_Ap.*c_in).*X_v-...
    (s_RA.*s_Ap-c_RA.*c_Ap.*c_in).*Y_v;
State(6:6:end)=s_Ap.*s_in.*X_v+c_Ap.*s_in.*Y_v;

end

