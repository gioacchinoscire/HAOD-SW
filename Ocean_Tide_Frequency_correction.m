function [Delta_Cnm,Delta_Snm]=Ocean_Tide_Frequency_correction(n,m,theta,CoeffOT,Coeff,ttt)
Corr=0;
% Delaunay Variables.......................................................

l=3600*134.96340251+1717915923.2178*ttt+31.8792*ttt^2+0.051635*ttt^3-0.00024470*ttt^4;
l_=3600*357.52910918+129596581.0481*ttt-0.5532*ttt^2+0.000136*ttt^3-0.00001149*ttt^4;
F=3600*93.27209062+1739527262.8478*ttt-12.7512*ttt^2-0.001037*ttt^3-0.00000417*ttt^4;
D=3600*297.85019547+1602961601.2090*ttt-6.3706*ttt^2+0.006593*ttt^3-0.00003169*ttt^4;
Om=3600*125.04455501-6962890.5431*ttt+7.4722*ttt^2+0.007702*ttt^3-0.00005939*ttt^4;
F_v=[l,l_,F,D,Om]./3600.*pi./180;

%..........................................................................

[row,~]=find(CoeffOT(:,2)==n & CoeffOT(:,3)==m);
Doodson=CoeffOT(row,2);
if isempty(Doodson)==0
    for k=1:length(Doodson)
        Cnm_p=CoeffOT(row(k),4)*10^-12;
        Snm_p=CoeffOT(row(k),5)*10^-12;
        Cnm_m=CoeffOT(row(k),6)*10^-12;
        Snm_m=CoeffOT(row(k),7)*10^-12;
        supp_vec=(Coeff(:,2)-Doodson(k));
        [rowfr,~]=find(supp_vec<0.01);
        for i=1:length(rowfr)
            N=Coeff(rowfr(i),9:13)*10^-12;
            Theta_f=m*(theta+pi)-dot(N,F_v);
            Corr=Corr+(Cnm_p-1j*Snm_p)*exp(1j*Theta_f)+(Cnm_m+1j*Snm_m)*exp(-1j*Theta_f);
        end
    end
end
Delta_Cnm=real(Corr);
Delta_Snm=imag(Corr);
end