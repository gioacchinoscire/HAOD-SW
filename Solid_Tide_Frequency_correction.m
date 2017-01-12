function [Corr]=Solid_Tide_Frequency_correction(n,m,theta,Coeff,ttt)
Corr=0;
% Delaunay Variables.......................................................

l=3600*134.96340251+1717915923.2178*ttt+31.8792*ttt^2+0.051635*ttt^3-0.00024470*ttt^4;
l_=3600*357.52910918+129596581.0481*ttt-0.5532*ttt^2+0.000136*ttt^3-0.00001149*ttt^4;
F=3600*93.27209062+1739527262.8478*ttt-12.7512*ttt^2-0.001037*ttt^3-0.00000417*ttt^4;
D=3600*297.85019547+1602961601.2090*ttt-6.3706*ttt^2+0.006593*ttt^3-0.00003169*ttt^4;
Om=3600*125.04455501-6962890.5431*ttt+7.4722*ttt^2+0.007702*ttt^3-0.00005939*ttt^4;
F_v=[l,l_,F,D,Om]./3600.*pi./180;

%..........................................................................
if n==2
    
    for i=1:length(Coeff(:,1))
        N=Coeff(i,9:13);
        if m==0
            ip=Coeff(i,15)*10^-12;
        else
            ip=Coeff(i,16)*10^-12;
            
        end
        op=Coeff(i,17)*10^-12;
        Corr=Corr+sqrt(ip^2+op^2)*exp(1j*(m*(theta+pi)-dot(N,F_v)));
        if m~=0
            Corr=(-(1j)^m)*Corr;
        end
    end
end
end