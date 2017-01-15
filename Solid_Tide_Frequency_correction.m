function [Corr]=Solid_Tide_Frequency_correction(n,m,theta,ttt)
Corr=0;
Coeff=ST_coeff(m);
% Delaunay Variables.......................................................

F_v=Delaunay_variables(ttt);

%..........................................................................
if n==2
    N=Coeff(:,9:13);
    if m==0
        ip=Coeff(:,15)*10^-12;
    else
        ip=Coeff(:,16)*10^-12;
    end
    op=Coeff(:,17)*10^-12;
    Corr=sqrt(ip.^2+op.^2).*exp(1j*(m*(theta+pi)-sum((N.*((ones(length(Coeff(:,1)),1)*F_v)))')'));
    if m~=0
        Corr=(-(1j)^m)*Corr;
    end
    Corr=sum(Corr);
end

end