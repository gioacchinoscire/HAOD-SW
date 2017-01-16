function [P]=Associate_Legendre(n_max,x)

%Legendre associate function calculation...................................

P=zeros(n_max+1,n_max+1);
for m=0:n_max-1
    row_m=m+1;
    if(m~=0)
        P(row_m,row_m)=(2*m-1)*sqrt(1-x^2)*P(row_m-1,row_m-1);
    else
        P(row_m,row_m)=1;
    end
    P(row_m+1,row_m)=(2*m+1)*x*P(row_m,row_m);
    col_m=m+1;
    for n=m+2:n_max
        row_n=n+1;
        P(row_n,col_m)=(1/(n-m))*((2*n-1)*x*P(row_n-1,col_m)-(n+m-1)*P(row_n-2,col_m));
    end
end
P(n_max+1,n_max+1)=(2*n_max-1)*sqrt(1-x^2)*P(n_max,n_max);