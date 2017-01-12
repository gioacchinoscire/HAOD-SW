function [Q] = cirs2gcrs(ttt,dX,dY)

F=zeros(1,14);
F(1:5)=Delaunay_variables(ttt);
F(6:14)=Planetary_nutation_arguments(ttt);

[X_coeff,Y_coeff,s_coeff]=CIP_coefficient();
X=-16617.0+2004191898.0*ttt-429782.9*ttt^2-198618.34*ttt^3+7.578*ttt^4+5.9285*ttt^5;
j=0;
for i=1:1600
    if(i==1307||i==1560||i==1596 ||i==1600)
        j=j+1;
    end
    N=X_coeff(i,4:end);
    arg=sum(F.*N);
    a_s=X_coeff(i,2);
    a_c=X_coeff(i,3);
    X=X+a_s*sin(arg)*ttt^j+a_c*cos(arg)*ttt^j;
end

Y=-6951.0-25896.0*ttt-22407274.7*ttt^2+1900.59*ttt^3+1112.526*ttt^4+0.1358*ttt^5;
j=0;
for i=1:1275
    if(i==963||i==1240||i==1270 ||i==1275)
        j=j+1;
    end
    N=Y_coeff(i,4:end);
    arg=sum(F.*N);
    a_s=Y_coeff(i,2);
    a_c=Y_coeff(i,3);
    Y=Y+a_s*sin(arg)*ttt^j+a_c*cos(arg)*ttt^j;
end

s=94.0+3808.65*ttt-122.68*ttt^2-72574.11*ttt^3+27.98*ttt^4+15.62*ttt^5;
j=0;
for i=1:66
    if(i==34||i==37||i==62 ||i==66)
        j=j+1;
    end
    N=s_coeff(i,4:end);
    arg=sum(F.*N);
    a_s=s_coeff(i,2);
    a_c=s_coeff(i,3);
    s=s+a_s*sin(arg)*ttt^j+a_c*cos(arg)*ttt^j;
end

X=deg2rad(X/1E6/3600);
Y=deg2rad(Y/1E6/3600);
s=deg2rad(s/1E6/3600);
s=s-(X*Y/2);

X=X+dX;
Y=Y+dY;

a=1/2+1/8*(X^2+Y^2);
R_s=[cos(s) sin(s) 0;-sin(s) cos(s) 0;0 0 1];
Q=[1-a*X^2 -a*X*Y X;-a*X*Y 1-a*Y^2 Y;-X -Y 1-a*(X^2+Y^2)]*R_s;


end

