function [fx]=lagrange_interpolation(a,fa,x)
fx=0;
for i=1:length(a)
    part_sum=0;
    for j=1:length(a)
        if(i~=j)
            part_sum=part_sum+(x-a(j))/(a(i)-a(j));
        end
    end
    
    fx=fx+fa(i)*part_sum;
end
end