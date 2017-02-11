function [fx]=lagrange_interpolation(a,fa,x)
fx=0;
[row,~]=find(a==x);
if(isempty(row)==1)
    for i=1:length(a)
        part_prod=1;
        for j=1:length(a)
            if(i~=j)
                part_prod=part_prod*(x-a(j))/(a(i)-a(j));
            end
        end
        fx=fx+fa(i)*part_prod;
    end
else
    fx=fa(row);
end

end