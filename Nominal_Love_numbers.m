function [knm,knm_p]=Nominal_Love_numbers(n,m)

if n<=3 && m<=3
    Number_matrix=[2,0,0.30190,-0.00089;
                   2,1,0.29830-0.00144*1j,-0.00080;
                   2,2,0.30102-0.00130*1j,-0.00057;
                   3,0,0.093,0;
                   3,1,0.093,0;
                   3,2,0.093,0;
                   3,3,0.094,0];
    [row,~] = find(Number_matrix(:,1)==n & Number_matrix(:,2)==m);
    knm=Number_matrix(row,3);
    knm_p=Number_matrix(row,4);
else
   knm=0;
   knm_p=0;
end

end