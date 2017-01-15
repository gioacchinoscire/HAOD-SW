function [knm,knm_p]=Nominal_Love_numbers(n,m)
knm=zeros(length(n),1);
knm_p=zeros(length(n),1);
for i=1:length(knm)
    if n(i)<=3 && m(i)<=3
        Number_matrix=[2,0,0.30190,-0.00089;
            2,1,0.29830-0.00144*1j,-0.00080;
            2,2,0.30102-0.00130*1j,-0.00057;
            3,0,0.093,0;
            3,1,0.093,0;
            3,2,0.093,0;
            3,3,0.094,0];
        [row,~] = find(Number_matrix(:,1)==n(i) & Number_matrix(:,2)==m(i));
        knm(i)=Number_matrix(row,3);
        knm_p(i)=Number_matrix(row,4);
    else
        knm(i)=0;
        knm_p(i)=0;
    end
end
end