function [curr_EOP]=find_EOP(JD,EOP,DAT)
eopMJD=0;
i=1;
curr_EOP=zeros(1,7);
[year,month,day]=invjday(JD);
Day_MJD=jday(year,month,day,0,0,0)-2400000.5;

while(Day_MJD~=eopMJD)
    sel_row=EOP{i};
    eopMJD=str2double(sel_row(13:19));
    if(Day_MJD==eopMJD)
        curr_EOP(1)=deg2rad(str2double(sel_row(20:30))/3600); %xp
        curr_EOP(2)=deg2rad(str2double(sel_row(31:41))/3600); %yp
        curr_EOP(3)=str2double(sel_row(42:53)); %dut1
        curr_EOP(4)=str2double(sel_row(54:65)); %lod
        curr_EOP(5)=deg2rad(str2double(sel_row(66:76))/3600); %dX
        curr_EOP(6)=deg2rad(str2double(sel_row(77:87))/3600); %dY
        break
    end
    i=i+1;
end
i=1;
while(i<length(DAT))
    sel_row=DAT{i};
    datMJD=str2double(sel_row(18:28))-2400000.5;
    if(datMJD<Day_MJD)
        curr_EOP(7)=str2double(sel_row(37:48));
    end
    i=i+1;
end

end