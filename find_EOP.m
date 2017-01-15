function [curr_EOP]=find_EOP(JD,EOP,DAT)
eopMJD=0;
i=1;
curr_EOP=zeros(1,7);
[year,month,day]=invjday(JD);
Day_MJD=jday(year,month,day,0,0,0)-2400000.5;
Sec_MJD=jday(year,month,day,0,0,0)-2400000.5;

while(Day_MJD~=eopMJD)
    sel_row=EOP{i};
    eopMJD=str2double(sel_row(13:19));
    if(Day_MJD==eopMJD)
        if(Sec_MJD==Day_MJD)
            curr_EOP(1)=deg2rad(str2double(sel_row(20:30))/3600); %xp
            curr_EOP(2)=deg2rad(str2double(sel_row(31:41))/3600); %yp
            curr_EOP(3)=str2double(sel_row(42:53)); %dut1
            curr_EOP(4)=str2double(sel_row(54:65)); %lod
            curr_EOP(5)=deg2rad(str2double(sel_row(66:76))/3600); %dX
            curr_EOP(6)=deg2rad(str2double(sel_row(77:87))/3600); %dY
        else
            if(i>=5 && i<=(length(EOP)-5))
                midMJD=zeros(10,1);
                mid_xp=zeros(10,1);
                mid_yp=zeros(10,1);
                mid_dut1=zeros(10,1);
                mid_lod=zeros(10,1);
                mid_dX=zeros(10,1);
                mid_dY=zeros(10,1);
                for j=i-4:i+5
                    mid_row=EOP{j};
                    midMJD(j)=str2double(sel_row(13:19));
                    mid_xp(j)=deg2rad(str2double(mid_row(20:30))/3600); %xp
                    mid_yp(j)=deg2rad(str2double(mid_row(31:41))/3600); %yp
                    mid_dut1(j)=str2double(mid_row(42:53)); %dut1
                    mid_lod(j)=str2double(mid_row(54:65)); %lod
                    mid_dX(j)=deg2rad(str2double(mid_row(66:76))/3600); %dX
                    mid_dY(j)=deg2rad(str2double(mid_row(77:87))/3600); %dY
                end
                curr_EOP(1)=lagrange_interpolation(midMJD,mid_xp,Sec_MJD); %xp
                curr_EOP(2)=lagrange_interpolation(midMJD,mid_yp,Sec_MJD); %yp
                curr_EOP(3)=lagrange_interpolation(midMJD,mid_dut1,Sec_MJD); %dut1
                curr_EOP(4)=lagrange_interpolation(midMJD,mid_lod,Sec_MJD); %lod
                curr_EOP(5)=lagrange_interpolation(midMJD,mid_dX,Sec_MJD); %dX
                curr_EOP(6)=lagrange_interpolation(midMJD,mid_dY,Sec_MJD); %dY
            else
                curr_EOP(1)=deg2rad(str2double(sel_row(20:30))/3600); %xp
                curr_EOP(2)=deg2rad(str2double(sel_row(31:41))/3600); %yp
                curr_EOP(3)=str2double(sel_row(42:53)); %dut1
                curr_EOP(4)=str2double(sel_row(54:65)); %lod
                curr_EOP(5)=deg2rad(str2double(sel_row(66:76))/3600); %dX
                curr_EOP(6)=deg2rad(str2double(sel_row(77:87))/3600); %dY
            end
        end
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