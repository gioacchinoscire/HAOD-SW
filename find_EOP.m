function [curr_EOP]=find_EOP(JD,EOP,DAT)
eopMJD=0;
curr_EOP=zeros(1,7);

[year, month, day,hr, minute,sec]=invjday(JD);
Day_MJD=jday(year,month,day,0,0,0)-2400000.5;
Sec_MJD=JD-2400000.5;

% Starting position for search in EOP file.................................

Month_day=[31,28,31,30,31,30,31,31,30,31,30,31];    % days in moth
leap_year=1964:4:year;                              % leap years
isleap=find(leap_year==year,1);
if(isempty(find(leap_year==year,1))==0)
    Month_day(2)=29;                                % add day for leap year
end

i=365*(year-1962)+(length(leap_year)-1);
if(month~=1)
    i=i+Month_day(month-1);
end

% IERS Earth Orientation Parameters (IAU2000/2006A)........................

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
                for j=-4:+5
                    mid_row=EOP{j+i};
                    midMJD(j+5)=str2double(mid_row(13:19));
                    mid_xp(j+5)=deg2rad(str2double(mid_row(20:30))/3600); %xp
                    mid_yp(j+5)=deg2rad(str2double(mid_row(31:41))/3600); %yp
                    mid_dut1(j+5)=str2double(mid_row(42:53));%dut1
                    mid_lod(j+5)=str2double(mid_row(54:65)); %lod
                    mid_dX(j+5)=deg2rad(str2double(mid_row(66:76))/3600); %dX
                    mid_dY(j+5)=deg2rad(str2double(mid_row(77:87))/3600); %dY
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
if(sum(curr_EOP)==0)
   error('Missing Earth Orientation Parameters, Please check the date'); 
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

% timezone=0;
% dut1=curr_EOP(3);
% dat=curr_EOP(7);
% [~, ~, jdut1, ~, ~, ~, ttt]= convtime (year, month, day, hr, minute,sec, timezone, dut1, dat);   % Time conversion

% % Libration Effects........................................................
% 
% [Deltaxp_lib,Deltayp_lib,DeltaUT1_lib,DeltaLOD_lib]=libration_effects(ttt,jdut1);
% 
% Effects of Tidal Deformation on Earth Rotation...........................

% [DeltaUT1,DeltaLOD,DeltaOmega]=Tidal_Effect_On_Earth_Rotation(ttt);         % Verified by IERS Example

% % Diurnal and SemiDiurnal variations due to ocean tides....................
% 
% [Deltaxp_subD,Deltayp_subD,DeltaUT1_subD,DeltaLOD_subD]=SubDiurnal_Earth_Rotation_Variation(ttt,jdut1);
% 
% curr_EOP(1)=curr_EOP(1)+Deltaxp_lib+Deltaxp_subD;
% curr_EOP(2)=curr_EOP(2)+Deltayp_lib+Deltayp_subD;
% curr_EOP(3)=curr_EOP(3)+DeltaUT1+DeltaUT1_lib+DeltaUT1_subD;
% curr_EOP(4)=curr_EOP(4)+DeltaLOD+DeltaLOD_lib+DeltaLOD_subD;

end