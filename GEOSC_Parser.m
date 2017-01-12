function [RA_degrees,DEC_degrees] = GEOSC_Parser(filename)
fileID = fopen('GEOSC_Measures.txt','w');
geosc_data=importdata(filename,'r');
RA_degrees=zeros(size(geosc_data));
DEC_degrees=zeros(size(geosc_data));
Month=[31,29,31,30,31,30,31,31,30,31,30,31];

for i=1:length(geosc_data)
    
    year{i}=geosc_data{i}(17:18);
    dayYear{i}=geosc_data{i}(19:21);
    curr_month=1;
    while(str2double(dayYear{i})>sum(Month(1:curr_month)))
       curr_month=curr_month+1;
    end
    curr_day=str2double(dayYear{i})-sum(Month(1:curr_month-1));
    secDay{i}=geosc_data{i}(22:26);
    curr_hour=floor(str2double(secDay{i})/3600);
    curr_minutes=floor((str2double(secDay{i})-curr_hour*3600)/60);
    curr_secs=str2double(secDay{i})-curr_hour*3600 -curr_minutes*60;
    fracSec{i}=geosc_data{i}(27:32);
    curr_secs=str2double(secDay{i})-curr_hour*3600 -curr_minutes*60+str2double(fracSec{i})/1E6;
    RA_h{i}=geosc_data{i}(37:38);
    RA_m{i}=geosc_data{i}(39:40);
    RA_s{i}=geosc_data{i}(41:45);
    RA_degrees(i)=rad2deg(hms2rad(str2double(RA_h{i}),str2double(RA_m{i}),str2double(RA_s{i})/1E3));
    DEC_sign{i}=geosc_data{i}(46);           % Dec sign
    DEC_d{i}=geosc_data{i}(47:48);
    DEC_m{i}=geosc_data{i}(49:50);
    DEC_s{i}=geosc_data{i}(51:54);
    DEC_degrees(i)=rad2deg(dms2rad(str2double(DEC_d{i}),str2double(DEC_m{i}),str2double(DEC_s{i})/1E2));
    if(DEC_sign{i}=='-')
        DEC_degrees(i)=-DEC_degrees(i);
    end
    RA_std{i}=geosc_data{i}(58:61);
    DEC_std{i}=geosc_data{i}(62:65);
    

    fprintf(fileID,'%4s/%2s/%2s %2s:%2s:%021.18s %22.18s %22.18s\n',...
        ['20',num2str(year{i},'%2i')],num2str(curr_month,'%02i'),num2str(curr_day,'%02i'),...
        num2str(curr_hour,'%02i'),num2str(curr_minutes,'%02i'),num2str(curr_secs,'%21.18f'),...
        num2str(RA_degrees(i),'%22.18f'),num2str(DEC_degrees(i),'%22.18f'));


end


end

