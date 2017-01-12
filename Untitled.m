Meas_data=importdata('GEOSC_Measures.txt',' ');

Meas=Meas_data.data;        % Optical measure matrix [RA (deg), DEC(deg)]

N=length(Meas);             % Measure Number

Time_str=Meas_data.textdata;% Measure Time vector YYYY/MM/DD h:m:s

Time_UTC=zeros(N,6);        % Measure Times UTC

Time_JD=zeros(N,1);         % Measure Times JD

for i=1:N
    
    data_str=Time_str{i,1};                          % Day
    day_str=Time_str{i,2};                           % Day instant(h:m:s)
    
    Time_UTC(i,1)=round(str2double(data_str(1:4)));  % year
    Time_UTC(i,2)=round(str2double(data_str(6:7)));  % month
    Time_UTC(i,3)=round(str2double(data_str(9:10))); % day
    Time_UTC(i,4)=round(str2double(day_str(1:2)));   % hours
    Time_UTC(i,5)=round(str2double(day_str(4:5)));   % minutes
    Time_UTC(i,6)=str2double(day_str(7:end));        % seconds
       
    Time_JD(i)=jday(Time_UTC(i,1),Time_UTC(i,2),Time_UTC(i,3),...
        Time_UTC(i,4),Time_UTC(i,5),Time_UTC(i,6))+3/86400;  % JD calculation   
end
