% Measures='Measure_object_3.txt'
Measures={'Measure_object_1.txt','Measure_object_2.txt','Measure_object_3.txt'};

for j=1:3
Meas_data=importdata(Measures{j},' ');

% Matrice contenente le misure ottiche [RA (deg), DEC(deg)]

Meas=Meas_data.data;
save Meas
% Vettore di celle contenente gli istanti temporali di misura espressi nel
% formatio g/m/y h:m:s

Time_str=Meas_data.textdata;

% Inizializzazione dei vettori per la scansione degli istanti temporali di
% misura

N=length(Meas); %Numero di misure

Time_UTC=zeros(N,6); % Vettore degli istanti temporali in UTC
Time_JD=zeros(N,1);  % Vettore degli istanti temporali in julian date

for i=1:N
    
    data_str=Time_str{i,1};      % Giorno
    day_str=Time_str{i,2};       % Istante del giorno
    
    Time_UTC(i,1)=round(str2double(data_str(1:4)));  % year
    Time_UTC(i,2)=round(str2double(data_str(6:7)));  % month
    Time_UTC(i,3)=round(str2double(data_str(9:10))); % day
    Time_UTC(i,4)=round(str2double(day_str(1:2)));   % hours
    Time_UTC(i,5)=round(str2double(day_str(4:5)));   % minutes
    Time_UTC(i,6)=str2double(day_str(7:end));        % seconds
    
    % Calcolo delle date juliane per ogni istante di misura
    
    Time_JD(i)=jday(Time_UTC(i,1),Time_UTC(i,2),Time_UTC(i,3),...
        Time_UTC(i,4),Time_UTC(i,5),Time_UTC(i,6));
end

plot(Meas(:,1),Meas(:,2),'.');hold on
end