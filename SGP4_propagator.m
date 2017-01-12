function [Dynamic_State] =SGP4_propagator(lg1,lg2,time_JD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function perform the dynamic integration using the SGP4 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% lg1: first tle line
% lg2 second tle line
% time_JD: vector of propagation julian dates
% OUTPUT
% Dynamic_State: Nx6 matrix that contains the dynamic state th the N time
% values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[satrec]=twoline2rv(84,lg1,lg2,'c');

[year_ril,mon_ril,day_ril,hr_ril,min_ril,sec_ril]=invjday(satrec.jdsatepoch );% qui calcola anno mese giorno ora del tle

jd_ril=jday(year_ril, mon_ril, day_ril, hr_ril, min_ril, sec_ril) ;

jd_prop=time_JD(1);

tempo_da_rilascio_dei_tle= jd_prop-jd_ril;

if tempo_da_rilascio_dei_tle>5
    
    warning('TechCorp:InconsistentDataType','The used TLE has released more than 5 days ago');
    
end

integ_time=(time_JD-jd_ril)*1440;

Dynamic_State=zeros(length(integ_time),6);

for i=1:length(integ_time)
    
    [satrec, r, v] = sgp4(satrec,integ_time(i));
    Dynamic_State(i,1)=r(1);
    Dynamic_State(i,2)=r(2);
    Dynamic_State(i,3)=r(3);
    Dynamic_State(i,4)=v(1);
    Dynamic_State(i,5)=v(2);
    Dynamic_State(i,6)=v(3);
    
end

end

