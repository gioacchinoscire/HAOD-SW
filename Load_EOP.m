function [EOP_vector] =Load_EOP(Data,Web)
%---------------------------------Load EOP--------------------------------%
% This function Load the EOP directly trough IERS FTP
%   INPUT:
%           Data:           Days in YYYY/MM/DD Format                 [Nx3]
%   OUTPUT:
%           EOP_matrix:     Earth Orientation Parameters (xp,yp,dut1,los,dPsi,dEps,dat)    [Nx7]
%   NOTE:   This Function Need an Active Internet connection
%-------------------------------------------------------------------------%

EOP_vector=zeros(length(Data(:,1)),7);          % EOP Vector Initialization
if(Web==1)
    f=ftp('ftp.iers.org');                          % Connect to the IERS FTP
    cd(f,'/products/eop/rapid/standard/');          % Change directory toward Standard Daily bullettin
    Bullettins=dir(f);                              % Identify all the bullettins in the folder
    if(isempty(f)==0)
       disp('connected to  ftp://ftp.iers.org')
    end
end
% Check which bullettin shall be read......................................
if(isempty(find(Data(:,1)<1992, 1)))
    bull_choice='finals.data';
else
    bull_choice='finals.all';
end
exsist=0;
if(Web==1)
    for i=1:length(Bullettins)
        if(strcmp(Bullettins(i).name,bull_choice))
            exsist=1;
        end
    end
    if(exsist==1)
        mget(f,bull_choice);
    end
end
%..........................................................................

%Read Bullettin and save EOP...............................................

Bullettin_read=importdata(bull_choice,'r');

for i=1:length(Bullettin_read)
    curr_row=Bullettin_read{i};
    curr_MJD=str2double(curr_row(8:15));
    for j=1:length(Data(:,1))
        JD= jday(Data(j,1),Data(j,2),Data(j,3),00,00,00);
        MJD=JD-2400000.5;
        if(MJD==curr_MJD)
            EOP_vector(j,1)=str2double(curr_row(19:27))/3600*pi/180;            % x polar motion    [arcsec->rad]
            EOP_vector(j,2)=str2double(curr_row(38:46))/3600*pi/180;            % y polar motion    [arcsec->rad]
            EOP_vector(j,3)=str2double(curr_row(59:68));                        % UT1-UTC           [sec]
            EOP_vector(j,4)=str2double(curr_row(80:86))/1000;                   % LOD               [millisec->sec]
            EOP_vector(j,5)=str2double(curr_row(100:106))/3600*pi/180/1000;     % dPsi              [milliarcsec->rad]
            EOP_vector(j,6)=str2double(curr_row(119:125))/3600*pi/180/1000;     % dEps              [milliarcsec->rad]
        end
    end
    
end
%..........................................................................
% Find TAI-UTC leap second.................................................
if(Web==1)
    cd(f,'/products/eop/rapid/bulletina/');
    Bull_directories=dir(f);
    index=1;
    i=1;
    while(i<=length(Bull_directories))
        cd(f,['/products/eop/rapid/bulletina/' Bull_directories(i).name '/']);
        Bull_files=dir(f);
        if(isempty(Bull_files)==0)
            j=1;
            while(j<=length(Bull_files))
                if(Bull_files(j).isdir==0)
                    mget(f,Bull_files(j).name);
                    Text=importdata(Bull_files(j).name,'\n');
                    k=1;
                    while(k<=length(Text))
                        [st_str,end_str]=regexp(Text{k},'\s+TAI-UTC(\(*)\w*(\)*)\s=\s\d+.\d+\s\d+ seconds');    % Regular Expression to find the TAI-UTC Field
                        if(isempty(st_str)==0)
                            data=Text{k}(st_str:end_str);
                            [st_num,end_num]=regexp(data,'\d+.\d+\s\d+');
                            data=Text{k}(st_num:end_num);
                            [sp_st,sp_end]=regexp(data,'\s+');
                            Delta_TAI=str2double([data(1:sp_st-1),data(sp_end+1:end)]);
                            [st_str,end_str]=regexp(Text{k-1},'Beginning \d+ \w+ \d+');    % Regular Expression to find the TAI-UTC Field
                            if(isempty(st_str)==0)
                                data=Text{k-1}(st_str:end_str);
                                [st_str,end_str]=regexp(data,'\d+ \w+ \d+');
                                file_date=datevec(data(st_str:end_str));
                                JD_file=jday(file_date(1),file_date(2),file_date(3),00,00,00);
                                TAI_Data(index,:)=[Delta_TAI,JD_file];
                                index=index+1;
                                k=length(Text);
                            end
                        end
                        k=k+1;
                    end
                    delete(Bull_files(j).name);
                end
                j=j+1;
            end
        end
        i=i+1;
    end
    save TAI_Data
else
    load TAI_Data
end
Old_JD=0;
for h=1:length(Data(:,1))
    Input_JD= jday(Data(h,1),Data(h,2),Data(h,3),00,00,00);   
    [~,index_dat]=min(abs(TAI_Data(:,2)-Input_JD));
    
    if(Old_JD~=Input_JD)
        if(TAI_Data(index_dat(1),2)>Input_JD)
            eop_h=TAI_Data(index_dat(1),1)-1;
        else
            eop_h=TAI_Data(index_dat(1),1);
        end
        EOP_vector(h,7)= eop_h;                         % dat   [seconds]
        Old_JD=Input_JD;
    else
        EOP_vector(h,7)=EOP_vector(h-1,7);              % dat   [seconds]
    end
end
close(f);
end

