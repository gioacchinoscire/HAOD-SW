function [] =TDM_Parser(dir_path,sat_ref)

files_name=dir(dir_path);
RA_degrees=0;
DEC_degrees=0;
curr_year=0;
curr_month=0;
curr_day=0;
curr_hour=0;
curr_minutes=0;
curr_secs=0;
RA_found=0;
DEC_found=0;

for i=1:length(files_name)
    kvnfile=files_name(i).name;
    lng_name=length(kvnfile);
    if(lng_name>4)
        if(strcmp(kvnfile(end-3:end),'.kvn')==1)
            if(dir_path(end)=='\')
                Data=importdata([dir_path,kvnfile],'r');
            else
                Data=importdata([dir_path,'\',kvnfile],'r');
            end
            
            if(isstruct(Data)==1)
                Data=Data.textdata;
            end
            sat_file=Data{9}(17:end);
            if(strcmp(sat_file,sat_ref)==1)
                fileID=fopen(['TDM_Measures_',sat_ref,'.txt'],'a+');
                for j=1:length(Data)
                    if(length(Data{j})>7)
                        if(strcmp(Data{j}(1:7),'ANGLE_1')==1)
                            
                            if(Data{j}(30)=='.')
                                RA_degrees=str2double(Data{j}(33:end));
                                curr_secs=Data{j}(28:32);
                            else
                                RA_degrees=str2double(Data{j}(31:end));
                                curr_secs=Data{j}(28:29);
                            end
                            curr_year=Data{j}(11:14);
                            curr_month=Data{j}(16:17);
                            curr_day=Data{j}(19:20);
                            curr_hour=Data{j}(22:23);
                            curr_minutes=Data{j}(25:26);
                            
                            RA_found=1;

                        elseif(strcmp(Data{j}(1:7),'ANGLE_2')==1)
                            if(Data{j}(30)=='.')
                                DEC_degrees=str2double(Data{j}(33:end));
                            else
                                DEC_degrees=str2double(Data{j}(31:end));
                            end
                            DEC_found=1;
                        end
                        if(RA_found==1 && DEC_found==1)
                            fprintf(fileID,'%4s/%2s/%2s %2s:%2s:%05.3f\t%08.5f\t%08.5f\n',...
                                curr_year,curr_month,curr_day,curr_hour,curr_minutes,str2double(curr_secs),...
                                RA_degrees,DEC_degrees);
                            RA_found=0;
                            DEC_found=0;
                        end
                    end
                end
            end
        end
    end
end
end

