function [state] =Planet_Ephemeris(JD,PLANET,REF,dut1,dat)
%-----------------------------Planet Ephemeris----------------------------%
% This function is used with the NASA Matlab Toolbox for the SPK format
% Ephemeris
%   INPUT:
%           JD:             Julian Date of propagation          [1xN]
%           PLANET:         Planet to identify wrt the Earth    [string]
%           SPKFILENAME:    Filename of the SPK-Format Ephemeris[string]
%   OUTPUT:
%           state:          Dynamic State of the Planet         [Nx6]
%   NOTE:   use the function  cspice_furnsh(SPKFILENAME) before This Func
%-------------------------------------------------------------------------%

size_JD=size(JD);
if(size_JD(2)==1 && size_JD(2)~=1)
    error('The JD Vector shall be a row vector 1xN')
end
timezone=0;

%Initialize values.........................................................
addpath( 'C:/mice/lib' )            % lib path mice
addpath( 'C:/mice/src/mice')        % src path mice
%..........................................................................
[year,mon,day,hr,min,sec] = invjday (JD);
[~, ~, ~, ~, ~, ~, ~, ~, ~, ttdb]=convtime(year,mon,day,hr,min,sec,timezone,dut1,dat);
[state, ~] = cspice_spkezr(PLANET,ttdb*36525.0*86400, 'J2000','LT',REF);
state=state';

end

