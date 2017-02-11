function [dRA,dDEC]=Aberration_Astrometry(JD,RA,DEC,dut1,dat)
%---------------------------Aberration_Astrometry--------------------------%
% This function is used with the NASA Matlab Toolbox for the SPK format
% Ephemeris
%   INPUT:
%           RA:             Right Ascension         [deg][1xN]            
%           DEC:            Declination             [deg][1xN]
%           lat:            latitude                [deg]            
%           long:           longitude               [deg]
%           H_alt:          Geoig Height+Altitude from geoid
%           JD:             Julian Data             [deg][1xN]
%   OUTPUT:
%           dRA_arcsec:     Delta Right Ascension   [arcsec][1xN]
%           dDEC_arcsec:    Delta Declination       [arcsec][1xN]
%   NOTE:   use the function  cspice_furnsh(SPKFILENAME) before This Func
%-------------------------------------------------------------------------%

vers_stelle_mis=[cosd(RA).*cosd(DEC),cosd(DEC).*sind(RA),sind(DEC)];
AU=149597870;

RA_mis_corr =zeros(size(RA));
DEC_mis_corr=zeros(size(DEC));
dRA=zeros(size(RA));
dDEC=zeros(size(RA));
for icorr=1:length(JD)
    Sun_Helio=Planet_Ephemeris(JD(icorr),'EARTH','SUN',dut1(icorr),dat(icorr));
    Sun_Bary=Planet_Ephemeris(JD(icorr),'EARTH','SOLAR SYSTEM BARYCENTER',dut1(icorr),dat(icorr));
    V_c=Sun_Bary(4:6)/299792.458;
    pdv=dot(V_c,vers_stelle_mis(icorr,:));
    bm1=sqrt(1-norm(V_c)^2);
    w1 = 1.0 + pdv/(1.0 + bm1);
    SRS = 1.97412574336e-8*AU;
    w2 = SRS/norm(Sun_Helio(1:3));
    p=vers_stelle_mis(icorr,:)*bm1 + w1*V_c+w2*(V_c-pdv*vers_stelle_mis(icorr,:));
    ppr=p/norm(p);
    RA_mis_corr(icorr) = rad2deg(atan2(ppr(2),ppr(1)));
    DEC_mis_corr(icorr)= asind(ppr(3));
    if RA_mis_corr(icorr)<0
        RA_mis_corr(icorr)=RA_mis_corr(icorr)+360;
    end
    dRA(icorr)=RA_mis_corr(icorr)-RA(icorr);
    dDEC(icorr)=DEC_mis_corr(icorr)-DEC(icorr);
end

end
