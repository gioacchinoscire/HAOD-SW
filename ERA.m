function [era_rad]=ERA(jdut1)

Tu=jdut1-2451545;
era_rad=2*pi*(0.7790572732640+1.00273781191135448*Tu);

end