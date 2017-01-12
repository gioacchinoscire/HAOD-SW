function [rasc, decl, rsun] = sun1 (jdate)

% solar ephemeris

% input
    
%  jdate = julian day

% output

%  rasc = right ascension of the sun (radians)
%         (0 <= rasc <= 2 pi)
%  decl = declination of the sun (radians)
%         (-pi/2 <= decl <= pi/2)
%  rsun = eci position vector of the sun (kilometers)

% note

%  coordinates are inertial, geocentric,
%  equatorial and true-of-date

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atr = pi / 648000;

rsun = zeros(3, 1);

% time arguments

djd = jdate - 2451545;

t = (djd / 36525) + 1;

% fundamental trig arguments (radians)

gs = r2r(0.993126 + 0.0027377785 * djd);
lm = r2r(0.606434 + 0.03660110129 * djd);
ls = r2r(0.779072 + 0.00273790931 * djd);
g2 = r2r(0.140023 + 0.00445036173 * djd);
g4 = r2r(0.053856 + 0.00145561327 * djd);
g5 = r2r(0.056531 + 0.00023080893 * djd);
rm = r2r(0.347343 - 0.00014709391 * djd);

% geocentric, ecliptic longitude of the sun (radians)

plon = 6910 * sin(gs) + 72 * sin(2 * gs) - 17 * t * sin(gs);
plon = plon - 7 * cos(gs - g5) + 6 * sin(lm - ls) ... 
       + 5 * sin(4 * gs - 8 * g4 + 3 * g5);
plon = plon - 5 * cos(2 * (gs - g2)) - 4 * (sin(gs - g2) ... 
       - cos(4 * gs - 8 * g4 + 3 * g5));
plon = plon + 3 * (sin(2 * (gs - g2)) - sin(g5) - sin(2 * (gs - g5)));
plon = ls + atr * (plon - 17 * sin(rm));

% geocentric distance of the sun (kilometers)

rsm = 149597870.691 * (1.00014 - 0.01675 * cos(gs) - 0.00014 * cos(2 * gs));

% obliquity of the ecliptic (radians)

obliq = atr * (84428 - 47 * t + 9 * cos(rm));

% geocentric, equatorial right ascension and declination (radians)

a = sin(plon) * cos(obliq);
b = cos(plon);
   
rasc = atan3(a, b);
decl = asin(sin(obliq) * sin(plon));

% geocentric position vector of the sun (kilometers)

rsun(1) = rsm * cos(rasc) * cos(decl);
rsun(2) = rsm * sin(rasc) * cos(decl);
rsun(3) = rsm * sin(decl);

