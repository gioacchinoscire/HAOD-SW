
% scaricato da:

%   https://it.mathworks.com/matlabcentral/fileexchange/59567-astronomy-toolbox/content/SOFA/iauAb.m

%  - - - - - -
%   i a u A b
%  - - - - - -
%
%  Apply aberration to transform natural direction into proper
%  direction.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%    pnat(3)       natural direction to the source (unit vector)
%    v(3)          observer barycentric velocity in units of c
%    s             distance between the Sun and the observer (au)
%    bm1           sqrt(1-|v|^2): reciprocal of Lorenz factor
%
%  Returned:
%    ppr(3)        proper direction to source (unit vector)
%
%  Notes:
%
%  1) The algorithm is based on Expr. (7.40) in the Explanatory
%     Supplement (Urban & Seidelmann 2013), but with the following
%     changes:
%
%     o  Rigorous rather than approximate normalization is applied.
%
%     o  The gravitational potential term from Expr. (7) in
%        Klioner (2003) is added, taking into account only the Sun's
%        contribution.  This has a maximum effect of about
%        0.4 microarcsecond.
%
%  2) In almost all cases, the maximum accuracy will be limited by the
%     supplied velocity.  For example, if the SOFA iauEpv00 function is
%     used, errors of up to 5 microarcseconds could occur.
%
%  References:
%
%     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
%     the Astronomical Almanac, 3rd ed., University Science Books
%     (2013).
%
%     Klioner, Sergei A., "A practical relativistic model for micro-
%     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
%
%  Called:
%     iauPdp       scalar product of two p-vectors
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ppr = iauAb(pnat, v, s, bm1)

constants

% pdv = iauPdp(pnat, v);
pdv = dot(pnat, v);
w1 = 1.0 + pdv/(1.0 + bm1);
w2 = SRS/s;
r2 = 0.0;
p = zeros(3,1);
for i = 1:3
    w = pnat(i)*bm1 + w1*v(i) + w2*(v(i) - pdv*pnat(i));
    p(i) = w;
    r2 = r2 + w*w;
end
r = sqrt(r2);
ppr = zeros(3,1);
for i = 1:3
    ppr(i) = p(i)/r;
end
