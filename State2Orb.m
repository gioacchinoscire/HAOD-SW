function [Orb_el,varargout] =State2Orb(State)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function allows the Keplerian osculating parameters determination
% from the ECI state.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% State: Dynamic state time x 6(position,velocity)
% OUTPUT
% Orb_el: Orbital osculating elements (time x 6) (a,e,i,RAAN,Ap,M0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLE
% mu: Earth planetary costant [Km^3/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=398600.4415;


%% ECI REFERECE FRAME......................................................

c1=[1;0;0]; % X-axis
c2=[0;1;0]; % Y-axis
c3=[0;0;1]; % Z-axis

%% ORBTITAL PARAMETERS DETERMINATION.......................................

num_orb=length(State)/6;

Orb_el=zeros(size(State));     % Orbital elements 
Ecc_an=zeros(1,num_orb);

for i=1:num_orb
    
    r=State(6*(i-1)+1:6*(i-1)+3);           % ECI position vector
    v=State(6*(i-1)+4:6*i);                 % ECI velocity vector
    h=cross(r,v);                           % angular momentum vector
    e=(1/mu)*(cross(v,h))-(r./norm(r));     % eccentricity vector
    
    % Semi-major axis......................................................
    
    a=-mu/(((norm(v).^2/2)-(mu/norm(r)))*2);
    
    % Eccentricity.........................................................
    
    ec=norm(e);
    
    % Inclination..........................................................
    
    in=acos(dot(h,c3)./norm(h));
    
    % Node Right Ascension-Argument of perigee-Mean Anomaly................
    
    if ec<=1E-15 && in<=1E-15
        
        % Circular equatorial orbit........................................
        
        RAAN=0;
        Ap=0;
        EA0=atan2(r(2),r(1));
        M0=(EA0-ec*sin(EA0));
        
    elseif ec<=1E-15 && in>=1E-15
        
        % Circular Inclined orbit..........................................
        
        N=(cross(c3,h)./(norm(cross(c3,h))));
        M=cross(h,N)./norm(cross(h,N));
        RAAN=atan2(N(2),N(1));
        Ap=0;
        EA0=atan2(dot(r,M),dot(r,N));
        M0=(EA0-ec*sin(EA0));
 
    elseif in<=1E-15 && ec>=1E-15
        
        % Equatorial Non-circular orbit....................................
        
        RAAN=0;
        Ap=atan2(dot(e,c2),dot(e,c1));
        p=cross(h,e);
        Theta0=atan2(dot(r,p)./(norm(r)*norm(p)),...
            dot(r,e)./(norm(r)*ec));
        if ec<1
            SinEA0=((sqrt(1-norm(e).^2)*sin(Theta0))/(1+norm(e)...
                *(cos(Theta0))));
            CosEA0=((norm(e)+cos(Theta0))/(1+norm(e).*(cos(Theta0))));
            EA0=atan2(SinEA0,CosEA0);
            M0=(EA0-ec*SinEA0);
        else
            EA0=atanh((tan(Theta0/2)*sqrt((ec-1)/(ec+1))));
            M0=-(EA0-ec*sinh(EA0));
        end
    else
        
        N=(cross(c3,h)./(norm(cross(c3,h))));
        RAAN=atan2(N(2),N(1));
        p=cross(h,e);
        SIN_Ap=-(dot(N,p)./(norm(p)));
        COS_Ap=(dot(N,e)./(ec));
        Ap=atan2(SIN_Ap,COS_Ap);
        Theta0=atan2(dot((r),p)./(norm(r)*norm(p)),...
            dot(r,e)./(norm(r)*norm(e)));
        if ec<1
            SinEA0=((sqrt(1-norm(e).^2)*sin(Theta0))/(1+norm(e)...
                *(cos(Theta0))));
            CosEA0=((norm(e)+cos(Theta0))/(1+norm(e).*(cos(Theta0))));
            EA0=atan2(SinEA0,CosEA0);
            M0=(EA0-ec*SinEA0);
        else
            EA0=atanh((tan(Theta0/2)*sqrt((ec-1)/(ec+1))));
            M0=-(EA0-ec*sinh(EA0));
        end
        
    end
    
    % Assigment............................................................
    
    Orb_el(6*(i-1)+1:6*i)=[a,ec,in,RAAN,Ap,M0];
    
    Ecc_an(i)=EA0;
    
end
if nargout>1
    varargout{1}=Ecc_an;
end