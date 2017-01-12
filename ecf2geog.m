function [geog]=ecf2geog(v_ecf)

geog(:,1)=rad2deg(atan2(v_ecf(:,2),v_ecf(:,1))); % longitude
normV=sqrt(dot(v_ecf',v_ecf'))';
geog(:,2)=rad2deg(asin(v_ecf(:,3)./normV));


end

