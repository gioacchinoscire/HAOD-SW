function [GMST]=Greenswich_MST(ttt,jdut1)

era=ERA(jdut1);
GMST =era+deg2rad((0.01450600 + 4612.15653400*ttt + 1.391581700*ttt^2-0.0000004400*ttt^3-0.00002995600*ttt^4-0.000000036800*ttt^5)/3600);
end