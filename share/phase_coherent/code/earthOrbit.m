function rE = earthOrbit(t)
%
% EARTHORBIT - calculates the position vector of the center of earth 
% relative to SSB in equatorial coordinates taking into account
% earth's yearly orbital motion. (assumes circular orbit of r=1 a.u. 
% and sun at SSB.)
%
% rE = earthOrbit(t, type)
% 
% t  - array of times (sec)
%
% rE - position vector of earth [3, numel(t)] relative to SSB in 
%      equatorial coordinates
%
% $Id:$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = physConstants('mks');
wES = 2*pi/const.sidYr; % angular freq (rad/sec)

% first calculate components in ecliptic coordinates
x = const.au*cos(wES*t);
y = const.au*sin(wES*t);
z = zeros(size(t));

% obliquity of ecliptic in radians
epsilon = const.obliquity;

% rotate around x-axis (by -epsilon) to go from ecliptic to 
% equatorial coordinates
rE(1,:) =  x;
rE(2,:) =  y*cos(epsilon) - z*sin(epsilon);
rE(3,:) =  y*sin(epsilon) + z*cos(epsilon); 

%%
return
