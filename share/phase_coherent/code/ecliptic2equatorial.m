function [ra, dec] = ecliptic2equatorial(elat, elon)
%
% convert from ecliptic latitude and longitude to equatorial (ra,dec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert (elat,elon) to standard (theta,phi) on 2-sphere
theta = pi/2-elat;
phi = elon;

% location in (x,y,z)
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);
z = cos(theta);

% rotate around x-axis (by -tilt) to go from ecliptic to
% equatorial coordinates
const = physConstants('mks');
tilt = const.obliquity;
x_equatorial =  x;
y_equatorial =  y*cos(-tilt) + z*sin(-tilt);
z_equatorial = -y*sin(-tilt) + z*cos(-tilt);

% convert to equatorial ra and dec
equatorial_theta = acos(z_equatorial);
equatorial_phi = atan2(y_equatorial, x_equatorial);
ra = equatorial_phi*(12/pi);
dec = pi/2-equatorial_theta;

return

