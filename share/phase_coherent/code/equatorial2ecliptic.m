function [elat, elon] = equatorial2ecliptic(ra, dec)
%
% convert from equatorial (ra,dec) to ecliptic latitude and
% longitude
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert (ra,dec) to standard (theta,phi) on 2-sphere
theta = pi/2-dec;
phi = (pi/12)*ra;

% location in (x,y,z)
x = sin(theta)*cos(phi);
y = sin(theta)*sin(phi);
z = cos(theta);

% rotate around x-axis (by tilt) to go from equatorial to
% ecliptic coordinates
const = physConstants('mks');
tilt = const.obliquity;
x_ecliptic =  x;
y_ecliptic =  y*cos(tilt) + z*sin(tilt);
z_ecliptic = -y*sin(tilt) + z*cos(tilt);

% convert to ecliptic latitude and longitude (in radians)
ecliptic_theta = acos(z_ecliptic);
elat = pi/2-ecliptic_theta;
elon = atan2(y_ecliptic, x_ecliptic);

return

