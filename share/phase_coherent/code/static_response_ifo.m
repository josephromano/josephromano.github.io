function [RG, RC] = static_response_ifo(l, m, u, v)
%
% calculate response of ifo at origin
%
% Inputs:
% l, m  - labels mode (l>=2, |m|<=l)
% u, v  - unit vectors along the arms of the interferometer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% extract angles
thetau = acos(u(3));
phiu = atan2(u(2), u(1));
thetav = acos(v(3));
phiv = atan2(v(2), v(1));

% calculate RG, RC (curl response = 0)
RC = 0;

if l==2 & abs(m)<=2
  Yu = sphericalharmonic(2, m, thetau, phiu);
  Yv = sphericalharmonic(2, m, thetav, phiv);
  RG = (4*pi/5)*sqrt(1/3)*(Yu - Yv);
else
  RG = 0;
end

return

