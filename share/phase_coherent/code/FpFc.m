function antenna = FpFc(theta, phi, psi, u, v, rd, f)
%
% FPFC - calculates antenna pattern functions and their derivatives 
% with respect to source position x=cos(theta), phi, and polarisation 
% angle psi in the long-wavelength limit.
%
% antenna = FpFc(theta, phi, psi, u, v)
%
% Note: Fp = d^ab ep_ab e^{-i2pifk.r/c}
%       Fc = d^ab ec_ab e^{-i2pifk.r/c}
%       d^ab = (1/2)(u^a u^b - v^a v^b) 
% 
% theta, phi - source location spherical polar angles (radians)
%              (theta, phi define unit vector n pting toward source; 
%               wave propagates in opposite direction k = -n) 
% psi        - polarisation angle
% u, v       - unit vectors pointing along interferometer arms
%              (can be arrays u(1,:), u(2,:), u(3,:), etc. 
%              corresponding to detector orientation at different t)
% rd         - position vector of detector
%              (can be arrays r(1,:), r(2,:), r(3,:), etc. 
%              corresponding to detector location at different t)
% f          - frequency
%
% antenna    - structure containing antenna pattern functions and
%              their derivatives
%              field    content
%              Fp       plus antenna pattern functions
%              Fc       cross antenna pattern functions
%
% $Id:$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = physConstants('mks');

% shorthand notation for variables
cphi = cos(phi);
sphi = sin(phi);
cpsi = cos(psi);
spsi = sin(psi);
x = cos(theta);
sx = sqrt(1-x^2);

% wave propagation vector
k = [sin(theta)*cphi, sin(theta)*sphi, x];

% calculate rotation matrices and their derivatives
% R1 = Rz(phi)
R1 = [ cphi, sphi, 0; ...
      -sphi, cphi, 0; ...
          0,    0, 1];

% R2 = Ry(theta) = Ry(x=cos(theta))
R2 = [  x,  0, -sx; ...
        0,  1,   0; ...
       sx,  0,   x];

% R3 = Rz(psi)
R3 = [ cpsi, spsi, 0; ...
      -spsi, cpsi, 0; ...
          0,    0, 1];

% combined rotation matrix and its derivatives
% R = Rz(psi) * Ry(theta) * Rz(phi) 
R = R3 * R2 * R1;

% rotate unit vectors along detector arms into source frame
%Ru = R * u;
%Rv = R * v;
Ru(1,:) = R(1,1)*u(1,:)+R(1,2)*u(2,:)+R(1,3)*u(3,:);
Ru(2,:) = R(2,1)*u(1,:)+R(2,2)*u(2,:)+R(2,3)*u(3,:);
Ru(3,:) = R(3,1)*u(1,:)+R(3,2)*u(2,:)+R(3,3)*u(3,:);

Rv(1,:) = R(1,1)*v(1,:)+R(1,2)*v(2,:)+R(1,3)*v(3,:);
Rv(2,:) = R(2,1)*v(1,:)+R(2,2)*v(2,:)+R(2,3)*v(3,:);
Rv(3,:) = R(3,1)*v(1,:)+R(3,2)*v(2,:)+R(3,3)*v(3,:);

% contract polarisation tensors with detector tensor to calculate
% Fp, Fc, and their derivatives

% calculate phase term
kdotr = k(1)*rd(1,:) + k(2)*rd(2,:) + k(3)*rd(3,:);
phase = exp(-sqrt(-1)*2*pi*f*kdotr/const.c);

% calculate Fp = (1/2)(u^a u^b - v^a v^b) ep_ab e^{-i2pifk.r/c}
uuep = Ru(1,:).^2 - Ru(2,:).^2;
vvep = Rv(1,:).^2 - Rv(2,:).^2;
Fp = (1/2)*(uuep - vvep).*phase;

% calculate Fc = (1/2)(u^a u^b - v^a v^b) ec_ab e^{-i2pifk.r/c}
uuec = 2*Ru(1,:).*Ru(2,:);
vvec = 2*Rv(1,:).*Rv(2,:);
Fc = (1/2)*(uuec - vvec).*phase;

% fill structure
antenna.Fp = Fp;
antenna.Fc = Fc;

return
