function detector = simulateIFOs(N, seed)
%
%  simulate random locations and orientations of ifos on the
%  surface of a spherical earth.
%
%  The output is in the form of a cell-array of detector structures
%  detector{i}  with the fields
%      r: [3x1 double] %  position vector (in units of meters)
%                         in Cartesian coordinates
%      u: [3x1 double] %  unit vector along x-arm of detector
%                         in Cartesian coordinates
%      v: [3x1 double] %  unit vector along y-arm of detector
%                         in Cartesian coordinates
%      T:              %  arm length measured in light propagation
%                         time (in units of seconds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try, seed; catch, seed=1; end;
rand('state',seed);

% constants
c = 299792458; % speed of light (m/s)
RE = 6.3725e6; % earth radius (m)

for ii=1:N
  % put ifo at random locations on the surface of the earth, uniformly 
  % distributed in x = cos(theta) and phi
  x = -1 + 2*rand;
  sx = sqrt(1-x.^2);

  phi = 2*pi*rand;
  cphi = cos(phi);
  sphi = sin(phi);

  % detector position vector
  detector{ii}.r = RE*[sx*cphi; sx*sphi; x];

  % unit vectors tangent to sphere
  thetahat = [x*cphi; x*sphi; -sx];
  phihat = [-sphi; cphi; 0];

  % random rotation of ifo arms (wrt thetahat, phihat unit vectors)
  psi = 2*pi*rand;
  cpsi = cos(psi);
  spsi = sin(psi);

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

  % R4 = Ry(-theta)
  R4 = [  x,  0,  sx; ...
          0,  1,   0; ...
        -sx,  0,   x];

  % R5 = Rz(-phi)
  R5 = [ cphi, -sphi, 0; ...
         sphi,  cphi, 0; ...
            0,     0, 1];

  % combined rotation matrix 
  % R = Rz(-phi) * Ry(-theta) * Rz(psi) * Ry(theta) * Rz(phi) 
  R = R5 * R4 * R3 * R2 * R1;

  detector{ii}.u = R*thetahat;
  detector{ii}.v = R*phihat;

  % all ifos have 4km long arms
  detector{ii}.T = 4000/c;

end

return

