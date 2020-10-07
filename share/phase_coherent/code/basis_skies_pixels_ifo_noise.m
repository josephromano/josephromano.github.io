function basis_skies_pixels_ifo_noise(f, res, N, Nt, type)
%
% calculate basis sky maps using pixel-based approach for a network 
% of ground-based interferometers in the small antenna limit allowing
% for simulated noise
%
% f   - frequency (Hz)
% res - resolution (e.g., r=7 -> nPix=3072; r=10 -> nPix = 768)
% N   - number of ifos  
% Nt  - number of times
% type - 'rotation', 'orbit' or 'both'
%
% NOTE: if N=1  use L1
%       if N=3  use H1, L1, V1
%       if N=6  use real ifo network (H1, L1, V1, K1, I1, A1) 
%       if N=12 use real ifo network + antipodal detectors
%       else simulate ifos on surface of earth
%
% example: basis_skies_ifo_noise(100, 10, 6, 200, 'both')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
DEBUG=0;
const = physConstants('mks');

% set flags
type
switch type
  case 'rotation'
    rotationonly = 1;    
    orbitonly = 0;    
  case 'orbit'
    rotationonly = 0;    
    orbitonly = 1;
  case 'both'
    rotationonly = 0;    
    orbitonly = 0;
  otherwise
    error('unrecognized type');
end

% interferometer network
if N==1

  %det{1} = getdetectorNew('default');
  det{1} = getdetectorNew('L1');

elseif N==3
  det{1} = getdetectorNew('H1');
  det{2} = getdetectorNew('L1');
  det{3} = getdetectorNew('V1');

elseif N==6
  det{1} = getdetectorNew('H1');
  det{2} = getdetectorNew('L1');
  det{3} = getdetectorNew('V1');
  det{4} = getdetectorNew('K1');
  det{5} = getdetectorNew('I1');
  det{6} = getdetectorNew('A1');

elseif N==12
  det{1} = getdetectorNew('H1');
  det{2} = getdetectorNew('L1');
  det{3} = getdetectorNew('V1');
  det{4} = getdetectorNew('K1');
  det{5} = getdetectorNew('I1');
  det{6} = getdetectorNew('A1');
  det{7} = getdetectorNew('aH1');
  det{8} = getdetectorNew('aL1');
  det{9} = getdetectorNew('aV1');
  det{10} = getdetectorNew('aK1');
  det{11} = getdetectorNew('aI1');
  det{12} = getdetectorNew('aA1');

else 
  % simulate ifo locations, orientations
  det = simulateIFOs(N);
end

% polarisation angle
psi = 0;

% angular frequency of Earth's daily rotation (rad/s)
if orbitonly
  wE = 0;
else
  wE = 2*pi/(const.sidDay);
end

% discrete times
t = (Nt/60)*[0:1/Nt:1-1/Nt]*const.sidDay;

% angular resolution (arcsec) for healpix pixels
r = res*60*60; 
nSide = res2nSide(r);
nPix = nSide2nPix(nSide); % number of modes

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);

% loop over interferometers
R = zeros(N*Nt, 2*nPix);
for kk=1:N

  % get detector information
  detector = det{kk};

  % extract detector parameters
  r = detector.r;
  u = detector.u;
  v = detector.v;

  % calculate (theta,phi) values for ifos at t=0
  thetaI(kk) = acos(r(3)/sqrt(sum(r.^2)));
  phiI(kk) = atan2(r(2),r(1));

  if DEBUG
    % display lat, long
    lat = (pi/2-thetaI(kk))*180/pi;
    lon = phiI(kk)*180/pi;
    fprintf('detector %d: lat = %f deg, lon = %f deg\n', kk, lat, lon);
  end

  % rotate detector unit vectors u,v keeping source fixed (equatorial coords)
  % (rotation by -wE*t around z-axis)
  u_t = zeros(3,Nt);
  u_t(1,:) =  cos(-wE*t)*u(1) + sin(-wE*t)*u(2);
  u_t(2,:) = -sin(-wE*t)*u(1) + cos(-wE*t)*u(2);
  u_t(3,:) =  u(3);
  v_t = zeros(3,Nt);
  v_t(1,:) =  cos(-wE*t)*v(1) + sin(-wE*t)*v(2);
  v_t(2,:) = -sin(-wE*t)*v(1) + cos(-wE*t)*v(2);
  v_t(3,:) =  v(3);

  % rotate position vector of detector from center of  earth
  % (rotation by -wE*t around z-axis)
  rdE = zeros(3,Nt);
  rdE(1,:) =  cos(-wE*t)*r(1) + sin(-wE*t)*r(2);
  rdE(2,:) = -sin(-wE*t)*r(1) + cos(-wE*t)*r(2);
  rdE(3,:) =  r(3);
  
  % calculate position vector of center of Earth relative to SSB
  % in equatorial coordinates due to earth's yearly orbital motion
  if rotationonly==1
      % fixed position vector of Earth 
      rES(1,:) = const.au*ones(size(t));
      rES(2,:) = zeros(size(t));
      rES(3,:) = zeros(size(t));
    else
      % orbiting Earth
      rES = earthOrbit(t);
  end

  % calculate position vector of detector relative to SSB
  if rotationonly==1
    %r_t = rdE; % rotation of earth at SSB
    r_t = rdE + rES; % rotation of earth at 1 AU from SSB
  else
    r_t = rdE + rES;
  end

  % loop over pixels on sphere
  Fp = zeros(Nt, nPix);
  Fc = zeros(Nt, nPix);
  for ii = 1:1:nPix

    % extract theta, phi (direction to source)
    theta = tp{ii}(1);
    phi = tp{ii}(2);

    % first convert to antipodal point (khat = - direction to source)
    theta = pi - theta;
    phi = pi + phi;

    % calculate antenna patterns
    antenna = FpFc(theta, phi, psi, u_t, v_t, r_t, f);
    Fp(:,ii) = antenna.Fp;
    Fc(:,ii) = antenna.Fc;

  end

  % calculate response matrix for pixel approach
  % R matrix is just +,x response functions
  R((kk-1)*Nt+1:kk*Nt,:) = [Fp, Fc];

end

% save all data to .mat file
filename = ['pix_' num2str(nPix) '_ifo_noise_' type '_' num2str(N) '_times_' num2str(Nt) '.mat'];
save(filename)

return
