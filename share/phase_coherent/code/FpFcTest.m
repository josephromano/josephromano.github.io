% 
% script for testing FpFc.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = physConstants('mks');

% zero frequency
f = 0;

% polarisation angle
psi = 0;

% get detector information
site = 'H1';
detector = getdetectorNew(site);

% extract detector parameters
r = detector.r;
u = detector.u;
v = detector.v;

% angular frequency of Earth's daily rotation (rad/s)
wE = 2*pi/const.sidDay;

% discrete times
t = [0:.02:1]*const.sidYr/52;
Nt = length(t);

% rotate detector unit vectors u,v keeping source fixed (equatorial coords)
% (rotation by -wE*t around z-axis)
u_t = zeros(3,Nt);
u_t(1,:) =  cos(-wE*t)*u(1) + sin(-wE*t)*u(2);
u_t(2,:) = -sin(-wE*t)*u(1) + cos(-wE*t)*u(2);
u_t(3,:) =  u(3);
v_t(1,:) =  cos(-wE*t)*v(1) + sin(-wE*t)*v(2);
v_t(2,:) = -sin(-wE*t)*v(1) + cos(-wE*t)*v(2);
v_t(3,:) =  v(3);
r_t(1,:) =  cos(-wE*t)*r(1) + sin(-wE*t)*r(2);
r_t(2,:) = -sin(-wE*t)*r(1) + cos(-wE*t)*r(2);
r_t(3,:) =  r(3);

% angular resolution (arcsec) for healpix pixels
res = 5*60*60; 
nSide = res2nSide(res);
nPix = nSide2nPix(nSide);

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);

% calculate antenna patterns as functions of (theta,phi)
fp = zeros(Nt, nPix);
fc = zeros(Nt, nPix);
for ii=1:1:nPix

  if mod(ii,50)==0
    fprintf('working on %d of %d\n', ii, nPix)
  end

  % extract theta, phi (direction to source)
  theta = tp{ii}(1);
  phi = tp{ii}(2);

  % first convert to antipodal point (khat = - direction to source)
  theta = pi - theta;
  phi = pi + phi;

  % calculate antenna patterns
  antenna = FpFc(theta, phi, psi, u_t, v_t, r_t, f);
  fp(:,ii) = antenna.Fp;
  fc(:,ii) = antenna.Fc;

end

% make plots
xyz = ang2vec(tp);
for jj=1:Nt

  figure(1)
  m_mollweide(xyz, real(fp(jj,:)))
  title('real(h+)')
  pause(.01)

end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
xyz = ang2vec(tp);
for jj=1:Nt

  figure
  subplot(2,2,1)
  m_mollweide(xyz, real(fp(jj,:)))
  title('real(h+)')
  subplot(2,2,2)
  m_mollweide(xyz, imag(fp(jj,:)))
  title('imag(h+)')
  subplot(2,2,3)
  m_mollweide(xyz, real(fc(jj,:)))
  title('real(hx)')
  subplot(2,2,4)
  m_mollweide(xyz, imag(fc(jj,:)))
  title('imag(hx)')
  pause(.01)

end
