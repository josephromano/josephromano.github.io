function makeprojection(a_Plm, res, filename, thetaI, phiI)
%
% make a mollweide projection map of a stochastic background 
% described by a_Plm
%
% Inputs
%    a_Plm    - mode coefficients packed as an lm array with grad 
%               and curl blocks
%    res      - angular resolution in degrees (e.g., 4 degrees)
%    filename - filename to save pixel data 
%               (defaults to 'map_mollweide.mat')
%
% NOTE: a_Plm is a structure with fields:
%
%   a_Plm.data = [a_Glm; a_Clm];
%   a_Plm.lvec = lvec;
%   a_Plm.mvec = mvec;
%   a_Plm.lmax = lmax;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angular resolution (arcsec) for healpix pixels
r = res*60*60; 
nSide = res2nSide(r);
nPix = nSide2nPix(nSide);

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);

% extract a_Glm, a_Clm, lmax
a_Glm = a_Plm.data(1:end/2);
a_Clm = a_Plm.data(end/2+1:end);
lmax = a_Plm.lmax;

% loop over pixels on sphere
vPix = zeros(1,nPix);
hpPix = zeros(1,nPix);
hcPix = zeros(1,nPix);
for ii = 1:1:nPix

  if mod(ii,50)==0
    fprintf('working on %d of %d\n', ii, nPix)
  end

  % extract theta, phi (direction to source)
  theta = tp{ii}(1);
  phi = tp{ii}(2);

  % first convert to antipodal point (khat = - direction to source)
  theta = pi - theta;
  phi = pi + phi;

  x = cos(theta);
  
  % calculate Nl
  Nl = calNl(lmax);
  
  % calculate Wlm, Xlm
  [Wlm, Xlm] = calWlmXlm(theta, phi, lmax);

  % calculate h+, hx
  hp = 0.5*sum(Nl.data.*a_Glm.*Wlm.data - Nl.data.*a_Clm.*Xlm.data);
  hc = 0.5*sum(Nl.data.*a_Glm.*Xlm.data + Nl.data.*a_Clm.*Wlm.data);
  
  hpPix(1,ii) = hp;
  hcPix(1,ii) = hc;
  vPix(1,ii) = abs(hp)^2+abs(hc)^2;

end

% make projections
xyz = ang2vec(tp);

figure
subplot(2,2,1)
m_mollweide(xyz, real(hpPix), 'colorbar', colorbar('southoutside'));
if exist('thetaI')
  hold on
  show_pulsars(thetaI, phiI)
end
title('real(h+)')

subplot(2,2,2)
m_mollweide(xyz, imag(hpPix), 'colorbar', colorbar('southoutside'));
if exist('thetaI')
  hold on
  show_pulsars(thetaI, phiI)
end
title('imag(h+)')

subplot(2,2,3)
m_mollweide(xyz, real(hcPix), 'colorbar', colorbar('southoutside'));
if exist('thetaI')
  hold on
  show_pulsars(thetaI, phiI)
end
title('real(hx)')

subplot(2,2,4)
m_mollweide(xyz, imag(hcPix), 'colorbar', colorbar('southoutside'));
if exist('thetaI')
  hold on
  show_pulsars(thetaI, phiI)
end
title('imag(hx)')

%figure
%m_mollweide(xyz, vPix, 'colorbar', colorbar('southoutside'));
%if exist('thetaI')
%  hold on
%  show_pulsars(thetaI, phiI)
%end
%title('|h+|^2+|hx|^2')

% save pixel data for later use
try
  filename;
catch  
  filename = 'map_mollweide.mat';
end

save(filename, 'vPix', 'hpPix', 'hcPix', 'a_Plm')

return

