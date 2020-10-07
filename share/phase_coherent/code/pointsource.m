function pointsource(thetaI, phiI, res, lmax)
%
% calculate point source using grad and curl spherical harmonics
% out to lmax.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% angular resolution (arcsec) for healpix pixels
r = res*60*60; 
nSide = res2nSide(r);
nPix = nSide2nPix(nSide);

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);
    
% loop over pixels on sphere
vPix = zeros(1,nPix);
hpPix = zeros(1,nPix);
hcPix = zeros(1,nPix);
for ii = 1:1:nPix

  if mod(ii,50)==0
    fprintf('working on %d of %d\n', ii, nPix)
  end

  % extract theta, phi
  theta = tp{ii}(1);
  phi = tp{ii}(2);
  x = cos(theta);
    
  % initialize variables
  hp = 0;
  hc = 0;
  hpG = 0;
  hcG = 0;
  hpC = 0;
  hcC = 0;

  for l=2:1:lmax

    % normalization
    Nl = sqrt(2*factorial(l-2)/factorial(l+2));

    for m=-l:1:l

      % calculate Wlm, Xlm
      Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );      
      [Gp, Gm] = Glm(l, m, x);
      Wlm = 2*Nlm*Gp*exp(i*m*phi);
      Xlm = 2*sqrt(-1)*Nlm*Gm*exp(i*m*phi);

      % calculate aGlm, aClm for point source in direction (thetaI,phiI)
      [Gp, Gm] = Glm(l, m, cos(thetaI));
      Wlmp = conj(2*Nlm*Gp*exp(i*m*phiI));
      Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiI));
      aGlm = Nl*(Wlmp + Xlmp);
      aClm = Nl*(Wlmp - Xlmp);

      % calculate h+, hx
      hp = hp + 0.5*Nl*(aGlm*Wlm - aClm*Xlm);
      hc = hc + 0.5*Nl*(aGlm*Xlm + aClm*Wlm);

      % grad component
      hpG = hpG + 0.5*Nl*aGlm*Wlm;
      hcG = hcG + 0.5*Nl*aGlm*Xlm;

      % curl component
      hpC = hpC - 0.5*Nl*aClm*Wlm;
      hcC = hcC + 0.5*Nl*aClm*Xlm;

    end

  end

  hpPix(ii) = hp;
  hcPix(ii) = hc;
  hpGPix(ii) = hpG;
  hcGPix(ii) = hcG;
  hpCPix(ii) = hpC;
  hcCPix(ii) = hcC;
  vPix(ii) = abs(hp)^2+abs(hc)^2;

end

% make projections

figure
xyz = ang2vec(tp);
subplot(2,2,1)
m_mollweide(xyz, real(hpPix), 'colorbar', colorbar('southoutside'))
title('real(h+)')
subplot(2,2,2)
m_mollweide(xyz, imag(hpPix), 'colorbar', colorbar('southoutside'))
title('imag(h+)')
subplot(2,2,3)
m_mollweide(xyz, real(hcPix), 'colorbar', colorbar('southoutside'))
title('real(hx)')
subplot(2,2,4)
m_mollweide(xyz, imag(hcPix), 'colorbar', colorbar('southoutside'))
title('imag(hx)')
print('-djpeg', 'point.jpg')

figure
xyz = ang2vec(tp);
subplot(2,2,1)
m_mollweide(xyz, real(hpGPix), 'colorbar', colorbar('southoutside'))
title('real(h+) - grad')
subplot(2,2,2)
m_mollweide(xyz, imag(hpGPix), 'colorbar', colorbar('southoutside'))
title('imag(h+) - grad')
subplot(2,2,3)
m_mollweide(xyz, real(hcGPix), 'colorbar', colorbar('southoutside'))
title('real(hx) - grad')
subplot(2,2,4)
m_mollweide(xyz, imag(hcGPix), 'colorbar', colorbar('southoutside'))
title('imag(hx) - grad')
print('-djpeg', 'pointG.jpg')

figure
xyz = ang2vec(tp);
subplot(2,2,1)
m_mollweide(xyz, real(hpCPix), 'colorbar', colorbar('southoutside'))
title('real(h+) - curl')
subplot(2,2,2)
m_mollweide(xyz, imag(hpCPix), 'colorbar', colorbar('southoutside'))
title('imag(h+) - curl')
subplot(2,2,3)
m_mollweide(xyz, real(hcCPix), 'colorbar', colorbar('southoutside'))
title('real(hx) - curl')
subplot(2,2,4)
m_mollweide(xyz, imag(hcCPix), 'colorbar', colorbar('southoutside'))
title('imag(hx) - curl')
print('-djpeg', 'pointC.jpg')

%figure
%xyz = ang2vec(tp);
%m_mollweide(xyz, vPix, 'colorbar', colorbar('southoutside'))
%title('|h+|^2 + |hx|^2')
%print('-djpeg', 'point_power.jpg')

% save data for later use
filename = 'pointsource.mat';
save(filename, 'vPix', 'hpPix', 'hcPix')

return

