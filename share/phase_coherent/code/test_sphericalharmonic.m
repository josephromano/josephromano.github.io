function test_sphericalharmonic(l, m, res)
%
% test script for spherical harmonic code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% angular resolution (arcsec) for healpix pixels
r = res*60*60; 
nSide = res2nSide(r);
nPix = nSide2nPix(nSide)

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);

% loop over pixels on sphere
vPix = zeros(1,nPix);
for ii = 1:1:nPix

  fprintf('working on %d of %d\n', ii, nPix)

  % extract theta, phi
  theta = tp{ii}(1);
  phi = tp{ii}(2);
    
  vPix(ii) = sphericalharmonic(l, m, theta, phi);

end 

% plot ylm on sphere
figure(1)
subplot(1,2,1)
hp3d(real(vPix))
view(0,0);
title('real(Ylm)')
subplot(1,2,2)
hp3d(real(vPix))
view(180,0);
title('real(Ylm)')

% make projection
figure(2)
xyz = ang2vec(tp);
m_mollweide(xyz, real(vPix), 'colorbar', colorbar('southoutside'))
title('real(Ylm)')
filename = ['Y_' num2str(l) num2str(m) '_pix' num2str(nPix) '.jpg'];
print('-djpeg', filename)

% save pixel data to file
hpPix = vPix;
hcPix = zeros(size(vPix));
filename = ['Y_' num2str(l) num2str(m) '_pix' num2str(nPix) '.mat'];
save(filename, 'vPix', 'hpPix', 'hcPix')

return

