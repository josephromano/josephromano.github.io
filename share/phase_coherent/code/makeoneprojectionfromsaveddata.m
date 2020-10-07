function makeoneprojectionfromsaveddata(filename, type, thetaI, phiI)
%
% make a mollweide projection map of real(h+), imag(h+), real(hx), or imag(hx)
% for a stochastic background using data previously saved in a file ('something.mat')
%
% the file contains healpix data saved in the variables
%   vPix, hpPix, hcPix
%
% type         - 'realhp', 'imaghp', 'realhc', 'imaghc', 'power'
% thetaI, phiI - optional argument containing array of pulsar locations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filename);
nPix = length(vPix);
nSide = nPix2nSide(nPix);
tp = pix2ang(nSide);
xyz = ang2vec(tp);

switch type

  case 'realhp'

    figure
    m_mollweide(xyz, real(hpPix), 'colorbar', colorbar('southoutside'));
    if exist('thetaI')
      hold on
      show_pulsars(thetaI, phiI)
    end
    title('real(h+)')

  case 'imaghp'

    figure
    m_mollweide(xyz, imag(hpPix), 'colorbar', colorbar('southoutside'));
    if exist('thetaI')
      hold on
      show_pulsars(thetaI, phiI)
    end
    title('imag(h+)')

  case 'realhc'

    figure
    m_mollweide(xyz, real(hcPix), 'colorbar', colorbar('southoutside'));
    if exist('thetaI')
      hold on
      show_pulsars(thetaI, phiI)
    end
    title('real(hx)')

  case 'imaghc'

    figure
    m_mollweide(xyz, imag(hcPix), 'colorbar', colorbar('southoutside'));
    if exist('thetaI')
      hold on
      show_pulsars(thetaI, phiI)
    end
    title('imag(hx)')

  case 'power'

    figure
    m_mollweide(xyz, vPix, 'colorbar', colorbar('southoutside'));
    if exist('thetaI')
      hold on
      show_pulsars(thetaI, phiI)
    end
    title('|h+|^2 + |hx|^2')

  otherwise
    error('unrecognized type')

end

return

