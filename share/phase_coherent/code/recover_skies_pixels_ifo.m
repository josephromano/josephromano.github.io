function recover_skies_pixels_ifo(sourcetype, detectortype, nPix, N, Nt)
%
% demo routine for recovering sky map using ifo data saved to a file
%
% Inputs:
%   sourcetype - 'pointsource', 'grad', 'curl', 'Y_10', 'Y_20'
%   detectortype - 'rotation', 'orbit', 'both'
%   N=number of ifos
%   Nt=number of discrete times
%
% Output:
%   real/imag(h+/hx) written to file 
%   'max_sourcetype_detectortype_pix_nPix_ifo_N_times_Nt.jpg (and .mat)'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% load simulated data
simfile = [sourcetype '_pix' num2str(nPix) '.mat'];
load(simfile)
h_sim = transpose([hpPix, hcPix]);

% load mapping matrix data for network of ifos
datafile = ['pix_' num2str(nPix) '_ifo_' detectortype '_' num2str(N) '_times_' num2str(Nt) '.mat'];
load(datafile)
Sigma = S;

% calculate response of ifos to the background
hI = H*h_sim;

% no noise injection
s = hI;

% calculate pseudo-inverse of Sigma
SigmaPlus = zeros(size(Sigma'));
M = min(size(Sigma,1), size(Sigma,2));
if M==1
  SigmaPlus(1,1)=1/Sigma(1,1);
else
  SigmaPlus(1:M,1:M) = diag(1./diag(Sigma));
end

% calculate maximum-likelihood estimate
h_ML = V*SigmaPlus*U'*s; 
hp_ML = h_ML(1:end/2);
hc_ML = h_ML(end/2+1:end);

% save all data to .mat file
filename = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_' num2str(N) '_times_' num2str(Nt) '.mat'];
save(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make maps of recovered background, including ifo locations
thetaI = zeros(N,1);
phiI = zeros(N,1);
for ii=1:N
  x = det{ii}.r(1);
  y = det{ii}.r(2);
  z = det{ii}.r(3);
  RE = sqrt(sum(r.^2));
  thetaI(ii) = acos(z/RE);
  phiI(ii) = atan2(y,x);
end

xyz = ang2vec(tp);
figure
subplot(2,2,1) 
m_mollweide(xyz, real(hp_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('real(h+)')

subplot(2,2,2) 
m_mollweide(xyz, imag(hp_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('imag(h+)')

subplot(2,2,3)
m_mollweide(xyz, real(hc_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('real(hx)')

subplot(2,2,4) 
m_mollweide(xyz, imag(hc_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('imag(hx)')

% print figure to file
outfile = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_' num2str(N) '_times_' num2str(Nt) '.jpg'];
print('-djpeg', outfile);

%%%%%%%%%%%%%%%%%%%
%figure
%power_ML = abs(hp_ML).^2 + abs(hc_ML).^2;
%m_mollweide(xyz, power_ML, 'colorbar', colorbar('southoutside'));
%if exist('thetaI')
%  hold on
%  show_ifos(thetaI, phiI)
%end
%title('|h+|^2+|hx|^2')

% print figure to file
%outfile = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_' num2str(N) '_times_' num2str(Nt) '_power.jpg'];
%print('-djpeg', outfile);

return

