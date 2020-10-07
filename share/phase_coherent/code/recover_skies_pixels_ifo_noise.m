function match_tot = ...
    recover_skies_pixels_ifo_noise(sourcetype, detectortype, regularizetype, nPix, N, Nt, noiseless)
%
% demo routine for recovering sky map in noise using ifo data saved to a file
%
% Inputs:
%   sourcetype - 'pointsource', 'grad', 'curl', 'Y_10', 'Y_20'
%   detectortype - 'rotation', 'orbit', 'both'
%   regularizetype - 'fisher', 'response'
%   N - number of ifos
%   Nt - number of discrete times
%   noiseless - flag for noiseless injections (default = 0, with noise) 
%
% Output:
%   real/imag(h+/hx) written to file 
%   'max_sourcetype_detectortype_pix_nPix_ifo_noise_N_times_Nt.jpg (and .mat)'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
try, noiseless; catch noiseless = 0; end
singleplots = 1; % flag for single plots for sky maps
seed = 1234; % random noise generation seed

% load simulated data
simfile = [sourcetype '_pix' num2str(nPix) '.mat'];
load(simfile)
h_sim = transpose([hpPix, hcPix]);

% load mapping matrix data for network of ifos
datafile = ['pix_' num2str(nPix) '_ifo_noise_' detectortype '_' num2str(N) '_times_' num2str(Nt) '.mat'];
load(datafile)

% calculate response of ifos to the background
hI = R*h_sim;

% calculate noise covariance matrix
deltaf = 0.25; % typical bin size for ground-based ifos
[Cn, invCn] = calNoiseCovariance(N, Nt, deltaf);

% calculate fisher information matrix
Fisher = R'*invCn*R;

% simulate noise 
nI = simulateNoise(Cn, seed);

% amplitude of injection
% noise level ~ 4e-24
%A = 1e-30; % no recovery
%A = 1e-26; % no recovery
%A = 1e-25; % weak recovery (barely noticeable)
%A = 3e-25; % medium recovery 
%A = 5e-25; % strong recovery 
%A = 1e-24; % very strong recovery
A = 4e-25; % 

% observed data equals ifo response plus noise
if noiseless
  s = A*hI;
else
  s = A*hI + nI;
end

% regularize either fisher matrix or response matrix
switch regularizetype
  case 'fisher'

    % do svd of Fisher matrix F    
    fprintf('working on svd\n');
    [U, S, V] = svd(Fisher);
    fprintf('finished with svd\n');

    % calculate pseudo-inverse of Fisher matrix Fisher = U S V'
    Splus = zeros(size(S'));
    M = min(size(S,1), size(S,2));
    if M==1
      Splus(1,1)=1/S(1,1);
    else
      evals = diag(S);
      figure
      semilogy(evals,'k-')
      xlabel('diagonal element')
      ylabel('singular value')
      outfile = ['evals_fisher_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.eps'];
      print('-depsc2', outfile);
 
      %ndx = floor(2*M/3); % keep only 2/3 of eigenvalues
      ndx = min(find(evals<1e-5*evals(1)))-1
      %temp = [evals(1:ndx); evals(ndx)*ones(M-ndx,1)];
      temp = [evals(1:ndx); inf(M-ndx,1)];
      Splus(1:M,1:M) = diag(1./temp);
      %Splus(1:M,1:M) = diag(1./diag(S)); % no regularization
    end
    invFisher = V*Splus*U';

    % calculate maximum-likelihood estimate
    h_ML = invFisher*(R'*invCn*s); 
    hp_ML = h_ML(1:end/2);
    hc_ML = h_ML(end/2+1:end);

    % calculate uncertainty
    sigma = sqrt(real(diag(invFisher*Fisher*invFisher)));
    sigmap = sigma(1:end/2);
    sigmac = sigma(end/2+1:end);

  case 'response'

    % do svd of response matrix R
    fprintf('working on svd\n');
    L = chol(invCn, 'lower');
    Rbar = L'*R;
    [U, S, V] = svd(Rbar);
    fprintf('finished with svd\n');
  
    % calculate pseudoinverse of pseudoinverse of response matrix
    Splus = zeros(size(S'));
    M = min(size(S,1), size(S,2));
    if M==1
      Splus(1,1)=1/S(1,1);
    else
      evals = diag(S);
      figure
      semilogy(evals,'k-')
      xlabel('diagonal element')
      ylabel('singular value')
      outfile = ['evals_response_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.eps'];
      print('-depsc2', outfile);

      %ndx = floor(2*M/3); % keep only 2/3 of eigenvalues
      ndx = min(find(evals<1e-5*evals(1)))-1
      if isempty(ndx)
        fprintf('no regularization needed\n')
        Splus(1:M,1:M) = diag(1./diag(S)); 
      else
        fprintf('regularization needed\n')
        %temp = [evals(1:ndx); evals(ndx)*ones(M-ndx,1)];
        temp = [evals(1:ndx); inf(M-ndx,1)];
        Splus(1:M,1:M) = diag(1./temp);
      end
    end
    invRbar = V*Splus*U';
 
    % calculate maximum-likelihood estimate
    h_ML = invRbar*L'*s;
    hp_ML = h_ML(1:end/2);
    hc_ML = h_ML(end/2+1:end);

    % calculate uncertainty
    var_ML = invRbar*invRbar';
    var_ML = 0.5*(var_ML+var_ML'); % force var_ML to be hermitian
    sigma = sqrt(diag(var_ML));
    sigmap = sigma(1:end/2);
    sigmac = sigma(end/2+1:end);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate match values
h_inj  = A*h_sim;
hp_inj = h_inj(1:end/2);
hc_inj = h_inj(end/2+1:end);

match_realhp = match(real(hp_ML), real(hp_inj));
match_imaghp = match(imag(hp_ML), imag(hp_inj));
match_realhc = match(real(hc_ML), real(hc_inj));
match_imaghc = match(imag(hc_ML), imag(hc_inj));
match_tot    = match(h_ML, h_inj);

fprintf('match (real h+) = %g\n', match_realhp);
fprintf('match (imag h+) = %g\n', match_imaghp);
fprintf('match (real hx) = %g\n', match_realhc);
fprintf('match (imag hx) = %g\n', match_imaghc);
fprintf('match (total)   = %g\n', match_tot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate network snr value (NOT USEFUL!!!!)
%r_ML = R*h_ML; 
%SNR = sqrt(r_ML'*invCn*r_ML) / sqrt(length(r_ML));
%fprintf('SNR = %g\n', SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save all data to .mat file
if noiseless
  filename = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.mat'];
else
  filename = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.mat'];
end
save(filename)

%return
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recovered map
if singleplots
  figure
else
  figure
  subplot(2,2,1) 
end
m_mollweide(xyz, real(hp_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('real(h+)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['max_realhp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['max_realhp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,2) 
end
m_mollweide(xyz, imag(hp_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('imag(h+)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['max_imaghp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['max_imaghp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,3)
end
m_mollweide(xyz, real(hc_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('real(hx)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['max_realhc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['max_realhc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,4) 
end
m_mollweide(xyz, imag(hc_ML), 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('imag(hx)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['max_imaghc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['max_imaghc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
else
  % print figure to file
  if noiseless
    outfile = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['max_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensitivity map

if singleplots
  figure
else
  figure
  subplot(2,1,1) 
end
m_mollweide(xyz, sigmap, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('sigma(h+)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['sigma_hp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['sigma_hp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,1,2) 
end
m_mollweide(xyz, sigmac, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('sigma(hx)')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['sigma_hc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['sigma_hc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
else
  % print figure to file
  if noiseless
    outfile = ['sigma_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['sigma_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% snr map

if singleplots
  figure
else
  figure
  subplot(2,2,1) 
end
m_mollweide(xyz, real(hp_ML)./sigmap, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('snr(real(h+))')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['snr_realhp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['snr_realhp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,2) 
end
m_mollweide(xyz, imag(hp_ML)./sigmap, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('snr(imag(h+))')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['snr_imaghp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['snr_imaghp_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,3)
end
m_mollweide(xyz, real(hc_ML)./sigmac, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('snr(real(hx))')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['snr_realhc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['snr_realhc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
end

if singleplots
  figure
else
  subplot(2,2,4) 
end
m_mollweide(xyz, imag(hc_ML)./sigmac, 'colorbar', colorbar('southoutside'))
if exist('thetaI')
  hold on
  show_ifos(thetaI, phiI)
end
title('snr(imag(hx))')

if singleplots
  % print figure to file
  if noiseless
    outfile = ['snr_imaghc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['snr_imaghc_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
print('-djpeg', outfile);
else
  % print figure to file
  if noiseless
    outfile = ['snr_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noiseless_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  else
    outfile = ['snr_' sourcetype '_' detectortype '_pix_' num2str(nPix) '_ifo_noise_' num2str(N) '_times_' num2str(Nt) '.jpg'];
  end
  print('-djpeg', outfile);
end

return

