function recover_skies_sph_ifo(sourcetype, detectortype, nPix, N, Nt)
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
simfile = [sourcetype '_sph' num2str(nPix) '.mat'];
load(simfile)
a_sim = a_Plm;

% load mapping matrix data for network of ifos
datafile = ['sph_' num2str(nPix) '_ifo_' detectortype '_' num2str(N) '_times_' num2str(Nt) '.mat'];
load(datafile)
Sigma = S;

% calculate response of ifos to the background
hI = H*a_sim.data;

% extract lvec, mvec for later use
lvec = a_sim.lvec;
mvec = a_sim.mvec;

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

% calculate maximum likelihood estimators
a_ML.data = zeros(size(lvec));
a_ML.data = V*SigmaPlus*U'*s; % only grad modes non-zero
a_ML.lvec = lvec;
a_ML.mvec = mvec;
a_ML.lmax = lmax;

% save all data to .mat file
filename = ['max_' sourcetype '_' detectortype '_sph_' num2str(nPix) '_ifo_' num2str(N) '_times_' num2str(Nt) '.mat'];
save(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make map of recovered background
makeprojection(a_ML, res)

% print figure to file
outfile = ['max_' sourcetype '_' detectortype '_sph_' num2str(nPix) '_ifo_' num2str(N) '_times_' num2str(Nt) '.jpg'];
print('-djpeg', outfile);

return

