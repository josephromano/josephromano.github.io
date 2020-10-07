function simulateBackground(res)
%
% routine for simulating total, grad, and curl components of a
% GW background
%
% res - angular resolution (e.g., res=7 -> 3072 pixels)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% set other parameter values
lG = 10; % lmax for isotropic grad mode
lC = 5;  % lmax for isotropic curl mode
lmax_simulation = 10; % lmax for simulation of background

seedG = 1; % seed for grad mode simulation
seedC = 2;  % seed for curl mode simulation

% angular resolution (arcsec) for healpix pixels
r = res*60*60; 
nSide = res2nSide(r);
nPix = nSide2nPix(nSide); % number of modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulate stochastic background
a_G = simulate_aPlm('isoG', lG, lmax_simulation, seedG);
a_C = simulate_aPlm('isoC', lC, lmax_simulation, seedC);

% NOTE!!!! make curl component stronger
a_C.data = 2*a_C.data; 

% extract lvec, mvec for later use
lvec = a_G.lvec;
mvec = a_G.mvec;

% combine a_G, a_C and fill structure
a_T.data = a_G.data + a_C.data;
a_T.lvec = lvec;
a_T.mvec = mvec;
a_T.lmax = lmax_simulation;
  
% make projections saving pixel data to files
filename = ['tot_pix' num2str(nPix) '.mat']; 
makeprojection(a_T, res, filename);
filename = ['tot_pix' num2str(nPix) '.jpg']; 
print('-djpeg', filename);

filename = ['curl_pix' num2str(nPix) '.mat']; 
makeprojection(a_C, res, filename);
filename = ['curl_pix' num2str(nPix) '.jpg']; 
print('-djpeg', filename);

filename = ['grad_pix' num2str(nPix) '.mat']; 
makeprojection(a_G, res, filename);
filename = ['grad_pix' num2str(nPix) '.jpg']; 
print('-djpeg', filename);

return
