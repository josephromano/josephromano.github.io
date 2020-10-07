% script to plot singular values for orbit-only and 
% rotation+orbital motion on the same plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
basis_skies_pixels_ifo_noise(100, 10, 1, 1200, 'both');
clear all
basis_skies_pixels_ifo_noise(100, 10, 1, 1200, 'orbit');

clear all
recover_skies_pixels_ifo_noise('grad', 'both', 'response', 768, 1, 1200);
clear all
recover_skies_pixels_ifo_noise('grad', 'orbit', 'response', 768, 1, 1200);

% rotation+orbit
clear all
load max_grad_both_pix_768_ifo_noise_1_times_1200.mat
svals1 = evals/evals(1);

% orbit-only
load max_grad_orbit_pix_768_ifo_noise_1_times_1200.mat
svals2 = evals/evals(1);

% make plot
close all
N = length(svals1);
semilogy([1:N], svals1, 'b-', [1:N], svals2, 'b--', 'linewidth', 3)
xlabel('diagonal element')
ylabel('singular value')
xlim([1 1200])
ylim([1e-20 10])

print('-depsc2', 'singvals.eps');

return

