% script for calculating match as a function of time for noiseless injections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndays = [1/12 1/6 1/3 1 2 3 4 5 6 7]; % total number of days of observation

% parameters defining search
f = 100; % frequency (Hz)
res = 10; % (corresponds to N_pix=768) 
Npix = 768;
Nd = 6; % number of detectors
noiseless = 1; % noiseless injection

% simulate
for ii=1:length(Ndays)
  fprintf('simulate: working on %d of %d\n', ii, length(Ndays))
  Nt = Ndays(ii)*60; % 60 measurements per sidereal day
  basis_skies_pixels_ifo_noise(f, res, Nd, Nt, 'both');
end

% recover
for ii=1:length(Ndays)
  fprintf('recovery: working on %d of %d\n', ii, length(Ndays))
  Nt = Ndays(ii)*60; % 60 measurements per sidereal day
  match_tot(ii) = recover_skies_pixels_ifo_noise('pointsource', 'both', 'response', Npix, Nd, Nt, noiseless);
end

% plot match as a function of total number of days of observation
figure
plot(Ndays, match_tot, 'k*-')
xlabel('number of days');
ylabel('match');
xlim([0 Ndays(end)])
ylim([0 1.001])
print -depsc2 'match_vs_time.eps'

return
