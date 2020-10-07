% test overlap_ifo.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
const = physConstants('mks');

site1 = 'H1';
site2 = 'L1';

% calculate overlap reduction functions for different values of lmax
[orf5, freq] = overlap_ifo(site1, site2, 5, 0, 1000);
[orf10, freq] = overlap_ifo(site1, site2, 10, 0, 1000);
[orf20, freq] = overlap_ifo(site1, site2, 20, 0, 1000);
[orf30, freq] = overlap_ifo(site1, site2, 30, 0, 1000);

% compare with overlap reduction function from standard calculation
load overlap_H1_L1_lw.mat;

% make plots
figure(1)
semilogx(freq, real(orf5*orf(1)/orf5(1)), ...
         freq, real(orf10*orf(1)/orf10(1)), ...
         freq, real(orf20*orf(1)/orf20(1)), ...
         freq, real(orf30*orf(1)/orf30(1)), ...
         f, real(orf), 'k')
xlim([0 1000])
xlabel('freq (Hz)')
ylabel('overlap(f)')
legend('lmax = 5', 'lmax = 10', 'lmax = 20', 'lmax = 30', 'standard', ...
       'location', 'SouthEast')

% print to file
outfile = ['overlap_comparison_' site1 '_' site2 '.eps'];
print('-depsc2', outfile);

figure(2)
semilogx(freq, real(orf5*orf(1)/orf5(1)), ...
         freq, real(orf10*orf(1)/orf10(1)), ...
         freq, real(orf20*orf(1)/orf20(1)), ...
         freq, real(orf30*orf(1)/orf30(1)), ...
         f, real(orf), 'k')
xlim([100 1000])
xlabel('freq (Hz)')
ylabel('overlap(f)')
legend('lmax = 5', 'lmax = 10', 'lmax = 20', 'lmax = 30', 'standard', ...
       'location', 'SouthEast')

% print to file
outfile = ['overlap_comparison_' site1 '_' site2 '_zoom.eps'];
print('-depsc2', outfile);

return

