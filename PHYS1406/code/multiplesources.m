% plot various quantities versus number of sources
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of sources
N = [1:.01:100];

% value for 1 source
p1 = 1; % pressure
I1 = 1; % intensity
L1 = 1; % sound intensity (or pressure) level
S1 = 1; % subjective loudness level

% value for multiples sources
ptot = sqrt(N)*p1;
Itot = N*I1;
Ltot = 10*log10(N)+L1; 
Stot = 2.^((Ltot-L1)/10)*S1; % Stot = 2.^log10(N)*S1;

% make plots
%semilogx(N, Ltot, 'k-', N, Itot, 'm-', N, ptot, '-b', N, Stot, 'r-');
loglog(N, Ltot, 'k-', N, Itot, 'm-', N, ptot, '-b', N, Stot, 'r-');
xlim([N(1) N(end)])
grid on
xlabel('N (number of sources)')
legend('L_L', 'I', 'p', 'S', 'location', 'northwest');
print('-depsc2', 'multiplesources.eps')

return

