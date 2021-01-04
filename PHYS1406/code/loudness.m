% plot loudness versus intensity curves
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold of hearing
I0 = 10^(-12); % W/m^2
p0 = 2*10^(-5); % Pa

% threshold of pain
Ipain = 1; % W/m^2
ppain = 20; % Pa

% discrete pressure and intensity values
p = logspace(log10(p0), log10(ppain), 100);
I = logspace(log10(I0), log10(Ipain), 100);

% Lp(=LI), LL, S
Lp = 20*log10(p/p0);
LI = 10*log10(I/I0);
LL = Lp; % at f=1000 Hz
S = 2.^((LL-40)/10);

% make plots
figure(1)
semilogx(I, LL, 'k-', I, S, 'r-');
xlim([I(1) I(end)])
ylim([S(1) S(end)])
grid on
xlabel('intensity at 1000 Hz(W/m^2)')
ylabel('loudness (phon or sone)')
legend('L_L', 'S', 'location', 'SouthEast');
print('-depsc2', 'loudness-semilog.eps')

figure(2)
loglog(I, LL, 'k-', I, S, 'r-');
xlim([I(1) I(end)])
ylim([S(1) S(end)])
grid on
xlabel('intensity at 1000 Hz (W/m^2)')
ylabel('loudness (phon or sone)')
legend('L_L', 'S', 'location', 'SouthEast');
print('-depsc2', 'loudness-loglog.eps')

return

