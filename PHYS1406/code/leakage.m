% data length, padding factor, max bins, FFT length
N = 1024;
pad = 128; % zero padding for increased frequency resolution
b_max = 5; % display leakage from -bmax to b_max
N_FFT = N*pad;

% discrete times, frequency bins
t = [0:N-1]'/N;
b = [0:pad*b_max]'/pad; 
b = union(-b, b); % union so 0 isn't double counted
 
% windows
w1 = rectwin(N);
w2 = triang(N);
w3 = 1 - (2*t-1).^2;
w4 = hann(N);

% ffts
y1=fft(w1, N_FFT);
y2=fft(w2, N_FFT);
y3=fft(w3, N_FFT);
y4=fft(w4, N_FFT);

% amplitude of leakage function (normalised to unity at bin=0)
a1=abs(y1)/max(abs(y1));
a2=abs(y2)/max(abs(y2));
a3=abs(y3)/max(abs(y3));
a4=abs(y4)/max(abs(y4));

% plot window functions
figure(1)
plot(t, w1, 'b', ...
     t, w2, 'g', ...
     t, w3, 'k', ...
     t, w4, 'r');
xlim([0 1]);
ylim([0 1.02]);
xlabel('time (t/T)', 'fontsize', 14);
ylabel('amplitude', 'fontsize', 14);
legend('rectangular', 'triangular', 'Welch', 'Hann', 'location', 'northeast');
print -depsc2 -f1 windows.eps

% plot leakage functions
figure(2)
Nb = (length(b)+1)/2;
semilogy(b, [a1(end-(Nb-2):end); a1(1:Nb)], 'b', ...
         b, [a2(end-(Nb-2):end); a2(1:Nb)], 'g', ...
         b, [a3(end-(Nb-2):end); a3(1:Nb)], 'k', ...
         b, [a4(end-(Nb-2):end); a4(1:Nb)], 'r');
grid on
ylim([1e-5 1])
xlabel('offset in units of frequency bins (bin size=1/T)', 'fontsize', 14);
ylabel('amplitude of leakage', 'fontsize', 14);
legend('rectangular', 'triangular', 'Welch', 'Hann', 'location', 'northwest');
print -depsc2 -f2 leakage.eps

