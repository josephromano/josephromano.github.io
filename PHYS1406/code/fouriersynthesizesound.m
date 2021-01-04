function fouriersynthesizesound(amplitude, phase)
%
% fourier synthesize sound using amplitude and phase 
% information from first N harmonics
%
% e.g, amplitude = [1 1/2 1/3 1/4 1/5]
%      phase     = [0 0 0 0 0]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

N = length(amplitude);

f = 220; % concert A2
tmin = 0;
tmax = 4; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
Nt = floor(tmax/deltaT)+1;
t = linspace(tmin, tmax, Nt);
T = 1/f;

% calculate harmonics and summed waveform
h = zeros(N, Nt);
y = zeros(1, Nt);

for n=1:N
  h(n,:) = amplitude(n)*sin(n*2*pi*f*t + phase(n)*pi/180);
  y(1,:) = y(1,:) + h(n,:);
end

% plot 4 cycles of summed waveform
figure
subplot(2,1,1);
ndx = find(t>=4*T);
plot(t(1:ndx(1))*1000, y(1:ndx(1)))
xlabel('time (ms)');
ylabel('y(t)');
grid on

[f, P] = fourieranalyze(t, y, Fs);
subplot(2,1,2);
loglog(f, P);
xlabel('freq (Hz)');
ylabel('power');
grid on

% play sound
p = audioplayer(y, Fs);
play(p); 
pause(5)

return

