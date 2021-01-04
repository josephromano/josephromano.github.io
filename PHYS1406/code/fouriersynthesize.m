function fouriersynthesize(amplitude, phase)
%
% fourier synthesizer using amplitude and phase 
% information from first N harmonics
%
% e.g, amplitude = [1 1/2 1/3 1/4 1/5]
%      phase     = [0 0 0 0 0]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

N = length(amplitude);

figure
grid on
xlimmin = 0;
xlimmax = 1;
ylimmin = -2;
ylimmax = 2;
xlim([xlimmin xlimmax])
ylim([ylimmin ylimmax])
xlabel('x or t');
ylabel('y');
grid on

% discrete times
tmin = 0;
tmax = 1; % sec
Nt = 1000;
t = linspace(tmin, tmax, Nt);
T = tmax;

% calculate harmonics and summed waveform
h = zeros(N,Nt);
y = zeros(1,Nt);
for n=1:N
  h(n,:) = amplitude(n)*sin(n*2*pi*t/T + phase(n)*pi/180);
  y(1,:) = y(1,:) + h(n,:);
  plot(t, h(n,:), 'LineWidth', 1);
  hold on
end

plot(t, y, 'k', 'LineWidth', 3);
xlim([xlimmin xlimmax])
ylim([ylimmin ylimmax])
xlabel('x or t');
ylabel('y');
grid on

return

