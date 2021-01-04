% script to illustrate aliasing for a sine wave sampled too slowly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for sine wave
T0 = 1; % period of sine wave (sec)
f0 = 1/T0; % frequency of sine wave (Hz)
numPeriods = 32;
T = numPeriods*T0; % total observation time 
deltaF = 1/T;

% nyquist and sampling periods
% (aliasing when f0>fN; power aliased back to 2*fN-f0)
fN = [2*f0 0.75*f0];
fS = 2*fN;
deltaT = 1./fS;

% 'smooth' reference sine wave
deltaT0 = T0/64; % 64 points per cycle
N0 = floor(T/deltaT0);
tlow = 0;
t0 = tlow + deltaT0*[0:N0-1]; 
x0 = sin(2*pi*f0*t0);

for ii=1:length(fS)

  % discrete times
  N= floor(T/deltaT(ii));
  t = tlow + deltaT(ii)*[0:N-1]';

  % discretised sine wave
  x= sin(2*pi*f0*t);

  % discrete frequencies
  f = -fN(ii)+deltaF*[0:N-1]';

  % plot time-series
  figure(1)
  subplot(2,2,2*ii-1)
  plot(t, x, 'b*-', t0, x0, 'r-');
  xlabel('Time (sec)', 'fontsize', 14)
  ylabel('Amplitude', 'fontsize', 14);
  xlim([0 4*T0]); % plot only 4 periods

  % plot dft
  dft_x = fftshift(fft(x));
  subplot(2,2,2*ii)
  plot(f, abs(deltaT(ii)*dft_x).^2, 'b-x');
  vline(f0, 'r-');
  vline(-f0, 'r-');
  xlabel('f (Hz)','fontsize',14);
  ylabel('Abs square of DFT','fontsize',14);
  xlim([-max(fN) max(fN)]); % put on same frequency scale

end

print -depsc2 -f1 aliasing.eps

