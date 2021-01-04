function periodicmotion(p)
%
% numerically integrate 1-D equation for periodic motion
% where F(x) = -sign(x) k |x|^p
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% some parameters
x0 = 1; % max displacement 
m = 1;  % mass

% chose k to give SHM frequency f0 (p=1)
f0 = 440;
T0 = 1/f0;
k = m*(2*pi*f0)^2;

% discrete times
tmin = 0;
tmax = 4; % sec
Fs = 22050; % sampling rate for audio file
deltaT = 1/Fs; %
Q = 10; % oversample when constructing waveform
dt = deltaT/Q;
N = tmax/dt;
t = linspace(tmin, tmax, N);

% initialize variables
x = zeros(N,1);
v = zeros(N,1);
x(1) = x0;
v(1) = 0;

% numerically integrate F = mxddot
for ii=2:N

  % evaluate acceleration from force
  a = -sign(x(ii-1))*(k/m)*abs(x(ii-1))^p;

  % update x and v
  v(ii) = v(ii-1) + a*dt;
  x(ii) = x(ii-1) + v(ii)*dt;

end

% determine actual period from data
ndx = find(x<0);
T = 4*t(ndx(1));

% rescale t and Fs so that period is T0 
t = t*T0/T;
Fs = Fs*T/T0;

% resample to Fs to play as audio file
ts = resample(t, 1, Q);
xs = resample(x, 1, Q);
p = audioplayer(xs, Fs);
play(p)

% display waveform (first 4 cycles)
figure(1)
subplot(2,1,1)
ndx = find(t>4*T0);
plot(t(1:ndx(1)),x(1:ndx(1)))
xlabel('time (sec)')
ylabel('displacement')

% display power spectrum
[f, P] = fourieranalyze(ts, xs, Fs);
subplot(2,1,2);
loglog(f, P);
xlabel('freq (Hz)');
ylabel('power');

pause(3);

return

