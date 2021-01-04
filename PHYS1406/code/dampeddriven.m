function dampeddriven(damping, A, f_ratio)
%
% numerically integrate 1-D equation for damped, driven
% oscillations
%
% damping: magnitude of damping force in units omega0
%          (<1, =1, >1 for under, critical, and over damping)
% A:       amplitude of driving force in units of initial displacement x0
% f_ratio: frequency of driving force in units of omega0
%
% e.g., damping = 0.1, A=0.2, f_ratio = 0.2 or 5 or 1 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% some parameters
x0 = 1; % max displacement (m)
f0 = 1; % natural frequency (Hz)
omega0 = 2*pi*f0; 
F0overm = A*omega0^2*x0; % driving force F0 over m (m/s^2)
beta = damping*omega0; % beta = b/2m (Hz); F = -bv (damping force)

% driving frequency (Hz)
f = f_ratio*f0;
omega = 2*pi*f;

% display resonant frequency and maximum displacement on resonance
if omega0^2>2*beta^2
  omegaR = sqrt(omega0^2 - 2*beta^2); % resonant frequency
  %DR = F0overm/(2*beta*sqrt(omega0^2 - beta^2)); 
  DR = F0overm/sqrt((omega0^2 - omegaR^2)^2 + 4*beta^2*omegaR^2);
  fprintf('resonant frequency = %f\n', omegaR/(2*pi));
  fprintf('steady state max displacement on resonance = %f\n', DR);
else
  fprintf('resonant frequency does not exist; damping too large\n')
end

% maximum displacement
D = F0overm/sqrt((omega0^2 - omega^2)^2 + 4*beta^2*omega^2);
fprintf('steady state max displacement = %f\n', D);

% discrete times
tmin = 0;
tmax = 25; % sec
Fs = 100*max(f0, f); % sample frequency (Hz)
dt = 1/Fs;
N = floor(tmax/dt);
t = linspace(tmin, tmax, N);

% driving force motion
xdriving = A*cos(omega*t);

% initialize variables
x = zeros(N,1);
v = zeros(N,1);
x(1) = x0;
v(1) = 0;

% numerically integrate F = mxddot
for ii=2:N

  % evaluate acceleration from force
  a = F0overm*cos(omega*t(ii)) - 2*beta*v(ii-1) - omega0^2*x(ii-1);
  v(ii) = v(ii-1) + a*dt;
  x(ii) = x(ii-1) + v(ii)*dt;

end

% display waveform 
figure
plot(t, x, 'b', t, xdriving, 'r')
xmax = max(abs(x));
ylim([-xmax xmax])
grid on
xlabel('time (sec)')
ylabel('displacement')
legend('mass motion', 'driving force', 'location', 'northeast')

return

