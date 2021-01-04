function pluckedstring(alpha)
%
% illustrate physics of a plucked string
%
% alpha = fractional distance from bridge where string is plucked
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% check that 0 < alpha < 1
if alpha <=0 | alpha>=1 
  error('improper location for plucking');
end

% assign certain parameter values
nmax = 25; % max number of harmonics
L = 1; % length of string
T = 1; % period 
v = 2*L/T; % velocity of wave on string
h = 0.1; % max displacement of string at location of plucking

% x-values 
N = 251;
x = linspace(0,L,N);
ndx = floor(alpha*N);
x0 = x(ndx);
amp = 2*h/(pi^2*alpha*(1-alpha));

% discrete time steps
Ttot = 2*T; % 2 periods
dt = T/(2*N+2); % 2(N+1) samples / period
numT = floor(Ttot/dt);
t = linspace(0, Ttot, numT);

% initial configuration of plucked string
y0 = zeros(length(x), 1);
m = h/x0;
y0(1:ndx) = m * x(1:ndx); % to left of x0
m = -h/(1-x(ndx+1));
y0(ndx+1:end) = h + m*(x(ndx+1:end) - x(ndx+1)); % to right of x0

% plot initial configuration of plucked string
figure(1)
subplot(2,1,1)
plot(x, y0, 'k');
axis equal
xlim([0 1])
ylim([-h h])
xlabel('x')
ylabel('y');

% calculate fourier coefficients for waveform and 
% sideways force on bridge
A = zeros(nmax);
B = zeros(nmax);
for n=1:nmax
  A(n) = amp * sin(n*pi*alpha) * (1/n^2);
  B(n) = sin(n*pi*alpha) * (1/n); % ignoring normalization
end

% plot fourier coefficients for plucked string
subplot(2,1,2)
bar(A/max(A));
xlabel('harmonic number')

pause(1.5)

% construct waveform and sideways force on bridge
force = zeros(numT);
for ii=1:numT
  
  y = zeros(1,N);
  z = 0;
  for n=1:nmax

    % waveform
    y = y + A(n) * sin(n*pi*x/L) * cos(n*2*pi*t(ii)/T);

    % sideways force on bridge (ignorning normalization)
    % (for small displacements, this is proportional to the slope = dy/dx at x=0)
    z = z + B(n) * cos(n*2*pi*t(ii)/T);

  end

  force(ii) = z;

  subplot(2,1,1)
  plot(x, y, 'k');
  axis equal
  xlim([0 1])
  ylim([-h h])
  xlabel('x')
  ylabel('y');
  pause(.0001)

end

figure(2)
subplot(2,1,1)
plot(t,force)
ylabel('force')
xlabel('time')
grid on

subplot(2,1,2)
bar(B/max(B));
xlabel('harmonic number')

return
