function bowedstring(alpha)
%
% illustrate the physics of a bowed violin string
%
% alpha = fractional distance from bridge where string is bowed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% check that 0 < alpha < 1
if alpha <=0 | alpha>=1 
  error('improper location for bowing');
end

% assign certain parameter values
nmax = 25; % max number of harmonics
L = 1; % length of string
T = 1; % period 
v = 2*L/T; % velocity of wave on string
h = 0.1; % max displacement of string (at center)
y0 = h*(1-(alpha*L-L/2).^2/((L/2)^2));
amp = 2*y0/(pi^2*alpha*(1-alpha));

% x-values 
N = 125;
x = linspace(0,L,N);
ndx = floor(alpha*N);
x0 = x(ndx);

% discrete time steps
Ttot = 2*T; % 2 periods
dt = T/(2*N+2); % 2(N+1) samples / period
numT = floor(Ttot/dt);
t = linspace(0, Ttot, numT);

% calculate envelope (parabola)
yp = h*(1-(x-L/2).^2/((L/2)^2));
ym = -yp;

% calculate fourier coefficients for waveform and 
% sideways force on bridge
A = zeros(nmax);
B = zeros(nmax);
for n=1:nmax
  A(n) = amp * (-1)^(n+1) * (1/n^2);
  B(n) = (-1)^(n+1) * (1/n); % ignoring normalization
end

% plot fourier coefficients for waveform
figure(1)
subplot(2,1,2)
bar(A/max(A));
xlabel('harmonic number')
pause(0.1)

% construct waveform, displacement, and sideways force on bridge
D = zeros(numT); % displacement of string at bow location
force = zeros(numT); % sideways force of string on bridge
m=0; % counter
for ii=1:numT

  y = zeros(size(x));
  z = 0;
  for n=1:nmax
   
    % waveform 
    y = y + A(n) * sin(n*pi*x/L) * sin(n*2*pi*t(ii)/T);

    % sideways force on bridge (ignorning normalization)
    % (for small displacements, this is proportional to the slope = dy/dx at x=0)
    z = z + B(n) * sin(n*2*pi*t(ii)/T);

  end;

  % displacement
  D(ii) = y(ndx);

  % sideways force
  force(ii) = z;

  % update integer counter if necessary
  if 2*t(ii)/T > (2*m+1)-alpha
    m=m+1;
  end

  subplot(2,1,1)
  if (2*m-1)+alpha <= 2*t(ii)/T & 2*t(ii)/T <= (2*m+1)-alpha
    v_up = 2*y0/((1-alpha)*T);
    plot(x0, v_up*(t(ii)-m*T), 'or', x, y, 'k', x, yp, 'g', x, ym, 'g');
  else
    v_down = -2*y0/(alpha*T);
    plot(x0, y0+v_down*(t(ii)-((2*(m-1)+1-alpha)/2)*T), 'ob', x, y, 'k', x, yp, 'g', x, ym, 'g');
  end
  vline(x0)

  xlabel('x')
  ylabel('y')
  axis equal
  xlim([0 1])
  ylim([-h h])
  pause(.0001)

end

figure(2)
plot(t,D)
ylabel('displacement')
xlabel('time')
grid on

figure(3)
subplot(2,1,1)
plot(t,force)
ylabel('force')
xlabel('time')
grid on

subplot(2,1,2)
bar(B/max(B));
xlabel('harmonic number')

return
