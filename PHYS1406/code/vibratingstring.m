function  vibratingstring(y0)
%
% small transverse vibrations of a stretched string fixed at both ends
% (approximate string as N discrete mass points)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
f=figure('Position',[0 400 850 640],...
         'Name','vibrating sting');

% check that first and last values of y0 are zero
tol = 1e-9;
if abs(y0(1))>tol | abs(y0(end))>tol
  error('string not fixed at both ends');
end

% assign certain parameter values
L = 1; % length of string
m = 1; % total mass of string
tau = 1; % tension in string

% x-values corresponding to y0
x0 = linspace(0,L,length(y0));

% number of mass points (doesn't include fixed ends)
N = length(y0)-2;
y = y0(2:end-1);
x = x0(2:end-1);
dx = x(2)-x(1); % spacing between masses
dm = m/(N+1); % mass of individual particles

% derived quantities
mu = m/L; % mass per unit length of string
c = sqrt(tau/mu); % wave velocity on the string
f1 = c/(2*L); % fundamental frequency 
T1 = 1/f1; % period

% discrete time steps
Ttot = 3*T1; % 3 periods
dt = T1/(2*N+2); % 2(N+1) samples / period
numT = floor(Ttot/dt);
t = linspace(0, Ttot, numT);

% plot initial configuration of string
subplot(4,1,1)
plot(x, y, 'or');
ylabel('y0');
ylim([-.15 .15])
xlim([0 L])
for j=1:N+1
  if j==1
    line([0 x(j)], [0 y(j)], 'Color', 'k');
  elseif j==N+1
    line([x(j-1) L], [y(j-1) 0], 'Color', 'k');   
  else
    line([x(j-1) x(j)], [y(j-1) y(j)], 'Color', 'k');
  end
end

% calculate and plot fourier coefficients
A = fourierdecompose(y0, x0);
subplot(4,1,2)
bar(A(1:10)/max(A)); % first 10 fourier coefficients
ylabel('A_n');

M = input('enter number of harmonics for fourier series approximation: ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize arrays
v = zeros(N,1);
a = zeros(N,1);
yold = y;
vold = v;

for ii=1:numT

  % determine new positions of masses
  for j=1:N
    % special case for j=1 and j=N (1st and last mass)
    if j==1
      a(j) = (tau/(dx*dm))*(0 - 2*yold(j) + yold(j+1));
    elseif j==N
      a(j) = (tau/(dx*dm))*(yold(j-1) - 2*yold(j) + 0);
    else
      a(j) = (tau/(dx*dm))*(yold(j-1) - 2*yold(j) + yold(j+1));
    end

    v(j) = vold(j) + a(j)*dt;
    y(j) = yold(j) + v(j)*dt;
  end

  pause(.0001)

  subplot(4,1,3)
  plot(x, y, 'or');
  ylabel('y');
  ylim([-.15 .15])
  xlim([0 L])
  for j=1:N+1
    if j==1
      line([0 x(j)], [0 y(j)], 'Color', 'k');
    elseif j==N+1
      line([x(j-1) L], [y(j-1) 0], 'Color', 'k');   
    else
      line([x(j-1) x(j)], [y(j-1) y(j)], 'Color', 'k');
    end

  end

  % update yold, vold
  yold = y;
  vold = v;

  % evolve approximate waveform
  yapprox = zeros(1,N+2);
  e = zeros(M,N+2);  % normal modes

  for n=1:M
    e(n,:) = sin(n*pi*x0/L).*cos(n*pi*c*t(ii)/L);
    yapprox = yapprox + A(n)*e(n,:);
  end

  subplot(4,1,4)
  plot(x0, yapprox, 'or');
  ylabel('yapprox');
  ylim([-.15 .15])
  xlim([0 L])
  for j=2:N+2
      line([x0(j-1) x0(j)], [yapprox(j-1) yapprox(j)], 'Color', 'k');
  end


end

return
