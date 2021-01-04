function forcedoscillator(F)
%
% illustrate chaotic motion of driven oscillator
%
% F - dimensionless driving torque (N_d/(m l^2 omega_0^2)
%
% differential equations:
%
%  dx/dt' = y
%  dy/dt' = -cy - sinx + F cos(omega t')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F = 0.4; % simple harmonic motion
%F = 0.5; % simple harmonic motion
%F = 0.6; % chaotic
%F = 0.7; % chaotic
%F = 0.8; % period 1
%F = 0.9; % period 2
%F = 1.0; % chaotic

close all

x0 = pi/6; % initial angle
y0 = -1; % initial angular velocity (y=dx/dt')

c = 0.05; % b/(m l^2 omega_0)
omega = 0.7; % omega_d/omega_0

N_cycles = 15000; % number of cycles 
N_steps = 1000; % number of steps/cycle
N_tot = N_cycles*N_steps;

dtp = 2*pi/(omega*N_steps);
tp = dtp*[0:1:N_tot-1]';
z = (2*pi/N_steps)*[0:1:N_tot-1]';  % z = omega*tp;

% initialize arrays
x = zeros(N_tot, 1);
y = zeros(N_tot, 1);
x(1) = x0;
y(1) = y0;

% evolve equations
for ii=0:N_cycles-1
  for jj=1:N_steps;
  
    kk = ii*N_steps + jj;

    dx = y(kk)*dtp;
    dy = (-c*y(kk) - sin(x(kk)) + F*cos(z(kk)))*dtp;

    % update variables
    x(kk+1) = x(kk) + dx;
    y(kk+1) = y(kk) + dy;

  end
end

% remove last value
x = x(1:end-1);
y = y(1:end-1);

% make angle variable run from -pi to pi
xmod = mod(x+pi,2*pi)-pi;

% extract poincare map
ind = [1:N_steps:N_tot];
x_poincare = xmod(ind);
y_poincare = y(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot variables
f1=figure('Position',[300 400 900 640]);

subplot(2,3,1)
plot(z(end-4*N_steps:end), y(end-4*N_steps:end))
xlabel('z = \omega tprime', 'FontSize', 12);
ylabel('y = dx/dtprime', 'FontSize', 12);
xlim([z(end-4*N_steps) z(end)]);
yl = 1.2*ylim;
ylim(yl);
title('Time series', 'FontSize', 12);

subplot(2,3,2)
plot(xmod(end-1000*N_steps:end), y(end-1000*N_steps:end), 'k.')
xlim([-pi pi])
ylim(yl);
xlabel('Angle x', 'FontSize', 12);
ylabel('y = dx/dtprime', 'FontSize', 12);
title('Phase space plot', 'FontSize', 12);

subplot(2,3,3)
plot(x_poincare(end-1000:end), y_poincare(end-1000:end), 'k.')
xlim([-pi pi])
ylim(yl);
xlabel('Angle x', 'FontSize', 12);
ylabel('y = dx/dtprime', 'FontSize', 12);
title('Poincare section', 'FontSize', 12);

subplot(2,3,5)
h5=plot(z(end-4*N_steps:end), xmod(end-4*N_steps:end));
xlim([z(end-4*N_steps) z(end)]);
ylim([-pi pi])
xlabel('z = \omega tprime', 'FontSize', 12);
ylabel('Angle x', 'FontSize', 12);
title('Time series', 'FontSize', 12);
view(-90, 90)

filename=['forcedoscillator_F_' num2str(F) '.eps'];
print('-depsc2', filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot motion (last 15 cycles)

figure(2)
for ii=N_tot-15*N_steps:floor(N_steps/100):N_tot
  xx = sin(x(ii));
  yy = -cos(x(ii)); 
  plot(xx,yy,'bo')
  line([0 0], [0 -1], 'Color', 'g');
  line([0 xx], [0 yy], 'Color', 'k');
  axis equal
  xlim([-1 1])
  ylim([-1 1])
  titstr = ['F = ' num2str(F) ', cycle = ' num2str(ceil(ii/N_steps))];
  title(titstr); 
  pause(0.00001);
end

return
