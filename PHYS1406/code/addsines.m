function addsines(phi)
%
% add two sine waves with the same amplitude, frequency, 
% but with a phase offset phi (in degrees)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% parameters for sine waves
A = 1; % amplitude
%f = 440; % frequency of sine wave (Hz)
f = 1; % frequency of sine wave (Hz)
T = 1/f; % period (sec)
numPeriods = 2;
Tobs = numPeriods*T; % total observation time 

% discrete times
N = 500; 
t = linspace(0, Tobs, N);

% sine waves
y1 = A*sin(2*pi*f*t);
y2 = A*sin(2*pi*f*t + phi*pi/180);
y  = y1+y2;
 
% plot functions
figure
plot(t, y1, '-k', t, y2, '-b', t, y, '-r');
grid on
xlabel('time (s)');
titlestr = ['Phase shift = ' num2str(phi) ' degrees'];
title(titlestr)
legend('y1', 'y2', 'y1+y2');

filename = 'addsines.eps';
print('-depsc2',  filename);

return

