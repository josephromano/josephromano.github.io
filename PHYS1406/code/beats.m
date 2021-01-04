function beats(f1, f2)
%
% add two sine waves with the same amplitude and phase,
% but different frequencies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
try, f1; catch, f1 = 9; end;
try, f2; catch, f2 = 11; end;

% parameters for sine waves
A = 1; % amplitude
%f = 440; % frequency of sine wave (Hz)
f0 = (f1 + f2)/2;

% beat frequency
fb = abs(f1-f2);
Tb = 1/fb;

% second-order beats
fb2 = abs(2*min(f1,f2) - max(f1,f2));
Tb2 = 1/fb2;

% total observation time
numPeriods = 2;
Tobs = numPeriods*max(Tb,Tb2); % total observation time 

% discrete times
N = 50*Tobs*f0; 
t = linspace(0, Tobs, N);

% sine waves
y1 = A*sin(2*pi*f1*t);
y2 = A*sin(2*pi*f2*t);
y  = y1+y2;
 
% plot functions
figure
subplot(2,1,1)
plot(t, y1, '-k', t, y2, '-b');
grid on
legend('y1', 'y2');
titlestr = ['f_1 = ' num2str(f1) ' Hz; f_2 = ' num2str(f2)  ' Hz'];
title(titlestr)

subplot(2,1,2)
plot(t, y, '-r');
grid on
legend('y1+y2');
xlabel('time (s)');

filename = 'beats.eps';
print('-depsc2',  filename);

return

