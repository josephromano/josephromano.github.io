function auralharmonics(type)
%
% illustrate non-linear response of the human ear
%
% type = 'pure' or 'complex'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% discrete times
N = 1e5;
T = 1; % period
tmax = 100*T;
t = linspace(0,tmax,N);
fs = 1/(t(2)-t(1));

% select input wave
A = 1;
switch type
  case 'pure'
    p = A*sin(2*pi*t/T);
  case 'complex'
    p = A*sin(3*2*pi*t/T)+A*sin(5*2*pi*t/T);
  otherwise
    error('unrecognized type\n')
end

% non-linear response
a0 = 0;
a1 = 1;
a2 = 1/2;
a3 = 1/3;
x = a0 + a1*p + a2*p.^2 + a3*p.^3;
 
figure(1)
subplot(2,1,1)
plot(t,p,'k')
xlim([0 4*T])
xlabel('t (sec)')
ylabel('p')
grid on

subplot(2,1,2)
plot(t,x,'k')
xlim([0 4*T])
xlabel('t (sec)')
ylabel('x')
grid on

filename = ['auralharmonicstime_' type '.eps'];
print('-depsc2', filename)

% fourier decompose
[f, Pp] = fourieranalyze(t, p, fs);
[f, Px] = fourieranalyze(t, x, fs);

figure(2)
subplot(2,1,1)
loglog(f, Pp, 'k');
xlim([1e-1 20])
xlabel('freq (Hz)');
ylabel('power in p');
grid on
  
subplot(2,1,2)
loglog(f, Px, 'k');
xlim([1e-1 20])
xlabel('freq (Hz)');
ylabel('power in x');
grid on
  
filename = ['auralharmonicsfreq_' type '.eps'];
print('-depsc2', filename)

return
