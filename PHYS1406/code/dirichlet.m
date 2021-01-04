% matlab script to compare sinc and dirichlet kernels for different
% values of N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% label for dirichlet kernel
N = 16;

% discrete x values (x=f*Delta t)
xmax = 1;
xmin = -xmax;
numX = 1024;
deltaX = (xmax-xmin)/numX;
x = xmin+deltaX*[0:numX-1];

% dirichlet kernel D_N(x) = sin(N pi x)/(N sin(pi x))
d = sin(N*pi*x)./(N*sin(pi*x));

% sinc function sinc(N*pi*x) = sin(N pi x)/(N pi x)
s = sin(N*pi*x)./(N*pi*x);

% plot
figure(1)
subplot(2,1,1)
plot(x, s, x, d);
title('N=16','fontsize', 16); 
legend('sinc', 'dirichlet');
xlabel('x=f\Delta t','fontsize',14)
grid on;
subplot(2,1,2)
plot(x, s, x, d);; 
legend('sinc', 'dirichlet');
xlabel('x=f\Delta t','fontsize',14)
xlim([-0.5 0.5]);
grid on;

print -depsc2 -f1 dirichlet
