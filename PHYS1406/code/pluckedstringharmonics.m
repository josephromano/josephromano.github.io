function pluckedstringharmonics(alpha)
%
% determine fourier coefficients for a plucked string
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
L = 1; % length of string
h = 0.1; % max displacement of string at location of plucking

% x-values 
N = 251;
x = linspace(0,L,N);
ndx = floor(alpha*N);
x0 = x(ndx);

% construct configuration of plucked string
y = zeros(length(x), 1);
m = h/x0;
y(1:ndx) = m * x(1:ndx); % to left of x0
m = -h/(1-x(ndx+1));
y(ndx+1:end) = h + m*(x(ndx+1:end) - x(ndx+1)); % to right of x0

% plot initial configuration of string
subplot(2,1,1)
plot(x, y, 'r-');
ylabel('y');
ylim([-.15 .15])
xlim([0 1])

% calculate and plot fourier coefficients
%A = fourierdecompose(y, x);
nmax = 25;
for n=1:nmax
  A(n) = sin(n*pi*alpha)/n^2;
end
subplot(2,1,2)
bar(A(1:10)/max(A)); % first 10 fourier coefficients
ylabel('A_n');

% print to file
filename = 'pluckedstringharmonics.eps';
print('-depsc2',filename)

return
