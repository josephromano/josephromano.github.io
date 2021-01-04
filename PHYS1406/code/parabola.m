% fourier series expansion of parabola
h = 0.1;
L = 1;
N = 1000;
x = linspace(0, 1, N);
yp = h*(1- (x-L/2).^2/((L/2)^2));

n_max = 10;
y = zeros(1,N);
for n=1:n_max
  if mod(n,2)==1
    An = 32*h/(n^3 * pi^3);
    y = y + An*sin(n*pi*x/L);
  end
end  

figure(1)
plot(x, yp, 'o', x, y, 'LineWidth', 2)
return  
