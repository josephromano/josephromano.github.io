function A = fourierdecompose(y0, x)
%
% calculate fourier coefficients for plucked guitar string
% fixed at both ends
%
% (initial velocity = 0, initial configuration = y0(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(y0);
L = x(end);
dx = x(2)-x(1);

A = zeros(N-2,1);
for n=1:N-2
  A(n) = 0;
  for jj=1:N
    A(n) = A(n) + (2/L)*dx*y0(jj)*sin(n*pi*x(jj)/L);
  end

end

return

