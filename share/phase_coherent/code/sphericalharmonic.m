function y = sphericalharmonic(l, m, theta, phi)
%
% calculate the value of the spherical harmonic
%
% Y_lm(theta, phi) = N_lm P_l^m(cos\theta) e^{im\phi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if l<0 
  error('l must be >= 0')
end

if abs(m)>l
  error('m must be <= l')
end

norm = Nlm(l,m);
x = cos(theta);
temp = legendre(l,x);
if m>=0 
  Plm = temp(m+1);
else
  Plm = (-1)^abs(m) * (factorial(l-abs(m))/factorial(l+abs(m))) * ...
        temp(abs(m)+1);
end
y = norm * Plm * exp(i*m*phi);

return
