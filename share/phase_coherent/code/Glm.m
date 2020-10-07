function [Gp,Gm] = Glm(l, m, x)
%
% calculate G^+_lm(x), G^-lm(x) of x=cos(theta)
% NOTE: l>=2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = cos(theta), y = sin(theta)
y = sqrt(1-x^2);

temp = legendre(l,x);
if m>=0 
  Plm = temp(m+1);
else
  Plm = (-1)^abs(m) * (factorial(l-abs(m))/factorial(l+abs(m))) * temp(abs(m)+1);
end

temp = legendre(l-1,x);
if (m>l-1) || (m<-(l-1))
  Plm1m =0;
else
  if m>=0 
    Plm1m = temp(m+1);
  else
    Plm1m = (-1)^abs(m) * (factorial(l-1-abs(m))/factorial(l-1+abs(m))) * temp(abs(m)+1);
  end
end

Gp = -((l-m^2)/(y^2) + 0.5*l*(l-1))*Plm + (l+m)*(x/(y^2))*Plm1m;
Gm = (m/(y^2))*((l-1)*x*Plm - (l+m)*Plm1m);

return

