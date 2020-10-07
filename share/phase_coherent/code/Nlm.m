function norm = Nlm(l,m)
%
% calculate normalization factor for spherical harmonics Ylm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if l<0 
  error('l must be >= 0')
end

if abs(m)>l
  error('m must be <= l')
end

norm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );

return

