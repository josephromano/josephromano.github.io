function mu = match(a,b)
%
% calculate the match (coherence) between two complex vectors
%
% NOTE: assume vectors are column vectors (Nx1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 0.5*(a'*b + b'*a)/(sqrt(a'*a) * sqrt(b'*b));

return

