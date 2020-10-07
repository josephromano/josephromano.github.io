function n = simulateNoise(Cn, seed)
%
% simulate a column vector of complex noise values given a noise
% covariance matrix (assume zero mean multivariate gaussian with
% statistically independent real and imag parts)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
  seed;
  randn('state',seed);
catch
  % do nothing
end

% simulate data a multivariate gaussian with zero mean
x = mvnrnd(zeros(1,size(Cn,1)), Cn);
y = mvnrnd(zeros(1,size(Cn,1)), Cn);
n = x + sqrt(-1)*y;

n = transpose(n);

return

