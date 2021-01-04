function harmonics(f0)
%
% calculate nearest equal-tempered frequencies for
% first 8 harmonics of f0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try, f0; catch, f0=110; end; % A2 = 110 Hz

N=8;
semitones = [0 12 19 24 28 31 34 36];

for n=1:N
  f(n) = n*f0;
  fET(n) = f0 * 2^(semitones(n)/12);
  c(n) = ratio2cents(f(n)/fET(n));
  fprintf('n=%d, f=%.2f, fET=%.2f, cents=%d\n', n, f(n), fET(n), round(c(n)));
end

return

