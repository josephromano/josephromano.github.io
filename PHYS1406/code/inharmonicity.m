function [f1, I] = inharmonicity(f0, n)
%
% calculates the deviation of the nth partial from the
% exact nth harmonic for a real piano string
%
% inputs:
%   f0 - ideal fundamental frequency 
%        (e.g., f0 = 261.63 Hz for equal-tempered C4)
%   n  - partial number (1, 2, ...)
%
% outputs:
%   f1 - fundamental frequency for real string
%   I  - interval in cents for fn/f1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign common parameter values
E = 2e11; % N/m^2 for steel
rho = 7840; % kg/m^3 for steel
r = 5e-4; % m (radius)
mu = rho*pi*r^2; % kg/m
tau = 150*4.45; % N (tension)

% calculate B 
B = pi^3 * E * r^4 * mu * f0^2/(tau^2)

% fundamental frequency for real piano string
f1 = f0*(1+0.5*B);

% calculate frequency ratio f_n/(n f_1) from the formula
% f_n = n f_1(1 + 0.5(n^2-1)B)
ratio = 1 + 0.5 * (n^2-1) * B;

% convert frequency ratio to an interval in cents
I = ratio2cents(ratio);

fprintf('fundamental freq (ideal string): f0 = %f Hz\n', f0);
fprintf('fundamental freq (real  string): f1 = %f Hz\n', f1);
fprintf('freq ratio fn/(n f1) for partial n=%d: I = %f cents\n', n, I)

return

