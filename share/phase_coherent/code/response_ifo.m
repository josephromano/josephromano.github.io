function [RG, RC] = response_ifo(l, m, f, u, v, x0)
%
% Inputs:
% l, m  - labels mode (l>=2, |m|<=l)
% f     - array of frequencies (Hz) 
% u, v  - unit vectors along the arms of the interferometer
% x0    - position vector of detector
%
% NOTE: u,v,x0 can be arrays u(1,:), u(2,:), u(3,:), etc
% corresponding to detector orientation at different times t
%
% Outputs:
% RG, RC - grad, curl response for given lm, array of frequencies,
%          array of u,v,x0 (Nt x Nf matrix)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
const = physConstants('mks');

% extract number of discrete time and discrete freqs (could = 1)
Nt = size(u,2);
Nf = length(f);

% extract angles for u, v
thetau = acos(u(3,:));
phiu = atan2(u(2,:), u(1,:));
thetav = acos(v(3,:));
phiv = atan2(v(2,:), v(1,:));

% extract unit vector in x0 direction
absX0 = sqrt(x0(1,:).^2 + x0(2,:).^2 + x0(3,:).^2);
tol = 1e-9;
if abs(absX0) < tol
  % doesn't matter in which direction unit vector points
  u0 = [zeros(1,Nt); zeros(1,Nt); ones(1,Nt)];
else
  u0 = [x0(1,:)./absX0; x0(2,:)./absX0; x0(3,:)./absX0];
end

% extract angles for u0  
theta0 = acos(u0(3,:));
phi0 = atan2(u0(2,:), u0(1,:));

% calculate alpha values (Nt x Nf array)
absX0 = reshape(absX0,Nt,1);
f = reshape(f,1,Nf);
alpha = 2*pi*absX0*f/const.c;

% calculate RG, RC
RG = zeros(size(alpha));
RC = zeros(size(alpha));

for L=l-2:l+2

  jL = sphericalbessel(L, alpha);
  fac = 4*pi*(-sqrt(-1))^L * jL * ...
        sqrt((2*2+1)*(2*l+1)*(2*L+1)/(4*pi)) * ...
        w3j(2, 2, l, -2, L, 0);

  for mp=-2:2

    % apply selection rule to eliminate sum_M=-L to L
    % but make sure that -L <= M <= M
    M = mp - m;

    if abs(M)<=L
      Y0 = sphericalharmonic(L, M, theta0, phi0);
      Y0 = transpose(Y0);
      Y0star = conj(Y0);

      % calculate Rbar^G_2m
      Yu = sphericalharmonic(2, mp, thetau, phiu);
      Yu = transpose(Yu);
      Yv = sphericalharmonic(2, mp, thetav, phiv);
      Yv = transpose(Yv);
      Rbar = (4*pi/5)*sqrt(1/3)*(Yu - Yv);
   
      temp = Rbar .* fac .* Y0star * w3j(2, -mp, l, m, L, M);

      RG = RG + temp * (-1)^mp * (1/2) * ((-1)^(l+L) + 1);
      RC = RC + temp * (-1)^mp * (1/(2*sqrt(-1))) * ((-1)^(l+L) - 1);
    end

  end
end

return

