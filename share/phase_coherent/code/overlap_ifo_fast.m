function [orf, f] = overlap_ifo_fast(site1, site2, flow, fhigh)
%
% calculate overlap reduction function for an isotropic uncorrelated
% stochastic background using simplified expression for orf(f)
%
% example:
%  overlap_ifo_fast('H1', 'L1', 0, 1000)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
const = physConstants('mks');

% get detector info
detector1 = getdetectorNew(site1);
detector2 = getdetectorNew(site2);
r1 = detector1.r;
u1 = detector1.u;
v1 = detector1.v;
r2 = detector2.r;
u2 = detector2.u;
v2 = detector2.v;

% calculate separation of ifo vertices
DeltaX = r2-r1;
absDeltaX = sqrt(sum(DeltaX.^2));

% extract angles of DeltaX
theta0 = acos(DeltaX(3)/absDeltaX);
phi0 = atan2(DeltaX(2),DeltaX(1));

% calculate rotation matrices 
% R1 = Rz(phi0)
R1 = [ cos(phi0), sin(phi0), 0; ...
      -sin(phi0), cos(phi0), 0; ...
          0,      0,         1];

% R2 = Ry(theta0) 
R2 = [  cos(theta0),  0, -sin(theta0); ...
        0,            1,   0; ...
        sin(theta0),  0,  cos(theta0)];

R = R2 * R1;

% rotate unit vectors along detector arms into computational frame
Ru1 = R * u1;
Rv1 = R * v1;
Ru2 = R * u2;
Rv2 = R * v2;

% discrete frequencies
Nf = 500;
f = linspace(flow, fhigh, Nf);

% discrete alpha values
alpha = 2*pi*f*absDeltaX/const.c;

% calculate spherical bessel functions
j0 = sphericalbessel(0,alpha);
j2 = sphericalbessel(2,alpha);
j4 = sphericalbessel(4,alpha);

%%%%%%%%%%%%%%
% calculate overlap function
orf = zeros(1, Nf);  
for m=-2:2

  % translated coordinates (works)
  %Rbar1 = static_response_ifo(2, m, u1, v1);
  %[R2, ignore] = response_ifo(2, m, f, u2, v2, DeltaX);
  %orf = orf + Rbar1.*conj(R2);

  % translated coordinates (works)
  %Rbar1 = static_response_ifo(2, m, u1, v1);
  %R2 = zeros(1, length(f));
  %for L=0:4

  %  jL = sphericalbessel(L, alpha);
  %  fac = 4*pi*(-sqrt(-1))^L * jL * ...
  %        sqrt((2*2+1)*(2*2+1)*(2*L+1)/(4*pi)) * ...
  %        w3j(2, 2, 2, -2, L, 0);

  %  for mp=-2:2
  %    % apply selection rule to eliminate sum_M=-L to L
  %    % but make sure that -L <= M <= M
  %    M = mp - m;

  %    if abs(M)<=L
  %      Y0 = sphericalharmonic(L, M, theta0, phi0);
  %      Y0star = conj(Y0);

  %      Rbarmp = static_response_ifo(2, mp, u2, v2);
  %      R2 = R2 + Rbarmp * fac * Y0star * w3j(2, -mp, 2, m, L, M)...
  %                * (-1)^mp * (1/2) * ((-1)^L + 1);
  %    end  

  %  end
  %end
  %orf = orf + Rbar1.*conj(R2);

  % translated and rotated coordinates (works)
  %Rbar1 = static_response_ifo(2, m, Ru1, Rv1);
  %[R2, ignore] = response_ifo(2, m, f, Ru2, Rv2, absDeltaX*[0; 0; 1]);
  %orf = orf + Rbar1.*conj(R2);

  % translated and rotated coordinates (works)
  Rbar1 = static_response_ifo(2, m, Ru1, Rv1);
  Rbar2 = static_response_ifo(2, m, Ru2, Rv2);
  R2 = Rbar2*(j0 ...
       - (-1)^m * (2*2+1)^2 * w3j(2,-m,2,m,2,0) * w3j(2,2,2,-2,2,0) * j2 ...
       + (-1)^m * (2*2+1)*(2*4+1) * w3j(2,-m,2,m,4,0) * w3j(2,2,2,-2,4,0) * j4);
               
  orf = orf + Rbar1.*conj(R2);

end

orf(1)

% make plot
semilogx(f,real(orf),'k-')
xlabel('freq (Hz)')
ylabel('overlap(f)')

outfile = ['overlap_' site1 '_' site2 '.eps'];
print('-depsc2', outfile);

return

