function overlap_ifo_orbit(f, T, deltaT)
%
% script for calculating overlap reduction function for an isotropic 
% uncorrelated stochastic background for an interferometer at different 
% times in earth's orbit around the sun.
%
% f = GW frequency (Hz) (e.g., 10 or 100)
% T = total time (sec)
% deltaT = time step (sec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% some constants
const = physConstants('mks');
R_au = const.au; % earth orbital radius (m)
T_yr = const.sidYr; % sidereal yr (sec)
T_day = const.sidDay; % sidereal day (sec)
norm = (1/2)*(8*pi/5); % normalization for colocated and coaligned detectors

% discrete times
N = floor(T/deltaT)+1;
t = deltaT*[0:1:N-1];

% interferomter arms always point in same direction
u = [0; 1; 0];
v = [0; 0; 1];
u1 = u;
v1 = v;
u2 = u;
v2 = v;

% calculate orf for ifo at different times in earth's orbit 
% around the sun
orf_orbit = zeros(N,1);
for ii=1:N

  if mod(ii,10)==0
    fprintf('working on %d of %d\n', ii, N)
  end

  % get detector info
  alpha = 2*pi*t(ii)/T_yr;
  r1 = [R_au; 0; 0];
  r2 = [R_au*cos(alpha); R_au*sin(alpha); 0];

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

  % discrete alpha values
  alpha = 2*pi*f*absDeltaX/const.c;

  % calculate spherical bessel functions
  j0 = sphericalbessel(0,alpha);
  j2 = sphericalbessel(2,alpha);
  j4 = sphericalbessel(4,alpha);

  %%%%%%%%%%%%%%
  % calculate overlap function
  orf = 0;
  for m=-2:2

    % translated and rotated coordinates (works)
    Rbar1 = static_response_ifo(2, m, Ru1, Rv1);
    Rbar2 = static_response_ifo(2, m, Ru2, Rv2);
    R2 = Rbar2*(j0 ...
         - (-1)^m * (2*2+1)^2 * w3j(2,-m,2,m,2,0) * w3j(2,2,2,-2,2,0) * j2 ...
         + (-1)^m * (2*2+1)*(2*4+1) * w3j(2,-m,2,m,4,0) * w3j(2,2,2,-2,4,0) * j4);
               
    orf = orf + Rbar1.*conj(R2);

  end

  % normalize
  orf = orf/norm;

  % check for colocated and coaligned detectors
  if ii==1
    orf_orbit(ii) = 1; 
  else
    orf_orbit(ii) = orf; 
  end

end

% caculate eigenvalues, eigenvectors of ORF matrix
HHdag = toeplitz(real(orf_orbit));
[V, D] = eig(HHdag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots

figure
plot(t/T_day, real(orf_orbit))
xlabel('time/1 day')
ylabel('overlap(t)')
outfile = ['overlap_ifo_orbit_f_' num2str(f) '_deltaT_' num2str(deltaT) '.eps'];
print('-depsc2', outfile);

figure
plot(diag(D))
ylabel('eigenvalues')
outfile = ['overlap_ifo_orbit_eigvalues_f_' num2str(f) '_deltaT_' num2str(deltaT) '.eps'];
print('-depsc2', outfile);

return

