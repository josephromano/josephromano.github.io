function rigidbody(rbtype)
%
% rbtype = 'sphere', 'football', or 'frisbee"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% ellipsoid parameters
a=1; b=1;

switch rbtype
  case 'sphere'
    c=1;
    fname = 'sphere';

  case 'frisbee'
    c=0.2;
    fname = 'frisbee';

  case 'football'
    c=2;
    fname = 'football';

  otherwise
    error('unknown rigidbody type');
end

% plot limits
plotlim = max([a b c]);
xmax=plotlim; ymax=plotlim; zmax=plotlim;

% fix orientation of spin axis theta0 and angular momentum L
theta0 = pi/2-pi/20; % approx horizontal
theta0 = pi/20; % approx vertical
L = 1; % along z-axis

% moments of inertia
M=1; 
I1 = (1/5)*M*(b^2 + c^2);
I2 = (1/5)*M*(c^2 + a^2);
I3 = (1/5)*M*(a^2 + b^2);

% relation between theta and alpha
% (theta = angle between n3 and \vec L) 
% (alpha = angle between n3 and \vec\omega)
alpha = atan(tan(theta0)*I3/I1);

% other angular frequencies 
omega3 = L*cos(theta0)/I3;
omega = omega3/cos(alpha);
Omega = omega3*(I3-I1)/I1;
phidot = (I3*omega3)/(I1*cos(theta0));
psidot = -Omega;

% display some results
fprintf('I1 = I2 = %f\n', I1)
fprintf('I3 = %f\n', I3)
fprintf('alpha = %f, theta = %f degrees\n', alpha*180/pi, theta0*180/pi);
%fprintf('omega :phidot = %f\n', I1*cos(theta0)/(I3*cos(alpha)));
fprintf('omega3:phidot = %f\n', I1*cos(theta0)/I3);
%fprintf('psidot:phidot = %f\n', I1*cos(theta0)/I3 - cos(theta0));

% angular coordinates on surface of ellipsoid
u = transpose(linspace(0, pi, 30));
v = transpose(linspace(0, 2*pi, 30));

% construct unit vector in direction of (u,v)
nx = sin(u) * cos(transpose(v));
ny = sin(u) * sin(transpose(v));
nz = cos(u) * ones(1,length(v));

% standard ellipse 
x = a*nx;
y = b*ny;
z = c*nz;

% discrete times (2 wobble periods)
Tmax = 2 * 2*pi/phidot;
t = transpose(linspace(0, Tmax, 100));

% euler angle time dependence
phi0 = 0;
psi0 = 0;
phi = phi0 + phidot*t;
psi = psi0 + psidot*t;
theta = theta0*ones(length(t),1);

% open VideoWriter for .mp4 file 
vw = VideoWriter (fname, 'MPEG-4');
open(vw);

for ii=1:length(t)
  % euler-angle rotation matrix (CW passive rotation = CCW active rotation)
  R1 = [cos(-phi(ii)) sin(-phi(ii)) 0; ...
       -sin(-phi(ii)) cos(-phi(ii)) 0; ...
        0 0 1];
  R2 = [cos(-theta(ii)) 0 -sin(-theta(ii)); ...
        0 1 0;
        sin(-theta(ii)) 0 cos(-theta(ii))];
  R3 = [cos(-psi(ii)) sin(-psi(ii)) 0; ...
       -sin(-psi(ii)) cos(-psi(ii)) 0; ...
        0 0 1];
  R = R1 * R2 * R3; % note opposite order since CW passive rotation
  
  xp = R(1,1)*x + R(1,2)*y + R(1,3)*z;
  yp = R(2,1)*x + R(2,2)*y + R(2,3)*z;
  zp = R(3,1)*x + R(3,2)*y + R(3,3)*z;

  % make plots
  figure(1)
  cm = ones(size(x));
  surf(xp,yp,zp,gradient(y))
  axis equal
  axis([-xmax xmax -ymax ymax -zmax zmax])
  xlabel('X-axis')
  ylabel('Y-axis')
  zlabel('Z-axis')
  camlight(120, 60) 
  lighting phong
  view(120, 20)
  pause(.01)

  frame = getframe;
  writeVideo(vw, frame);
  

end

close(vw)

return
