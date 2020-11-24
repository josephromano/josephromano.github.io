% ellipsoid parameters
a=2;
b=1;
plotlim = max([a b]);
xmax=plotlim; ymax=plotlim;

% spin axis
theta0 = pi/6;

% angular coordinates on surface of ellipsoid
v = transpose(linspace(0, 2*pi, 30));

% construct unit vector in direction of (u,v)
nx = cos(transpose(v));
ny = sin(transpose(v));

% standard ellipse 
x = a*nx;
y = b*ny;

psi(1)=pi/6;

for ii=1:1

  % euler-angle rotation matrix (CW passive rotation = CCW active rotation)
  R = [cos(-psi(ii)) sin(-psi(ii)); ...
      -sin(-psi(ii)) cos(-psi(ii))];
  
  xp = R(1,1)*x + R(1,2)*y;
  yp = R(2,1)*x + R(2,2)*y;

  % make plots
  figure(1)
  plot(xp,yp)
  axis equal
  axis([-xmax xmax -ymax ymax])
  xlabel('X-axis')
  ylabel('Y-axis')
 
end

return
