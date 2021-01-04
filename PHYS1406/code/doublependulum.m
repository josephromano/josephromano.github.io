%function  doublependulum(thetai, thetadoti, r)
%
% Inputs:
%
%   thetai: initial positions (degrees)
%   thetadoti: initial velocities (degrees per second)
%   ratio = m1/m2 (m1 assumed to be >> m2)
%
% small oscillations of coplanar double pendulum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample inputs
thetai = [5 -5];
thetadoti = [0 0];
r = 50; 

close all

% assign certain parameter values
l = 1; % pendula length (assumed equal)
g = 1;  % acceleration due to gravity
m2 = 1;
m1 = r*m2;
 
% convert degrees to radians
thetai = thetai* pi/180;
thetadoti = thetadoti* pi/180;

% convert ratio to small parameter epsilon
epsilon = sqrt(1/r);

% eigenfrequencies
w1 = sqrt(g/l)*(1+0.5*epsilon);
w2 = sqrt(g/l)*(1-0.5*epsilon);
wbeat = sqrt(g/l)*epsilon;

% eigenvectors
norm = 1/(l*epsilon*sqrt(2*m1));
a1 = norm * transpose([-epsilon 1]);
a2 = norm * transpose([ epsilon 1]);

% matrices
A = [a1 a2];
B = [w1*a1 w2*a2];

% initial displacements and velocities (from equilibrium position)
thetai = transpose(thetai);
thetadoti = transpose(thetadoti);

% solve for complex coefficients C using initial conditions
ReC = inv(A)*thetai;
ImC = inv(B)*thetadoti;
C = ReC + i * ImC;

% discrete times
% choose Tmax = 2 periods of beat freq = sqrt(g/l)*epsilon
numT = 500;
Tmax = 2*(2*pi/wbeat);
t = linspace(0, Tmax, numT);

% normal coords
zeta1 = C(1)*exp(-i*w1*t);
zeta2 = C(2)*exp(-i*w2*t);

% angles (need to take real part)
theta1 = real(A(1,1)*zeta1 + A(1,2)*zeta2);
theta2 = real(A(2,1)*zeta1 + A(2,2)*zeta2);

% convert to x-y coordinates
x1 =  l*sin(theta1);
y1 = -l*cos(theta1);
x2 =  x1 + l*sin(theta2);
y2 =  y1 - l*cos(theta2);

figure(1)
% make movie
for ii=1:numT

  pause(.01)

  plot(x1(ii), y1(ii), 'bo', x2(ii), y2(ii), 'bo');
  line([0 0], [0 -2.5], 'Color', 'g');
  line([0 x1(ii)], [0 y1(ii)], 'Color', 'k')
  line([x1(ii) x2(ii)], [y1(ii) y2(ii)], 'Color', 'k')
  axis equal
  ylim([-2.5 0])
  xlim([-0.75 0.75]);
  F(ii) = getframe;

end

% plot angles
figure(2)
plot(t, theta1*180/pi, 'b', t, theta2*180/pi, 'r');
grid on

fname = ['doublependulum.avi'];
movie2avi(F,fname)

