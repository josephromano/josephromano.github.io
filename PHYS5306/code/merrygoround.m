function merrygoround(x0, y0, v0, theta0, T)

% motion of ball on merry-go-round wrt merry-go-round frame
%
% x0, y0      - initial position measured from center of merry-go-round
% v0, theta0  - intial velocity and direction wrt merry-go-round frame
% T           - duration of motion
%
% e.g., from Marion & Thornton p.389: 
% merrygoround(-0.5, 0, 1.5, pi/2, 0.86)
% merrygoround(-0.5, 0, 0.8, pi/2, 2.9)
% merrygoround(-0.5, 0, 0.45, pi/2, 17.3)
% merrygoround(-0.5, 0, 0.328, pi/2, 5)
% merrygoround(-0.5, 0, 0.47, pi/4, 3.83)
% merrygoround(-0.5, 0, 0.283, pi/4, 3.3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% parameters
R = 1; % radius of merrry-go-round (m/s)
w = 1; % angular velocity of merry-go-round (rad/s)

% discrete times
t = linspace(0,T,100);

% initial velocity wrt inertial frame 
% (1st terms wrt mgr frame, 2nd terms are from w x r)
v0x = v0*cos(theta0) - w*y0;
v0y = v0*sin(theta0) + w*x0;

% motion in inertial frame
x = x0 + v0x*t;
y = y0 + v0y*t;

% motion in merry-go-round frame
xp =  cos(w*t).*x + sin(w*t).*y;
yp = -sin(w*t).*x + cos(w*t).*y;

% make plot
figure
circle(1,[0 0]);
hold on
plot(0,0,'+k');
plot(xp, yp, 'linewidth', 2)
xlim([-1 1]);
ylim([-1 1]);
axis square
title(['v0 = ' num2str(v0) '; theta0 = ' num2str(theta0) '; T = ' num2str(T)]);
print -depsc2 merrygoround.eps

return
