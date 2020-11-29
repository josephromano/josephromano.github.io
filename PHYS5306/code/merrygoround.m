function merrygoround(x0, y0, v0, theta0, T)

% motion of ball on merry-go-round wrt merry-go-round frame
%
% x0, y0      - initial position measured from center of merry-go-round
% v0, theta0  - intial velocity and direction wrt merry-go-round frame
% T           - duration of motion
%
% e.g., from Marion & Thornton p.389: 
% merrygoround(-0.5, 0, 1.5, pi/2, 0.86)
% merrygoround(-0.5, 0, 0.47, pi/4, 3.83)
% merrygoround(-0.5, 0, 0.8, pi/2, 2.9)
% merrygoround(-0.5, 0, 0.45, pi/2, 17.3)
% merrygoround(-0.5, 0, 0.328, pi/2, 5)
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

% motion of line of sight in rotating frame
xlos = -sin(w*t);
ylos = cos(w*t);

% Prepare the new file.
filename = ['merrygoround_' num2str(v0) '_' num2str(theta0) '_' num2str(T) '.avi'];
vidObj = VideoWriter(filename);

% Set and view the frame rate (frames/sec)
vidObj.FrameRate = 25;
open(vidObj);

% animation
for ii=1:length(t)

  % make plots
  figure(1)
  plot(x(ii), y(ii), 'ko', 0, 0, 'k+', 'MarkerSize', 12);
  circle(1,[0 0]);
  line([0 xlos(ii)], [0 ylos(ii)])
  axis equal
  axis([-1 1 -1 1])
  xlabel('X-axis (inertial)')
  ylabel('Y-axis (inertial)')
  pause(.01)
  hold off
  currFrame = getframe(gca);
  writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);

% make plot in rotating frame
figure(2)
circle(1,[0 0]);
hold on
line([0 xlos(1)], [0 ylos(1)])
plot(0,0,'+k');
plot(xp, yp, 'linewidth', 2);
xlim([-1 1]);
ylim([-1 1]);
xlabel('x-axis (rotating)')
ylabel('y-axis (rotating)')

axis square
title(['v0 = ' num2str(v0) '; theta0 = ' num2str(theta0) '; T = ' num2str(T)]);
filename = ['merrygoround_' num2str(v0) '_' num2str(theta0) '_' num2str(T) '.eps'];
print('-depsc2', filename)

return
