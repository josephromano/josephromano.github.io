%
% script for testing detector orbit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = physConstants('mks');

% discrete times
t = [0:.0001:1]*const.sidDay;

% get detector information
site = 'default';
detector = getdetectorNew(site);

% calculate orbit of center of earth
rE = earthOrbit(t);
 
% calculate orbit of detector
rd = detectorOrbit(t, detector);

% make plots
figure(1);
plot3(rE(1,:), rE(2,:), rE(3,:), 'b');
hold on
plot3(rd(1,:), rd(2,:), rd(3,:), 'r');

return
