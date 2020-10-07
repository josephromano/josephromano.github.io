function rd = detectorOrbit(t, detector)
%
% DETECTORORBIT - calculates the position vector of the detector
% relative to the SSB in equatorial coordinates taking into account
% earth's yearly orbital and daily rotational motions.  (assumes
% circular orbit of r=1 a.u. and that the sun is at the SSB.)
%
% rd = detectorOrbit(t)
% 
% t        - array of times (sec)
% detector - structure containing information about the detector
%            field      content
%            site       string identifier (e.g., 'H1')
%            r          position vector relative to center of earth
%            u          unit vector along arm 1 (geographic coords)
%            v          unit vector along arm 2 (geographic coords)
%
% rd       - position vector of detector [3, numel(t)] relative
%            to SSB in equatorial coordinates
%
% $Id:$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const = physConstants('mks');
wE = 2*pi/const.sidDay; % angular freq (rad/sec)

% extract position vector of detector relative to center of earth
% in earth-fixed geographic coordinates
r = detector.r;

% rotate position vector of detector to take into account earth's 
% daily rotational motion 
% NOTE: active rotation by wE*t = passive rotation by -wE*t
rdE(1,:) =  r(1)*cos(-wE*t) + r(2)*sin(-wE*t);
rdE(2,:) = -r(1)*sin(-wE*t) + r(2)*cos(-wE*t);
rdE(3,:) =  r(3)*ones(size(t));

% calculate position vector of center of earth relative to SSB 
% in equatorial coordinates due to earth's yearly orbital motion 
rES = earthOrbit(t);

% calculate position vector of detector relative to SSB
rd = rdE + rES;

%%
return
