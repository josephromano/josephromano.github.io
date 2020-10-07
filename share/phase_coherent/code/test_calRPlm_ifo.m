% test calRPlm_ifo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% set some parameters
lmax = 3;
Nt = 1;
f = 100;

const = physConstants('mks');

% calculate number of modes Nm 
Nm = (lmax+1)^2 - 4;

% angular frequency of Earth's daily rotation (rad/s)
wE = 2*pi/(const.sidDay);

% discrete times
t = [0:1/Nt:1-1/Nt]*const.sidDay;

% get detector information
site = 'H1';
detector = getdetectorNew(site);

% extract detector parameters
r = detector.r;
u = detector.u;
v = detector.v;

% rotate detector unit vectors u,v keeping source fixed (equatorial coords)
% (rotation by -wE*t around z-axis)
u_t = zeros(3,Nt);
u_t(1,:) =  cos(-wE*t)*u(1) + sin(-wE*t)*u(2);
u_t(2,:) = -sin(-wE*t)*u(1) + cos(-wE*t)*u(2);
u_t(3,:) =  u(3);
v_t = zeros(3,Nt);
v_t(1,:) =  cos(-wE*t)*v(1) + sin(-wE*t)*v(2);
v_t(2,:) = -sin(-wE*t)*v(1) + cos(-wE*t)*v(2);
v_t(3,:) =  v(3);

% rotate position vector detector from center of  earth
% (rotation by -wE*t around z-axis)
r_t = zeros(3,Nt);
r_t(1,:) =  cos(-wE*t)*r(1) + sin(-wE*t)*r(2);
r_t(2,:) = -sin(-wE*t)*r(1) + cos(-wE*t)*r(2);
r_t(3,:) =  r(3);

% calculate integrated repsonse for each ifo
R_Plm = calRPlm_ifo(f, u_t, v_t, r_t, lmax);

% debugging
R_Plm.mvec'
R_Plm.lvec'
R_Plm.data'
size(R_Plm.data)

return

