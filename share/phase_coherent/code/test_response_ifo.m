function test_response_ifo(l, m)
%
% test function for response_ifo.m
% l >= 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

thetau = pi/4;
phiu = 0;
thetav = pi/3;
phiv = pi/2;
theta0 = pi/4;
phi0 = pi/4;
absX0 = 1;

u = [sin(thetau)*cos(phiu); sin(thetau)*sin(phiu); cos(thetau)];
v = [sin(thetav)*cos(phiv); sin(thetav)*sin(phiv); cos(thetav)];
x0 = absX0*[sin(theta0)*cos(phi0); sin(theta0)*sin(phi0); cos(theta0)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test 1
f = 0;
[RG, RC] = response_ifo(l, m, f, u, v, x0);
      
% compare with static detector at origin
[RG_static, RC_static] = static_response_ifo(l, m, u, v);

RG
RG_static
RC
RC_static

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test 2
Nf = 200;
alpha = linspace(0, 2*pi*4, Nf);
f = alpha/(2*pi*absX0/(3e8));
[RG, RC] = response_ifo(l, m, f, u, v, x0);

% make plots
figure
subplot(2,1,1)
plot(alpha, real(RG), alpha, imag(RG)),
legend('real', 'imag')
xlabel('alpha')
ylabel('RG')
titlestr = ['grad response: l=' num2str(l) ', m=' num2str(m)];
title(titlestr)

subplot(2,1,2)
plot(alpha, real(RC), alpha, imag(RC)),
legend('real', 'imag')
xlabel('alpha')
ylabel('RC')
titlestr = ['curl response: l=' num2str(l) ', m=' num2str(m)];
title(titlestr)

% print to file
filename = ['response_ifo_l_' num2str(l) '_m_' num2str(m) '.jpg'];
print('-djpeg', filename);

return

