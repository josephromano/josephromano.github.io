function  triatomic(qi, vi, r, k)
%
% Inputs:
%
%   qi: initial positions
%   vi: initial velocities
%   r = M/m (mass ratio)
%   k: spring constant
%
% linear triatomic molecule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample inputs
%qi = [-1.5 0 1.5];
%qi = [-0.5 0 1.5];
%qi = [-1 0.5 1];
%qi = [-1 0.5 1];
%vi = [0 0 0];

%qi = [-1 0 1];
%vi = [-0.1 0 0.1];
%vi = [0 0.5 0];

%r = 1; 
%k= 1;

% assign certain parameter values
m = 1;
M = m*r;

% equilibrium positions
q1_0 = -1;
q2_0 = 0; 
q3_0 = 1;

% eigenfrequencies
w1 = 0;
w2 = sqrt(k/m);
w3 = sqrt((k/m)*(1+2*m/M));

% eigenvectors
N1 = 1/sqrt(2*m + M);
a1 = N1 * transpose([1 1 1]);
N2 = 1/sqrt(2*m);
a2 = N2 * transpose([1 0 -1]);
N3 = 1/sqrt(2*m*(1 + 2*m/M));
a3 = N3 * transpose([-1 2*m/M -1]);

% matrices
A = [a1 a2 a3];
B = [w1*a1 w2*a2 w3*a3];

% initial displacements and velocities (from equilibrium position)
etai = transpose(qi) - transpose([q1_0 q2_0 q3_0]); 
etadoti = transpose(vi);

% solve for complex coefficients C using initial conditions
ReC = inv(A)*etai;
% B is singular so ImC = inv(B)*etadoti won't work 
% restrict to submatrix instead
temp = inv(B(2:3,2:3))*etadoti(2:3);
ImC = [0; temp];
C = ReC + i * ImC;

% discrete times
% choose Tmax = 4 periods of w3 (this is arbitrary)
numT = 100;
Tmax = 4 * 2*pi/w3; 
t = linspace(0, Tmax, numT);

% normal coords
zeta1 = C(1)*exp(-i*w1*t);
zeta2 = C(2)*exp(-i*w2*t);
zeta3 = C(3)*exp(-i*w3*t);

% displacements from equilibrium
eta1 = A(1,1)*zeta1 + A(1,2)*zeta2 +A(1,3)*zeta3;
eta2 = A(2,1)*zeta1 + A(2,2)*zeta2 +A(2,3)*zeta3;
eta3 = A(3,1)*zeta1 + A(3,2)*zeta2 +A(3,3)*zeta3;

% positions (need to take real part)
q1 = q1_0 + real(eta1);
q2 = q2_0 + real(eta2);
q3 = q3_0 + real(eta3);

% make movie
for ii=1:numT

  pause(.01)

  plot(q1(ii), 0, 'ro', q2(ii), 0, 'bo', q3(ii), 0, 'go');
  axis equal
  ylim([-1 1])
  xlim([-2 2]);
  F(ii) = getframe;

end

fname = ['triatomic.avi'];
movie2avi(F,fname)

