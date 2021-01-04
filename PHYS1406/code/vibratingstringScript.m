%
% script for running vibratingstring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = input('enter of initial configuration: fundamental, second, third, fourth, sinepulse, trianglecenter, trianglequarter\n');
 
switch type
  case 'simple'
    y = [0 .1  .15 .1 .05 0];

  case 'fundamental'
    x = linspace(0,1,51);
    y = 0.1*sin(pi*x);

  case 'second'
    n = 2;
    x = linspace(0,1,51);
    y = 0.1*sin(n*pi*x);

  case 'third'
    n = 3;
    x = linspace(0,1,51);
    y = 0.1*sin(n*pi*x);

  case 'fourth'
    n = 4;
    x = linspace(0,1,51);
    y = 0.1*sin(n*pi*x);

  case 'sinepulse'
    n = 4;
    x = linspace(0,1,100);
    z = 0.1*sin(n*pi*x);
    y = zeros(length(x),1);
    y(1:floor(end/4)) = z(1:floor(end/4));

  case 'trianglecenter'
    x = linspace(0,1,51);
    y = 0.05+0.05*sawtooth(2*pi*x, 0.5);

  case 'trianglequarter'
    x = linspace(0,1,51);
    y = zeros(length(x), 1);
    m = 0.1/x(floor(end/4));
    y(1:floor(end/4)) = m * x(1:floor(end/4));
    m = -0.1/(1-x(floor(end/4)+1));
    y(floor(end/4)+1:end) = 0.1 + m * (x(floor(end/4)+1:end) - x(floor(end/4)+1));

  otherwise
    error('unknown type')

end

vibratingstring(y);

