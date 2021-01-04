function playsound(reverse)
%
% play sound (if reverse=1, reverse sound in time domain)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

try, reverse; catch reverse=0; end

type = input('enter type of sound wave:\n  sine, triangle, square, sawtooth, pulse, noise, beats\n');

%type = input('enter type of sound wave:\n  handel, sine, triangle, square, sawtooth, pulse, doublepulse, seebecksiren, noise, beats\n');

switch type
  case('handel')

    data = load('handel');
    y = data.y;
    Fs = data.Fs;

    deltaT = 1/Fs;
    N = length(y);
    tmin = 0;
    tmax = N*deltaT;
    t = linspace(tmin, tmax, N);
    y = double(y);

  case('sine')
    
    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y = sin(2*pi*f*t);  

  case('square')

    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y = 0.5*(1+square(2*pi*f*t));

  case('triangle')

    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y = sawtooth(2*pi*f*t, 0.5);

  case('sawtooth')

    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y = sawtooth(2*pi*f*t, 1);

  case('pulse')
    
    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    d = tmin : 1/f: tmax;
    %y = pulstran(t, d, 'tripuls', 0.2/f);
    y = pulstran(t, d, 'rectpuls', 0.2/f);

  case('doublepulse')
    
    f = 2*440; % A5
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    d = tmin : 1/f: tmax;
    y = pulstran(t, d, 'rectpuls', 0.2/f);

  case('seebecksiren')
    
    f = 440; % concert A4
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    d = tmin : 1/f: tmax;
    y1 = pulstran(t, d, 'rectpuls', 0.2/f);

    % time shift y1 by 0.3*period
    ndx = find(t>0.3/f);
    y2 = [y1(ndx(1):end) y1(1:ndx(1)-1)];

    y = y1+y2;

  case('noise')

    f = 440; % reference freq (for display)
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y = randn(N,1);

  case('beats')
    
    f1 = 440; % concert A4
    f2 = 442; 
    tmin = 0;
    tmax = 4; % sec
    Fs = 22050; % sampling rate
    deltaT = 1/Fs;
    N = tmax/deltaT;
    t = linspace(tmin, tmax, N);

    y1 = sin(2*pi*f1*t);  
    y2 = sin(2*pi*f2*t);  
    y = y1+y2;

  otherwise

    error('unknown sound type');
  
end

% reverse sound in time-domain if desired
if reverse==1
  y = flipud(y);
end

figure
subplot(2,1,1);
if strncmp(type,'beats',5)
  ndx = length(t);
else
  T = 1/f;
  t0 = 10*T;
  temp = find(t>=t0);
  ndx = temp(1);
end
plot(t(1:ndx), y(1:ndx))
xlabel('time (sec)');
ylabel('y(t)');
grid on

[f, P] = fourieranalyze(t, y, Fs);
subplot(2,1,2);
loglog(f, P);
xlabel('freq (Hz)');
ylabel('power');
grid on

% play sound
p = audioplayer(y, Fs);
play(p); 
pause(4);

return
