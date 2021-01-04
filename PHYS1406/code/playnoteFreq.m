function playnoteFreq(f, dur)
%
% play note (note specified by frequency)
%
% Input:
%    f  - frequency in Hz (default = 440 Hz)
%    dur   - duration in seconds (default = 1 s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try, f; catch f=440; end
try, dur; catch dur=1; end

% construct pure tone (sine wave)
tmin = 0;
tmax = dur; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
N = floor(tmax/deltaT);
t = linspace(tmin, tmax, N);
y = sin(2*pi*f*t);

% simulate attack and decay transients
Na = round(N*0.1);
Nd = round(N*0.2);
wa = sin((pi/2)*[0:Na-1]/Na);
wd = exp(-10*[0:Nd-1]/Nd);
y(1:Na) = y(1:Na).*wa;
y(end-(Nd-1):end) = y(end-(Nd-1):end).*wd;

% play sound
p = audioplayer(y, Fs);
play(p); 
pause(dur);

return

