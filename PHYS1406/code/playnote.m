function playnote(note, level, temp, dur)
%
% play note (note specified by name)
%
% Input:
%    note  - e.g., 'C', 'C#', 'Eb', ...
%    level - e.g., 4 for C4
%    temp  - temperament e.g., equal, just, pyth(agorean), mean(tone) (default = equal)
%    dur   - duration in seconds (default = 1 s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try, note; catch note='C'; end
try, level; catch level=4; end
try, temp; catch temp='equal'; end
try, dur; catch dur=1; end

f = note2freq(note, level, temp);

% construct pure tone (sine wave)
tmin = 0;
tmax = dur; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
N = floor(tmax/deltaT);
t = linspace(tmin, tmax, N);
y = sin(2*pi*f*t);

% simulate attack and decay transients
if dur>=0.2
  % absolute attack and decay times
  Na = round(0.08/deltaT);
  Nd = round(0.10/deltaT);
else
  Na = round(N*0.1);
  Nd = round(N*0.2);
end

% apply transients
wa = sin((pi/2)*[0:Na-1]/Na);
wd = exp(-10*[0:Nd-1]/Nd);
y(1:Na) = y(1:Na).*wa;
y(end-(Nd-1):end) = y(end-(Nd-1):end).*wd;

% play sound
p = audioplayer(y, Fs);
play(p); 
pause(dur);

return

