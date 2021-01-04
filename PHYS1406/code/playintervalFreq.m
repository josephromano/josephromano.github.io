function playintervalFreq(f1, f2, dur)
%
% play two notes in unison (notes specified by frequency)
%
% Input:
%    f1, f2 - frequencies in Hz 
%    dur    - duration in seconds (default = 4 s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try, dur; catch dur=1; end
%playnoteFreq(f1, dur);
%playnoteFreq(f2, dur);

%return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try, dur; catch dur=4; end

% construct pure tones (sine waves)
tmin = 0;
tmax = dur; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
N = tmax/deltaT;
t = linspace(tmin, tmax, N);
y1 = sin(2*pi*f1*t);  
y2 = sin(2*pi*f2*t); 
y = 0.5*(y1+y2);
 
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
