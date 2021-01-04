function playchord(note1, level1, temp1, note2, level2, temp2, ...
                   note3, level3, temp3)
%
% play three notes in unison 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur = 4; % duration in second

% convert notes to frequencies
f1 = note2freq(note1, level1, temp1);
f2 = note2freq(note2, level2, temp2);
f3 = note2freq(note3, level3, temp3);

% construct pure tones (sine waves)
tmin = 0;
tmax = dur; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
N = tmax/deltaT;
t = linspace(tmin, tmax, N);
y1 = sin(2*pi*f1*t);  
y2 = sin(2*pi*f2*t); 
y3 = sin(2*pi*f3*t); 
y = (1/3)*(y1+y2+y3);
 
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
