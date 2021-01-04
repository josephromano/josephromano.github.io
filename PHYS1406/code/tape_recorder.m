function tape_recorder
%
% record sound and fourier analyse 
% (allow forward/reverse and speed adjustment)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% define audio recorder
Fs = 22050; % sample rate (Hz)
r = audiorecorder(Fs, 16, 1);

% begin recording
cont = input('hit any key to start recording');
record(r);

% stop recording
cont = input('hit any key to stop recording');
stop(r);

% get data
z = getaudiodata(r, 'int16'); % get data as int16 array
y = double(z); % convert to double floating point
y = y/(max(abs(y))); % renormalize so that i can use audioplayer to play it later
 
while 1

  direction = input('input forward (1) or reverse (-1): ');

  % reverse sound in time-domain if desired
  if direction == -1;
    ynew = flipud(y);
  else
    ynew = y;
  end

  % playback at different speed if desired
  speed = input('input speed factor (normal=1) : ');
  Fsnew = Fs*speed;

  % play sound
  p = audioplayer(ynew, Fsnew);
  play(p);

  N = length(ynew);
  tmin = 0;
  tmax = N/Fsnew;
  t = linspace(tmin, tmax, N);

  figure
  subplot(2,1,1);
  plot(t, ynew)
  grid on
  xlabel('time (sec)');
  ylabel('y(t)');

  [f, P] = fourieranalyze(t, ynew, Fsnew);
  subplot(2,1,2);
  loglog(f, P);
  grid on
  xlabel('freq (Hz)');
  ylabel('power');
  
end

return
