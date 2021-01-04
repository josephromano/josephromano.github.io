function recordsound(reverse)
%
% record sound and fourier analyse
% (if reverse=1, reverse sound in time-domain)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

try, reverse; catch reverse=0; end

% define audio recorder
Fs = 22050; % sample rate (Hz)
r = audiorecorder(Fs, 16, 1);

% begin recording
cont = input('hit any key to start recording');
record(r);

% stop recording
cont = input('hit any key to stop recording');
stop(r);
%pause(r);
%resume(r);     

% playback and plot time series and power spectrum
cont = input('hit any key to play recording and display time series');

%p = play(r);
z = getaudiodata(r, 'int16'); % get data as int16 array
y = double(z); % convert to double floating point
y = y/(max(abs(y))); % renormalize so that i can use audioplayer to play it later
 
% reverse sound in time-domain if desired
if reverse==1
  y = flipud(y);
end

% play sound
p = audioplayer(y, Fs);
play(p);

N = length(y);
tmin = 0;
tmax = N/Fs;
t = linspace(tmin, tmax, N);

figure
subplot(2,1,1);
plot(t, y)
xlabel('time (sec)');
ylabel('y(t)');

[f, P] = fourieranalyze(t, y, Fs);
subplot(2,1,2);
loglog(f, P);
xlabel('freq (Hz)');
ylabel('power');
  
% zoom in
while 1

  % get input from mouse
  [u,v]=ginput(2);
  P = [u(1);v(1)];
  Q = [u(2);v(2)];

  if tmin <= min(P(1),Q(1)) & max(P(1),Q(1)) <= tmax

    % redefine grid of complex plane
    tmin = min(P(1),Q(1));
    tmax = max(P(1),Q(1));

    %ymin = min(P(2),Q(2));
    %ymax = max(P(2),Q(2));

    ndx1 = find(tmin <= t);
    ndx2 = find(t <= tmax);
    ndx = intersect(ndx1, ndx2);

    tnew = t(ndx);
    ynew = y(ndx);

  else

    return
    
    % revert to original time series
    %tmin = t(1);
    %tmax = t(end);
    %tnew = t;
    %ynew = y;

  end

  figure
  subplot(2,1,1);
  plot(tnew, ynew)
  xlabel('time (sec)');
  ylabel('y(t)');

  [f, P] = fourieranalyze(tnew, ynew, Fs);
  subplot(2,1,2);
  loglog(f, P);
  xlabel('freq (Hz)');
  ylabel('power');

end

