function playrecordedsound(filename, reverse)
%
% play recorded sound that was saved to a file
% (if reverse=1, reverse sound in time domain)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

try, reverse; catch reverse=0; end

% add .mat extension
filename = [filename '.mat'];

% extract saved data (r, Fs)
load(filename);
deltaT = t(2)-t(1);

% reverse sound in time-domain if desired
if reverse==1
  y = flipud(y);
end

% play recorded sound
p = audioplayer(y, Fs);
play(p);
pause(tmax)

% plot time series
figure
subplot(2,1,1);
plot(t, y)
xlabel('time (sec)');
ylabel('y(t)');
grid on

[f, P] = fourieranalyze(t, y, Fs);
subplot(2,1,2);
loglog(f, P);
xlabel('freq (Hz)');
ylabel('power');
grid on
  
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

    % return
    return
    
    % revert to original time series
    %tmin = t(1);
    %tmax = t(end);
    %tnew = t;
    %ynew = y;

  end

  % play sound
  pnew = audioplayer(ynew, Fs);
  play(pnew);
  pause(deltaT*length(ynew));

  figure
  subplot(2,1,1);
  plot(tnew, ynew)
  xlabel('time (sec)');
  ylabel('y(t)');
  grid on

  [f, P] = fourieranalyze(tnew, ynew, Fs);
  subplot(2,1,2);
  loglog(f, P);
  xlabel('freq (Hz)');
  ylabel('power');
  grid on

end

