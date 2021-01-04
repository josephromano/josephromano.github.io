function sound_spectrogram
%
% record sound and calculate spectrogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% create figure and resize
f1=figure();
set(f1,'Position',[0 50 750 630])

% define audio recorder
Fs = 22050; % sample rate (Hz)
r = audiorecorder(Fs, 16, 1);

while 1

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
 
  % play sound
  p = audioplayer(y/max(abs(y)), Fs);
  play(p);

  N = length(y);
  tmin = 0;
  tmax = N/Fs;
  t = linspace(tmin, tmax, N);

  % plot time series
  subplot(2,1,1)
  plot(t, y/max(abs(y)))
  xlabel('time (sec)');
  ylabel('y(t)');

  % calculate spectrogram and make plot
  Nfft = floor(max(256,floor(length(y)/256)));
  Noverlap = floor(0.5*Nfft);
  window = hamming(Nfft);
  subplot(2,1,2)
  spectrogram(y, window, Noverlap, Nfft, Fs, 'yaxis');
  %ylim([0 8000])

  % repeat??
  repeat = input('repeat (1=y, 0=n): ');
  if repeat~=1
    return
  end

end

return
