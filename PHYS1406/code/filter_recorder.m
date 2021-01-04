function filter_recorder
%
% record sound, filter, and calculate spectrogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% create figure and resize
f1=figure();
set(f1,'Position',[0 50 750 630])

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
%y = y/(max(abs(y))); % renormalize so that i can use audioplayer to play it later
 
while 1

  type = input('input filter type: h, l, or b: ', 's');

  switch type

    case 'h'
      fc = input('input high pass corner frequency in Hz: ');
      ynew = hipass(y, Fs, fc, 6, 1, 'butter');
      ynew = ynew/(max(abs(ynew)));

    case 'l'
      fc = input('input low pass cutoff frequency in Hz: ');
      ynew = lowpass(y, Fs, fc, 6, 1, 'butter');
      ynew = ynew/(max(abs(ynew)));

    case 'b'
      fl = input('input bandpass low  frequency in Hz: ');
      fh = input('input bandpass high frequency in Hz: ');
      ynew = bandpass(y, Fs, fl, fh, 6, 1, 'butter');
      ynew = ynew/(max(abs(ynew)));

    otherwise
      error('unrecognized filter type')
  end

  % play sound
  p = audioplayer(ynew, Fs);
  play(p);

  N = length(ynew);
  tmin = 0;
  tmax = N/Fs;
  t = linspace(tmin, tmax, N);

  subplot(2,2,1);
  plot(t, y/max(abs(y)))
  grid on
  xlabel('time (sec)');
  ylabel('y(t)');
  title('original')

  % calculate spectrogram and make plot
  Nfft = floor(max(256,floor(length(y)/256)));
  Noverlap = floor(0.5*Nfft);
  window = hamming(Nfft);
  subplot(2,2,2)
  spectrogram(y, window, Noverlap, Nfft, Fs, 'yaxis');
  ylim([0 8000])
  title('original')
  
  subplot(2,2,3);
  plot(t, ynew, 'r')
  grid on
  xlabel('time (sec)');
  ylabel('y(t)');
  title('filtered')

  % calculate spectrogram and make plot
  Nfft = floor(max(256,floor(length(ynew)/256)));
  Noverlap = floor(0.5*Nfft);
  window = hamming(Nfft);
  subplot(2,2,4)
  spectrogram(ynew, window, Noverlap, Nfft, Fs, 'yaxis');
  ylim([0 8000])
  title('filtered')
  
end

return
