function recordsoundandsave(filename)
%
% record sound and write to .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define audio recorder
Fs = 22050; % sample rate (Hz)
r = audiorecorder(Fs, 16, 1);

% begin recording
cont = input('hit any key to start recording');
record(r);

% stop recording
cont = input('hit any key to stop recording');
stop(r);

% extract sound data array and convert to double with proper normalization
z = getaudiodata(r, 'int16'); % get data as int16 array
y = double(z); % convert to double floating point
y = y/(max(abs(y))); % renormalize so that i can use audioplayer to play it later

% construct time series
N = length(y);
tmin = 0;
tmax = N/Fs;
t = linspace(tmin, tmax, N);

% save to .mat file
save(filename)

% save to .wav file
audiowrite([filename '.wav'], y, Fs)

return

