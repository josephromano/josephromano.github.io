function mat2wav(matfile, reverse)
%
% convert .mat file to .wav file
%
% (if reverse=1, reverse sound in time domain)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

try, reverse; catch reverse=0; end

% add .mat extension
filename_mat = [matfile '.mat'];

% add .wav
filename_wav = [matfile '.wav'];

% extract saved data (r, Fs)
load(filename_mat);
deltaT = t(2)-t(1);

% reverse sound in time-domain if desired
if reverse==1
  % add "-reverse.wav"
  filename_wav = [matfile '-reverse.wav'];
  y = flipud(y);
end

% play recorded sound
%p = audioplayer(y, Fs);
%play(p);
%pause(tmax)

% save sound to .wav file
audiowrite(filename_wav, y, Fs)

end

