function [f, P] = fourieranalyze(t, y, fs)
%
% fourier analyze a time series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seglength = floor(length(t)/4);
H=spectrum.welch('Hann', seglength);
%H=spectrum.periodogram;
hpsd = psd(H, y, 'Fs', fs);
f = hpsd.Frequencies;
P = hpsd.Data;

return
