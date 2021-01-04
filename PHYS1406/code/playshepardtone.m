function playshepardtone(note, dur);
%
% play notes an octave apart in equal temperament
% with intensities weighted by a gaussian centered at C4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gaussian intensity distribution for partials
% (peak at C4, sigma in octaves)
C4 = note2freq('C', 4, 'equal');
sigma = 1.5; 
maxlevel = 9; % from level 0 to 8

% calculate frequencies and intensities for notes an octave apart
for ii=1:maxlevel
  switch note
    case {'C', 'C#', 'D', 'D#', 'Eb', 'E', 'F', 'F#'}
      f(ii) = note2freq(note, ii-1, 'equal');
    case {'G', 'G#', 'Ab', 'A', 'A#', 'Bb', 'B'}
      f(ii) = note2freq(note, ii-2, 'equal');
  end

  ratio = C4/f(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensity(ii) = exp(-0.5*(diff/sigma)^2);
end

% construct pure tones (sine waves)
tmin = 0;
tmax = dur; % sec
Fs = 22050; % sampling rate
deltaT = 1/Fs;
N = floor(tmax/deltaT);
t = linspace(tmin, tmax, N);

y = zeros(maxlevel,N);
for ii=1:maxlevel
  A = sqrt(intensity(ii));
  y(ii,:) = A*sin(2*pi*f(ii)*t);  
end

% superimpose sine waves 
ytot = zeros(1,N);
for ii=1:maxlevel
  ytot = ytot + y(ii,:);
end
ytot = ytot/maxlevel;

% simulate attack and decay transients
Na = round(N*0.1);
Nd = round(N*0.2);
wa = sin((pi/2)*[0:Na-1]/Na);
wd = exp(-10*[0:Nd-1]/Nd);
ytot(1:Na) = ytot(1:Na).*wa;
ytot(end-(Nd-1):end) = ytot(end-(Nd-1):end).*wd;

% play sound
p = audioplayer(ytot, Fs);
play(p); 
pause(dur);

return
