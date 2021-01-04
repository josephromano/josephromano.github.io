% plot gaussian intensity distribution for shepard tones
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% gaussian intensity distribution for partials
% (peak at C4, sigma in octaves)
C4 = note2freq('C', 4, 'equal');
sigma = 1.5; % in octaves
maxlevel = 9; % levels 0 to 8

% smooth gaussian
Nf = 1e5;
freqs = linspace(10,1e4,Nf);
p = exp(-0.5*((log10(freqs)-log10(C4))/(sigma*log10(2))).^2);

figure
% calculate frequencies and intensities for notes an octave apart
for ii=1:maxlevel

  % shepard C
  fC(ii) = note2freq('C', ii-1, 'equal');
  ratio = C4/fC(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityC(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fC(ii)) log10(fC(ii))], [0, intensityC(ii)], 'color', 'k');
  
  % shepard Cs
  fCs(ii) = note2freq('C#', ii-1, 'equal');
  ratio = C4/fCs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityCs(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fCs(ii)) log10(fCs(ii))], [0, intensityCs(ii)], 'color', 'k');
 
  % shepard D
  fD(ii) = note2freq('D', ii-1, 'equal');
  ratio = C4/fD(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityD(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fD(ii)) log10(fD(ii))], [0, intensityD(ii)], 'color', 'b');

  % shepard D#
  fDs(ii) = note2freq('D#', ii-1, 'equal');
  ratio = C4/fDs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityDs(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fDs(ii)) log10(fDs(ii))], [0, intensityDs(ii)], 'color', 'r');

  % shepard E
  fE(ii) = note2freq('E', ii-1, 'equal');
  ratio = C4/fE(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityE(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fE(ii)) log10(fE(ii))], [0, intensityE(ii)], 'color', 'g');

  % shepard F
  fF(ii) = note2freq('F', ii-1, 'equal');
  ratio = C4/fF(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityF(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fF(ii)) log10(fF(ii))], [0, intensityF(ii)], 'color', 'g');

  % shepard F#
  fFs(ii) = note2freq('F#', ii-1, 'equal');
  ratio = C4/fFs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityFs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fFs(ii)) log10(fFs(ii))], [0, intensityFs(ii)], 'color', 'r');

  % shepard G
  fG(ii) = note2freq('G', ii-2, 'equal');
  ratio = C4/fG(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityG(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fG(ii)) log10(fG(ii))], [0, intensityG(ii)], 'color', 'g');

  % shepard G#
  fGs(ii) = note2freq('G#', ii-2, 'equal');
  ratio = C4/fGs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityGs(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fGs(ii)) log10(fGs(ii))], [0, intensityGs(ii)], 'color', 'r');

  % shepard A
  fA(ii) = note2freq('A', ii-2, 'equal');
  ratio = C4/fA(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityA(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fA(ii)) log10(fA(ii))], [0, intensityA(ii)], 'color', 'g');

  % shepard A#
  fAs(ii) = note2freq('A#', ii-2, 'equal');
  ratio = C4/fAs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityAs(ii) = exp(-0.5*(diff/sigma)^2);
  %line([log10(fAs(ii)) log10(fAs(ii))], [0, intensityAs(ii)], 'color', 'r');

  % shepard B
  fB(ii) = note2freq('B', ii-2, 'equal');
  ratio = C4/fB(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityB(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fB(ii)) log10(fB(ii))], [0, intensityB(ii)], 'color', 'g');

end

% make plot
hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fC), intensityC, 'k*')
hold on
plot(log10(fD), intensityD, 'b*')
hold on
plot(log10(fFs), intensityFs, 'r*')
hold on
plot(log10(fB), intensityB, 'g*')

xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0','C1','C2','C3','C4','C5','C6','C7','C8'});

print('-depsc2', 'shepardintensity.eps')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find max intensities values and frequencies
[maxC,ndx] = max(intensityC);
maxfC = fC(ndx)

[maxCs,ndx] = max(intensityCs);
maxfCs = fCs(ndx)

[maxD,ndx] = max(intensityD);
maxfD = fD(ndx)

[maxDs,ndx] = max(intensityDs);
maxfDs = fDs(ndx)

[maxE,ndx] = max(intensityE);
maxfE = fE(ndx)

[maxF,ndx] = max(intensityF);
maxfF = fF(ndx)

[maxFs,ndx] = max(intensityFs);
maxfFs = fFs(ndx)

[maxG,ndx] = max(intensityG);
maxfG = fG(ndx)

[maxGs,ndx] = max(intensityGs);
maxfGs = fGs(ndx)

[maxA,ndx] = max(intensityA);
maxfA = fA(ndx)

[maxAs,ndx] = max(intensityAs);
maxfAs = fAs(ndx)

[maxB,ndx] = max(intensityB);
maxfB = fB(ndx)

