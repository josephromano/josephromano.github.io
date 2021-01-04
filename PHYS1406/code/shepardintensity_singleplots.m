% plot gaussian intensity distribution for shepard tones
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% gaussian intensity distribution for partials
% (peak at C4, sigma in octaves)
C4 = note2freq('C', 4, 'equal');
sigma = 1.5; % in octaves
maxlevel = 9;

% smooth gaussian
Nf = 1e5;
freqs = linspace(10,1e4,Nf);
p = exp(-0.5*((log10(freqs)-log10(C4))/(sigma*log10(2))).^2);

% calculate frequencies and intensities for notes an octave apart
figure
for ii=1:maxlevel

  % shepard C
  fC(ii) = note2freq('C', ii-1, 'equal');
  ratio = C4/fC(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityC(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fC(ii)) log10(fC(ii))], [0, intensityC(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fC), intensityC, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard C Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardC.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard C#
  fCs(ii) = note2freq('C#', ii-1, 'equal');
  ratio = C4/fCs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityCs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fCs(ii)) log10(fCs(ii))], [0, intensityCs(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fCs), intensityCs, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard C# Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardCs.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard D
  fD(ii) = note2freq('D', ii-1, 'equal');
  ratio = C4/fD(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityD(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fD(ii)) log10(fD(ii))], [0, intensityD(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fD), intensityD, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard D Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardD.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard D#
  fDs(ii) = note2freq('D#', ii-1, 'equal');
  ratio = C4/fDs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityDs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fDs(ii)) log10(fDs(ii))], [0, intensityDs(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fDs), intensityDs, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard D# Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardDs.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard E
  fE(ii) = note2freq('E', ii-1, 'equal');
  ratio = C4/fE(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityE(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fE(ii)) log10(fE(ii))], [0, intensityE(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fE), intensityE, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard E Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardE.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard F
  fF(ii) = note2freq('F', ii-1, 'equal');
  ratio = C4/fF(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityF(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fF(ii)) log10(fF(ii))], [0, intensityF(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fF), intensityF, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard F Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardF.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard F#
  fFs(ii) = note2freq('F#', ii-1, 'equal');
  ratio = C4/fFs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityFs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fFs(ii)) log10(fFs(ii))], [0, intensityFs(ii)], 'color', 'r');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fFs), intensityFs, 'r*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard F# Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardFs.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard G
  fG(ii) = note2freq('G', ii-2, 'equal');
  ratio = C4/fG(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityG(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fG(ii)) log10(fG(ii))], [0, intensityG(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fG), intensityG, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard G Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardG.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard G#
  fGs(ii) = note2freq('G#', ii-2, 'equal');
  ratio = C4/fGs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityGs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fGs(ii)) log10(fGs(ii))], [0, intensityGs(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fGs), intensityGs, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard G# Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardGs.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard A
  fA(ii) = note2freq('A', ii-2, 'equal');
  ratio = C4/fA(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityA(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fA(ii)) log10(fA(ii))], [0, intensityA(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fA), intensityA, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard A Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardA.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard A#
  fAs(ii) = note2freq('A#', ii-2, 'equal');
  ratio = C4/fAs(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityAs(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fAs(ii)) log10(fAs(ii))], [0, intensityAs(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fAs), intensityAs, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard A# Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardAs.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii=1:maxlevel

  % shepard B
  fB(ii) = note2freq('B', ii-2, 'equal');
  ratio = C4/fB(ii);
  diff = ratio2cents(ratio)/1200; % difference in octaves
  intensityB(ii) = exp(-0.5*(diff/sigma)^2);
  line([log10(fB(ii)) log10(fB(ii))], [0, intensityB(ii)], 'color', 'b');

end

hold on
plot(log10(freqs), p, 'k');
hold on
plot(log10(fB), intensityB, 'b*')
xlabel('log10(frequency)')
ylabel('intensity (arbitrary units)')
title('Shepard B Tone')

% label octaves
for ii=1:maxlevel
  C = note2freq('C', ii-1, 'equal');
  octave(ii) = log10(C);
end
set(gca,'XTick',octave);
set(gca,'XTickLabel',{'C0', 'C1','C2','C3','C4','C5','C6','C7', 'C8'});
print -depsc2 'shepardB.eps'

return
