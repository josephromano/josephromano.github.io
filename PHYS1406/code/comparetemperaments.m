% compare equal, pythagorean, and just temperaments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% note names
note{1} = 'C';
note{2} = 'C#';
note{3}  = 'D';
note{4} = 'Eb';
note{5}  = 'E';
note{6}  = 'F';
note{7} = 'F#';
note{8}  = 'G';
note{9} = 'Ab';
note{10}  = 'A';
note{11} = 'Bb';
note{12}  = 'B';
note{13} = 'C';

%f=figure('Position',[300 400 1100 200])
figure
xlim([0 1200]);
ylim([0 1.25]);
xlabel('freq ratio (cents)');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 4;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[300 400 xSize*50 ySize*50])

% calculate C4
C4_equal = note2freq(note{1}, 4, 'equal');
C4_pyth = note2freq(note{1}, 4, 'pyth');
C4_just = note2freq(note{1}, 4, 'just');

% first note (all temperaments agree here)
equal(1) = round(ratio2cents(1));
pyth(1) =  round(ratio2cents(1));
just(1) =  round(ratio2cents(1));
line([equal(1) equal(1)], [0.00 1.00], 'color', 'k', 'linestyle', ':');
line([pyth(1) pyth(1)],   [0.50 1.00], 'color', 'b', 'linestyle', '-');
line([just(1) just(1)],   [0.00 0.50], 'color', 'r', 'linestyle', '-');
text(equal(1)-10, 1.15, note{1}, 'fontsize', 18)

% calculate freq and cent value of notes
for ii=2:12
  ratio = note2freq(note{ii}, 4, 'equal')/C4_equal;
  equal(ii) = round(ratio2cents(ratio));
  line([equal(ii) equal(ii)], [0.00 1.00], 'color', 'k', 'linestyle', ':');

  ratio = note2freq(note{ii}, 4, 'pyth')/C4_pyth;
  pyth(ii) = round(ratio2cents(ratio));
  line([pyth(ii) pyth(ii)], [0.50 1.00], 'color', 'b', 'linestyle', '-');

  ratio = note2freq(note{ii}, 4, 'just')/C4_just;
  just(ii) = round(ratio2cents(ratio));
  line([just(ii) just(ii)], [0.00 0.50], 'color', 'r', 'linestyle', '-');

  text(equal(ii)-10, 1.15, note{ii}, 'fontsize', 18)
end

% final note one octave higher (all temperaments agree here)
equal(13) = round(ratio2cents(2));
pyth(13)  = round(ratio2cents(2));
just(13)  = round(ratio2cents(2));
line([equal(13) equal(13)], [0.00 1.00], 'color', 'k', 'linestyle', ':');
line([pyth(13) pyth(13)],   [0.50 1.00], 'color', 'b', 'linestyle', '-');
line([just(13) just(13)],   [0.00 0.50], 'color', 'r', 'linestyle', '-');
text(equal(13)-10, 1.15, note{13}, 'fontsize', 18)

print('-depsc2', 'comparetemperaments.eps')

return
