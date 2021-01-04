function fifths(level)
% compare two different fifths (C-G) and (C#-Ab) in
% different temperaments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate freqs for the different notes
C_equal = note2freq('C',level,'equal');
C_pyth  = note2freq('C',level,'pyth');
C_just  = note2freq('C',level,'just');

G_equal = note2freq('G',level,'equal');
G_pyth  = note2freq('G',level,'pyth');
G_just  = note2freq('G',level,'just');

Cs_equal = note2freq('C#',level,'equal');
Cs_pyth  = note2freq('C#',level,'pyth');
Cs_just  = note2freq('C#',level,'just');

Ab_equal = note2freq('Ab',level,'equal');
Ab_pyth  = note2freq('Ab',level,'pyth');
Ab_just  = note2freq('Ab',level,'just');

% calculate difference from a perfect fifth of 3/2
CG_equal_diff = ratio2cents(G_equal/C_equal*2/3);
CG_pyth_diff = ratio2cents(G_pyth/C_pyth*2/3);
CG_just_diff = ratio2cents(G_just/C_just*2/3);

CsAb_equal_diff = ratio2cents(Ab_equal/Cs_equal*2/3);
CsAb_pyth_diff = ratio2cents(Ab_pyth/Cs_pyth*2/3);
CsAb_just_diff = ratio2cents(Ab_just/Cs_just*2/3);

% display results
fprintf('C-G   equal = %.3f, diff = %d cents\n', G_equal/C_equal, round(CG_equal_diff));
fprintf('C-G   pyth  = %.3f, diff = %d cents\n', G_pyth/C_pyth, round(CG_pyth_diff));
fprintf('C-G   just  = %.3f, diff = %d cents\n', G_just/C_just, round(CG_just_diff));

fprintf('C#-Ab equal = %.3f, diff = %d cents\n', Ab_equal/Cs_equal, round(CsAb_equal_diff));
fprintf('C#-Ab pyth  = %.3f, diff = %d cents\n', Ab_pyth/Cs_pyth, round(CsAb_pyth_diff));
fprintf('C#-Ab just  = %.3f, diff = %d cents\n', Ab_just/Cs_just, round(CsAb_just_diff));

% play intervals
playinterval('C', level, 'equal', 'G', level, 'equal');
playinterval('C', level, 'pyth', 'G', level, 'pyth');
playinterval('C', level, 'just', 'G', level, 'just');

playinterval('C#', level, 'equal', 'Ab', level, 'equal');
playinterval('C#', level, 'pyth', 'Ab', level, 'pyth');
playinterval('C#', level, 'just', 'Ab', level, 'just');

return

