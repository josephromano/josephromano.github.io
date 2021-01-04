% calculate frequency ratios for musical intervals 
% and compare to equal temperament

ET_semitone = 2^(1/12);

octave.name = 'octave';
fifth.name = 'fifth';
fourth.name = 'fourth';
Mthird.name = 'major third';
mthird.name = 'minor third';
msixth.name = 'minor sixth';
Msixth.name = 'major sixth';

octave.ratio = 2;
fifth.ratio = 3/2;
fourth.ratio = 4/3;
Mthird.ratio = 5/4;
mthird.ratio = 6/5;
msixth.ratio = 8/5;
Msixth.ratio = 5/3;

octave.ET = ET_semitone^12;
fifth.ET = ET_semitone^7;
fourth.ET = ET_semitone^5;
Mthird.ET = ET_semitone^4;
mthird.ET = ET_semitone^3;
msixth.ET = ET_semitone^8;
Msixth.ET = ET_semitone^9;

octave.diff = ratio2cents(octave.ratio/octave.ET);
fifth.diff = ratio2cents(fifth.ratio/fifth.ET);
fourth.diff = ratio2cents(fourth.ratio/fourth.ET);
Mthird.diff = ratio2cents(Mthird.ratio/Mthird.ET);
mthird.diff = ratio2cents(mthird.ratio/mthird.ET);
msixth.diff = ratio2cents(msixth.ratio/msixth.ET);
Msixth.diff = ratio2cents(Msixth.ratio/Msixth.ET);

fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        octave.name, octave.ratio, octave.ET, round(octave.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        fifth.name, fifth.ratio, fifth.ET, round(fifth.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        fourth.name, fourth.ratio, fourth.ET, round(fourth.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        Mthird.name, Mthird.ratio, Mthird.ET, round(Mthird.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        mthird.name, mthird.ratio, mthird.ET, round(mthird.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        msixth.name, msixth.ratio, msixth.ET, round(msixth.diff));
fprintf('%11s: ratio = %.3f, ET ratio = %.3f, diff = %d cents\n', ...
        Msixth.name, Msixth.ratio, Msixth.ET, round(Msixth.diff));

return

