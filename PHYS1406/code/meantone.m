% calculate frequency ratios in mean-tone temperament 
% and compare to equal temperament

delta = 81/80; % syntonic comma
fifth = 3/2; % fifth
dfifth = fifth*delta^(-1/4); % diminished fifth
semitone = 2^(1/12);
octave = 2;

% note names
C1.name = 'C1';
Cs.name = 'C#';
D.name  = 'D';
Eb.name = 'Eb';
E.name  = 'E';
F.name  = 'F';
Fs.name = 'F#';
G.name  = 'G';
Ab.name = 'Ab';
A.name  = 'A';
Bb.name = 'Bb';
B.name  = 'B';
C2.name = 'C2';

% mean freq ratios
C1.mean = dfifth^0*octave^(0);
G.mean  = dfifth^1*octave^(0);
D.mean  = dfifth^2*octave^(-1);
A.mean  = dfifth^3*octave^(-1);
E.mean  = dfifth^4*octave^(-2);
B.mean  = dfifth^5*octave^(-2);
Fs.mean = dfifth^6*octave^(-3);
Cs.mean = dfifth^7*octave^(-4);
F.mean  = dfifth^(-1)*octave^1;
Bb.mean = dfifth^(-2)*octave^2;
Eb.mean = dfifth^(-3)*octave^2;
Ab.mean = dfifth^(-4)*octave^3;
C2.mean = 2;

% equal-temperament freq ratios
C1.equal = 1;
Cs.equal = semitone^1;
D.equal  = semitone^2;
Eb.equal = semitone^3;
E.equal  = semitone^4;
F.equal  = semitone^5;
Fs.equal = semitone^6;
G.equal  = semitone^7;
Ab.equal = semitone^8;
A.equal  = semitone^9;
Bb.equal = semitone^10;
B.equal  = semitone^11;
C2.equal = semitone^12;

% difference
C1.diff = ratio2cents(C1.mean/C1.equal);
Cs.diff = ratio2cents(Cs.mean/Cs.equal);
D.diff  = ratio2cents(D.mean/D.equal);
Eb.diff = ratio2cents(Eb.mean/Eb.equal);
E.diff  = ratio2cents(E.mean/E.equal);
F.diff  = ratio2cents(F.mean/F.equal);
Fs.diff = ratio2cents(Fs.mean/Fs.equal);
G.diff  = ratio2cents(G.mean/G.equal);
Ab.diff = ratio2cents(Ab.mean/Ab.equal);
A.diff  = ratio2cents(A.mean/A.equal);
Bb.diff = ratio2cents(Bb.mean/Bb.equal);
B.diff  = ratio2cents(B.mean/B.equal);
C2.diff = ratio2cents(C2.mean/C2.equal);

fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        C1.name, C1.mean, C1.equal, round(C1.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        Cs.name, Cs.mean, Cs.equal, round(Cs.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        D.name, D.mean, D.equal, round(D.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        Eb.name, Eb.mean, Eb.equal, round(Eb.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        E.name, E.mean, E.equal, round(E.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        F.name, F.mean, F.equal, round(F.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        Fs.name, Fs.mean, Fs.equal, round(Fs.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        G.name, G.mean, G.equal, round(G.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        Ab.name, Ab.mean, Ab.equal, round(Ab.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        A.name, A.mean, A.equal, round(A.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        Bb.name, Bb.mean, Bb.equal, round(Bb.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        B.name, B.mean, B.equal, round(B.diff));
fprintf('%3s: mean = %.3f, equal = %.3f, diff = %d cents\n', ...
        C2.name, C2.mean, C2.equal, round(C2.diff));

return

