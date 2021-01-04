% calculate frequency ratios in just temperament 
% and compare to equal temperament

fifth = 3/2;
third = 5/4;
octave = 2;
semitone = 2^(1/12);

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

% just freq ratios
C1.just = 1;
E.just  = third;
Ab.just = third^(-1)*octave;

G.just  = fifth;
B.just  = fifth*third;
Eb.just = fifth*third^(-1);

D.just  = fifth^2*octave^(-1);
Fs.just = fifth^2*third*octave^(-1);
Bb.just = fifth^2*third^(-1);

F.just  = fifth^(-1)*octave;
A.just  = fifth^(-1)*third*octave;
Cs.just = fifth^(-1)*third^2;

C2.just = 2;

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
C1.diff = ratio2cents(C1.just/C1.equal);
Cs.diff = ratio2cents(Cs.just/Cs.equal);
D.diff  = ratio2cents(D.just/D.equal);
Eb.diff = ratio2cents(Eb.just/Eb.equal);
E.diff  = ratio2cents(E.just/E.equal);
F.diff  = ratio2cents(F.just/F.equal);
Fs.diff = ratio2cents(Fs.just/Fs.equal);
G.diff  = ratio2cents(G.just/G.equal);
Ab.diff = ratio2cents(Ab.just/Ab.equal);
A.diff  = ratio2cents(A.just/A.equal);
Bb.diff = ratio2cents(Bb.just/Bb.equal);
B.diff  = ratio2cents(B.just/B.equal);
C2.diff = ratio2cents(C2.just/C2.equal);

fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        C1.name, C1.just, C1.equal, round(C1.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        Cs.name, Cs.just, Cs.equal, round(Cs.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        D.name, D.just, D.equal, round(D.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        Eb.name, Eb.just, Eb.equal, round(Eb.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        E.name, E.just, E.equal, round(E.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        F.name, F.just, F.equal, round(F.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        Fs.name, Fs.just, Fs.equal, round(Fs.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        G.name, G.just, G.equal, round(G.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        Ab.name, Ab.just, Ab.equal, round(Ab.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        A.name, A.just, A.equal, round(A.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        Bb.name, Bb.just, Bb.equal, round(Bb.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        B.name, B.just, B.equal, round(B.diff));
fprintf('%3s: just = %.3f, equal = %.3f, diff = %d cents\n', ...
        C2.name, C2.just, C2.equal, round(C2.diff));

return

