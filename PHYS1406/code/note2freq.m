function f = note2freq(note, level, temp)
%
% convert note to frequency 
%
% Input:
%    note  - e.g., 'C', 'C#', 'Eb', ...
%    level - e.g., 4 for C4
%    temp  - temperament e.g., equal, pyth, just, mean
%
% Output:
%    f    - frequency (Hz)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% just intervals
octave = 2/1;
fifth  = 3/2;
fourth = 4/3;
third = 5/4;

% for mean-tone temperament
delta = 81/80; % syntonic comma
dfifth = fifth*delta^(-1/4); % diminished fifth
    
% for equal temperament
semitone = 2^(1/12);

% reference frequency (A4 = 440 Hz for all tuning systems)
f_A4 = 440; % Hz

switch temp

  case 'equal'

    % calculate f_C4 given A4
    f_C4 = f_A4/(semitone^9);

    % calculate frequency ratio 
    switch note
      case 'C'
         ratio = semitone^0;
      case {'C#','Db'} 
         ratio = semitone^1;
      case 'D'
         ratio = semitone^2;
      case {'Eb','D#'} 
         ratio = semitone^3;
      case 'E'
         ratio = semitone^4;
      case {'F', 'E#'}
         ratio = semitone^5;
      case {'F#','Gb'}
         ratio = semitone^6;
      case 'G'
         ratio = semitone^7;
      case {'Ab','G#'} 
         ratio = semitone^8;
      case 'A'
         ratio = semitone^9;
      case {'Bb','A#'}
         ratio = semitone^10;
      case 'B'
         ratio = semitone^11;
      case 'B#'
         ratio = semitone^12;
      otherwise
        error('unknown note\n')
    end

  case 'pyth'

    % calculate f_C4 given A4
    f_C4 = f_A4/(fifth^3*octave^(-1));

    % calculate frequency ratio 
    switch note
      case 'C'
         ratio = fifth^0*octave^(0);
      case 'G'
         ratio = fifth^1*octave^(0);
      case 'D'
         ratio = fifth^2*octave^(-1);
      case 'A'
         ratio = fifth^3*octave^(-1);
      case 'E'
         ratio = fifth^4*octave^(-2);
      case 'B'
         ratio = fifth^5*octave^(-2);
      case 'F#'
         ratio = fifth^6*octave^(-3);
      case 'C#' 
         ratio = fifth^7*octave^(-4);
      case 'G#' 
         ratio = fifth^8*octave^(-4);
      case 'D#' 
         ratio = fifth^9*octave^(-5);
      case 'A#' 
         ratio = fifth^10*octave^(-5);
      case 'E#' 
         ratio = fifth^11*octave^(-6);
      case 'B#' 
         ratio = fifth^12*octave^(-6);
      case 'F'
         ratio = fifth^(-1)*octave^1;
      case 'Bb'
         ratio = fifth^(-2)*octave^2;
      case 'Eb' 
         ratio = fifth^(-3)*octave^2;
      case 'Ab' 
         ratio = fifth^(-4)*octave^3;
      case 'Db' 
         ratio = fifth^(-5)*octave^3;
      case 'Gb' 
         ratio = fifth^(-6)*octave^4;
      case 'Cb' 
         ratio = fifth^(-7)*octave^5;
      case 'Fb' 
         ratio = fifth^(-8)*octave^5;
      case 'Bbb' 
         ratio = fifth^(-9)*octave^6;
      case 'Ebb' 
         ratio = fifth^(-10)*octave^6;
      case 'Abb' 
         ratio = fifth^(-11)*octave^7;
      case 'Dbb' 
         ratio = fifth^(-12)*octave^7;
      otherwise
        error('unknown note\n')
    end

  case 'mean'

    % calculate f_C4 given A4
    f_C4 = f_A4/(dfifth^3*octave^(-1));

    % calculate frequency ratio 
    switch note
      case 'C'
         ratio = dfifth^0*octave^(0);
      case 'G'
         ratio = dfifth^1*octave^(0);
      case 'D'
         ratio = dfifth^2*octave^(-1);
      case 'A'
         ratio = dfifth^3*octave^(-1);
      case 'E'
         ratio = dfifth^4*octave^(-2);
      case 'B'
         ratio = dfifth^5*octave^(-2);
      case {'F#','Gb'}
         ratio = dfifth^6*octave^(-3);
      case {'C#','Db'} 
         ratio = dfifth^7*octave^(-4);
      case 'F'
         ratio = dfifth^(-1)*octave^1;
      case {'Bb','A#'}
         ratio = dfifth^(-2)*octave^2;
      case {'Eb','D#'} 
         ratio = dfifth^(-3)*octave^2;
      case {'Ab','G#'} 
         ratio = dfifth^(-4)*octave^3;
      otherwise
        error('unknown note\n')
    end

  case 'just'

    % calculate f_C4 given A4
    f_C4 = f_A4/(fifth^(-1)*third*octave);

    % calculate frequency ratio 
    switch note
      case 'C'
         ratio = 1;
      case 'E'
         ratio = third;
      case {'Ab','G#'} 
         ratio = third^(-1)*octave;
      case 'G'
         ratio = fifth;
      case 'B'
         ratio = fifth*third;
      case {'Eb','D#'} 
         ratio = fifth*third^(-1);
      case 'D'
         ratio = fifth^2*octave^(-1);
      case {'F#','Gb'}
         ratio = fifth^2*third*octave^(-1);
      case {'Bb','A#'}
         ratio = fifth^2*third^(-1);
      case 'F'
         ratio = fifth^(-1)*octave;
      case 'A'
         ratio = fifth^(-1)*third*octave;
      case {'C#','Db'} 
         ratio = fifth^(-1)*third^2;
      otherwise
        error('unknown note\n')
    end

  otherwise
    error('unknown temperament')

end

% calculate absolute frequency
f = f_C4*ratio*2^(level-4);

return

