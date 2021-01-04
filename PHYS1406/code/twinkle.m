function little_lamb(key, level, temp)
%
% play "twinkle, twinkle"
% and temperaments
%
% key = 'C' or 'C#'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% duration of notes (in sec)
q = 0.35;
h = 0.70;

switch key

  case 'C'

    playnote('C', level, temp, q);
    playnote('C', level, temp, q);
    playnote('G', level, temp, q);
    playnote('G', level, temp, q);
    playnote('A', level, temp, q);
    playnote('A', level, temp, q);
    playnote('G', level, temp, h);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('D', level, temp, q);
    playnote('D', level, temp, q);
    playnote('C', level, temp, h);

    playnote('G', level, temp, q);
    playnote('G', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('D', level, temp, h);
    playnote('G', level, temp, q);
    playnote('G', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('D', level, temp, h);

    playnote('C', level, temp, q);
    playnote('C', level, temp, q);
    playnote('G', level, temp, q);
    playnote('G', level, temp, q);
    playnote('A', level, temp, q);
    playnote('A', level, temp, q);
    playnote('G', level, temp, h);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('D', level, temp, q);
    playnote('D', level, temp, q);
    playnote('C', level, temp, h);

  case 'C#'

    playnote('C#', level, temp, q);
    playnote('C#', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('Bb', level, temp, q);
    playnote('Bb', level, temp, q);
    playnote('Ab', level, temp, h);
    playnote('F#', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('C#', level, temp, h);

    playnote('Ab', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('Eb', level, temp, h);
    playnote('Ab', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('Eb', level, temp, h);

    playnote('C#', level, temp, q);
    playnote('C#', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('Bb', level, temp, q);
    playnote('Bb', level, temp, q);
    playnote('Ab', level, temp, h);
    playnote('F#', level, temp, q);
    playnote('F#', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('C#', level, temp, h);

end

return

