function little_lamb(key, level, temp)
%
% play "mary had a little lamb" in different major keys
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

    playnote('E', level, temp, q);
    playnote('D', level, temp, q);
    playnote('C', level, temp, q);
    playnote('D', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, h);
    playnote('D', level, temp, q);
    playnote('D', level, temp, q);
    playnote('D', level, temp, h);
    playnote('E', level, temp, q);
    playnote('G', level, temp, q);
    playnote('G', level, temp, h);

    playnote('E', level, temp, q);
    playnote('D', level, temp, q);
    playnote('C', level, temp, q);
    playnote('D', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, q);
    playnote('E', level, temp, h);
    playnote('D', level, temp, q);
    playnote('D', level, temp, q);
    playnote('E', level, temp, q);
    playnote('D', level, temp, q);
    playnote('C', level, temp, h);

  case 'C#'

    playnote('F', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('C#', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, h);
    playnote('Eb', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('Eb', level, temp, h);
    playnote('F', level, temp, q);
    playnote('Ab', level, temp, q);
    playnote('Ab', level, temp, h);

    playnote('F', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('C#', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, q);
    playnote('F', level, temp, h);
    playnote('Eb', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('F', level, temp, q);
    playnote('Eb', level, temp, q);
    playnote('C#', level, temp, h);

end

return

