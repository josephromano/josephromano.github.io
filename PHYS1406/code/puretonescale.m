function puretonescale(type, dur)
%
% play pure tone scale
%
% type = 'ascending' or 'descending'
% dur  = duration of notes in seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch type

  case 'ascending'
    for ii=1:4
      playnote('C', 4, 'equal', dur);
      playnote('C#', 4, 'equal', dur);
      playnote('D', 4, 'equal', dur);
      playnote('Eb', 4, 'equal', dur);
      playnote('E', 4, 'equal', dur);
      playnote('F', 4, 'equal', dur);
      playnote('F#', 4, 'equal', dur);
      playnote('G', 4, 'equal', dur);
      playnote('Ab', 4, 'equal', dur);
      playnote('A', 4, 'equal', dur);
      playnote('Bb', 4, 'equal', dur);
      playnote('B', 4, 'equal', dur);
    end

  case 'descending'
    for ii=1:4
      playnote('C', 4, 'equal', dur);
      playnote('B', 4, 'equal', dur);
      playnote('Bb', 4, 'equal', dur);
      playnote('A', 4, 'equal', dur);
      playnote('Ab', 4, 'equal', dur);
      playnote('G', 4, 'equal', dur);
      playnote('F#', 4, 'equal', dur);
      playnote('F', 4, 'equal', dur);
      playnote('E', 4, 'equal', dur);
      playnote('Eb', 4, 'equal', dur);
      playnote('D', 4, 'equal', dur);
      playnote('C#', 4, 'equal', dur);
    end

  otherwise
    error('unknown type')

end

return
