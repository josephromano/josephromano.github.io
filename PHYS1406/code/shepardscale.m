function shepardscale(type, dur)
%
% play shepard tone scale
%
% type = 'ascending' or 'descending'
% dur  = duration of notes in seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch type

  case 'ascending'
    for ii=1:4
      playshepardtone('C',dur);
      playshepardtone('C#',dur);
      playshepardtone('D',dur);
      playshepardtone('Eb',dur);
      playshepardtone('E',dur);
      playshepardtone('F',dur);
      playshepardtone('F#',dur);
      playshepardtone('G',dur);
      playshepardtone('Ab',dur);
      playshepardtone('A',dur);
      playshepardtone('Bb',dur);
      playshepardtone('B',dur);
    end

  case 'descending'
    for ii=1:4
      playshepardtone('C',dur);
      playshepardtone('B',dur);
      playshepardtone('Bb',dur);
      playshepardtone('A',dur);
      playshepardtone('Ab',dur);
      playshepardtone('G',dur);
      playshepardtone('F#',dur);
      playshepardtone('F',dur);
      playshepardtone('E',dur);
      playshepardtone('Eb',dur);
      playshepardtone('D',dur);
      playshepardtone('C#',dur);
    end

  otherwise
    error('unknown type')

end

return
