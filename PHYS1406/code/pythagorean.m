% calculate frequency ratios in pythagorean temperament
% and compare to equal temperament

fifth = 3/2;
octave = 2;
ET_fifth = 2^(7/12);

up_note{1}='C1';
up_note{2}='G';
up_note{3}='D';
up_note{4}='A';
up_note{5}='E';
up_note{6}='B';
up_note{7}='F#';
up_note{8}='C#';
up_note{9}='G#';
up_note{10}='D#';
up_note{11}='A#';
up_note{12}='E#';
up_note{13}='B#';

down_note{1}='C2';
down_note{2}='F';
down_note{3}='Bb';
down_note{4}='Eb';
down_note{5}='Ab';
down_note{6}='Db';
down_note{7}='Gb';
down_note{8}='Cb';
down_note{9}='Fb';
down_note{10}='Bbb';
down_note{11}='Ebb';
down_note{12}='Abb';
down_note{13}='Dbb';

% calculate frequency ratios
for ii=1:13

  % calculate exponent appropriately
  if ii<=7
    p = floor((ii+1)/2);
  else
    p = 1+floor(ii/2);
  end

  pyth_up(ii) = fifth^(ii-1)*octave^(1-p);
  pyth_down(ii) = fifth^(1-ii)*octave^p;

  ET_up(ii) = ET_fifth^(ii-1)*octave^(1-p);
  ET_down(ii) = ET_fifth^(1-ii)*octave^p;

  diff_up(ii) = ratio2cents(pyth_up(ii)/ET_up(ii));
  diff_down(ii) = ratio2cents(pyth_down(ii)/ET_down(ii));

end

% display upward notes, freqs
for ii=1:13
  fprintf('%3s: pyth = %.3f, ET = %.3f, diff = %d cents\n', ...
          up_note{ii}, pyth_up(ii), ET_up(ii), round(diff_up(ii)))
end
fprintf('\n');

% display downward notes, freqs
for ii=1:13
  fprintf('%3s: pyth = %.3f, ET = %.3f, diff = %d cents\n', ...
          down_note{ii}, pyth_down(ii), ET_down(ii), round(diff_down(ii)))
end
fprintf('\n');

% sort frequencies and notes
pyth_freqs = [pyth_up pyth_down];
[pyth_freqs_sorted, ndx] = sort(pyth_freqs);

ET_freqs = [ET_up ET_down];
ET_freqs_sorted = ET_freqs(ndx);

diff = [diff_up diff_down];
diff_sorted = diff(ndx);

notes = [up_note down_note];
notes_sorted = notes(ndx);

% display sorted notes, freqs
for ii=1:26
  fprintf('%3s: pyth = %.3f, ET = %.3f, diff = %d cents\n', ...
          notes_sorted{ii}, pyth_freqs_sorted(ii), ET_freqs_sorted(ii), round(diff_sorted(ii)))
end

return

