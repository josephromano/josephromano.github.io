% script for fouriersynthesizesound

type = input('enter type of wave:\n  sine, triangle, square, sawtooth, sawtoothphase, missing\n');

switch type
  case 'sine'
    % pure sine
    amplitude = [1];
    phase = [0];

  case 'triangle'
    % triangle wave
    amplitude = [1 0 1/9 0 1/25 0 1/49 0 1/81 0];
    phase = [0 0 180 0 0 0 180 0 0 0];
     
  case 'square' 
    % square wave
    amplitude = [1 0 1/3 0 1/5 0 1/7 0 1/9 0];
    phase = [0 0 0 0 0 0 0 0 0 0];

  case 'sawtooth' 
    % sawtooth wave
    amplitude = [1 1/2 1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10];
    phase = [0 0 0 0 0 0 0 0 0 0];

  case 'sawtoothphase' 
    % sawtooth wave with alternating phase
    amplitude = [1 1/2 1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10];
    phase = [0 90 0 90 0 90 0 90 0 90];

  case 'missing'
    % missing fundamental
    amplitude = [0 1/2 1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10];
    phase = [0 0 0 0 0 0 0 0 0 0];

  otherwise
    error('unrecognized type');

end

fouriersynthesizesound(amplitude, phase);

