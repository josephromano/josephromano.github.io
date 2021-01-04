function c = ratio2cents(r)
%
% convert ratio of frequencies to an interval in cents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 100 cents = 1 semitone
% R = 2^(n/1200)

c = (1200/log10(2))*log10(r);

return

