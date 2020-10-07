function N = calNl(lmax)
%
% calculate Nl values packed an lm array with m running from 
% -lmax to +lmax.
%
% output structure has fields
%
%   N.data = Nl;
%   N.lvec = lvec;
%   N.mvec = mvec;
%   N.lmax = lmax;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize variables for m>=0, l>=2
numRows = (lmax+1)*(lmax+2)/2 - 3;
lvec = zeros(numRows,1);
mvec = zeros(numRows,1);
Nl = zeros(numRows,1);

% loop over m>=0 values
ii = 1;
for m=0:lmax

  for l = max(2,m):lmax
    mvec(ii) = m;
    lvec(ii) = l;

    % calculate Nl
    Nl(ii) = sqrt(2*factorial(l-2)/factorial(l+2));

    % increment counter
    ii = ii+1;
  end

end

% extend arrays to have values for negative m
Nl = [ flipud(Nl(lmax:end)); Nl];

% extend vectors of l and m values to negative values of m
lvec = [ flipud(lvec(lmax:end)); lvec];
mvec = [-flipud(mvec(lmax:end)); mvec];

% fill structure
N.data = Nl;
N.lvec = lvec;
N.mvec = mvec;
N.lmax = lmax;

return

