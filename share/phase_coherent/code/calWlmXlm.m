function [W, X] = calWlmXlm(theta, phi, lmax)
%
% calculate Wlm, Xlm values for single sky position packed as an
% lm array with m running from -lmax to +lmax
%
% output structures have fields
%   W.data = Wlm;
%   W.lvec = lvec;
%   W.mvec = mvec;
%   W.lmax = lmax;
%
%   X.data = Xlm;
%   X.lvec = lvec;
%   X.mvec = mvec;
%   X.lmax = lmax;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = cos(theta);

% initialize variables for m>=0, l>=2
numRows = (lmax+1)*(lmax+2)/2 - 3;
lvec = zeros(numRows,1);
mvec = zeros(numRows,1);
Wlm = zeros(numRows,1);
Xlm = zeros(numRows,1);

% loop over m>=0 values
ii = 1;
for m=0:lmax

  for l = max(2,m):lmax
    mvec(ii) = m;
    lvec(ii) = l;

    % calculate Wlm, Xlm
    Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );      
    [Gp, Gm] = Glm(l, m, x);
    Wlm(ii) = 2*Nlm*Gp*exp(i*m*phi);
    Xlm(ii) = 2*sqrt(-1)*Nlm*Gm*exp(i*m*phi);

    % increment counter
    ii = ii+1;
  end

end

% extend arrays to have values for negative m
Wlmm = conj(Wlm).*(-1).^mvec;
Wlm = [flipud(Wlmm(lmax:end,:)); Wlm]; 
Xlmm = conj(Xlm).*(-1).^mvec;
Xlm = [flipud(Xlmm(lmax:end,:)); Xlm]; 

% extend vectors of l and m values to negative values of m
lvec = [ flipud(lvec(lmax:end)); lvec];
mvec = [-flipud(mvec(lmax:end)); mvec];

% fill structures
W.data = Wlm;
W.lvec = lvec;
W.mvec = mvec;
W.lmax = lmax;

X.data = Xlm;
X.lvec = lvec;
X.mvec = mvec;
X.lmax = lmax;

return
