function R_Plm = calRPlm_ifo(f, u, v, x0, lmax)
%
% calculates R^P_lm values for a single frequency and a single ifo 
% packed as an lm array with grad and curl blocks (Nt x Nm)
%
% Inputs:
%   f     - array of frequencies (Hz)
%   u, v  - unit vectors along the arms of the interferometer
%           (these can be arrays u(1,:), u(2,:), u(3,:), etc.
%           corresponding to different discrete times t)
%   x0    - position vector of detector
%   lmax  - max value of L
%
% Output:
%
%   R_Plm is a structure containing the following fields:
%
%   R_Plm.data = [R_Glm, R_Clm];
%   R_Plm.lvec = lvec 
%   R_Plm.mvec = mvec    
%   R_Plm.lmax = lmax
%
%   An lm array is indexed with m running from -lmax to +lmax
%   example: lmax = 5
%
%   m=-5, l=5
%   m=-4, l=5
%   m=-4, l=4
%   m=-3, l=5
%   m=-3, l=4
%   m=-3, l=3
%   ...
%   m=-1, l=2
%   m=0,  l=2
%   m=0,  l=3
%   m=0,  l=4
%   m=0,  l=5
%   m=1,  l=2
%   ...
%   m=4,  l=4
%   m=4,  l=5
%   m=5,  l=5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract number of discrete times (might = 1)
Nt = size(u,2);

% initialize variables for m>=0, l>=2
numCols = (lmax+1)*(lmax+2)/2 - 3;
lvec = zeros(1,numCols);
mvec = zeros(1,numCols);
R_Glm = zeros(Nt,numCols);
R_Clm = zeros(Nt,numCols);
R_Glmm = zeros(Nt,numCols);
R_Clmm = zeros(Nt,numCols);

% loop over m>=0 values
ii = 1;
for m=0:lmax

  for l = max(2,m):lmax
    mvec(ii) = m;
    lvec(ii) = l;

    [R_Glm(:,ii), R_Clm(:,ii)] = response_ifo(l, m, f, u, v, x0);

    % values for m<=0
    R_Glmm(:,ii) = conj(R_Glm(:,ii))*(-1)^m;
    R_Clmm(:,ii) = conj(R_Clm(:,ii))*(-1)^m;

    % increment counter
    ii = ii+1;
  end

end

% extend arrays to have values for negative m
R_Glm = [fliplr(R_Glmm(:,lmax:end)), R_Glm]; 
R_Clm = [fliplr(R_Clmm(:,lmax:end)), R_Clm]; 

% extend vectors of l and m values to negative values of m
lvec = [ fliplr(lvec(lmax:end)), lvec];
mvec = [-fliplr(mvec(lmax:end)), mvec];

% extend vectors by adding on curl modes 
R_Plm.data = [R_Glm, R_Clm];
R_Plm.lvec = [lvec, lvec];
R_Plm.mvec = [mvec, mvec];
R_Plm.lmax = lmax;

return

