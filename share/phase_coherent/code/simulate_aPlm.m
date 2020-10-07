function a_Plm = simulate_aPlm(type, l0, lmax, seed, theta, phi)
%
% simulate a_Plm for different types of statistically isotropic backgrounds
% 
% Inputs:
%   type:
%     'statiso'  - statistically isotropic with C_l=1 for l=l0, 
%                                               C_l=0 otherwise
%     'statisoG' - same as above but only grad mode
%     'statisoC' - same as above but only curl mode
%     'iso'      - statistically isotropic with C_l=1 for l=2, ..., l0, 
%                                               C_l=0 otherwise
%     'isoG'     - same as above but only grad mode
%     'isoC'     - same as above but only curl mode
%     'pointP'   - deterministic +-polarised point source at theta, phi
%     'pointC'   - deterministic x-polarised point source at theta, phi
%     'point'    - deterministic point source (both +,x) at theta, phi
%
%   l0 - single l-value for 'statiso' type backgrounds
%   lmax - max l-value 
%   seed - seed for random number generation
%   theta, phi - direction to deterministic point source
%  
% Outputs:
%   a_Plm is a structure with fields:
%
%     a_Plm.data = [a_Glm; a_Clm];
%     a_Plm.lvec = lvec;
%     a_Plm.mvec = mvec;
%     a_Plm.lmax = lmax;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state', seed);

% initialize variables for m>=0, l>=2
numRows = (lmax+1)*(lmax+2)/2 - 3;
lvec = zeros(numRows,1);
mvec = zeros(numRows,1);
a_Glm = zeros(numRows,1);
a_Clm = zeros(numRows,1);

% loop over m>=0 values
ii = 1;
for m=0:lmax

  for l = max(2,m):lmax
    mvec(ii) = m;
    lvec(ii) = l;

    switch type

      case 'pointP'
        % first convert to antipodal point 
        % (direction of wave propagation = - direction to source)
        thetaA = pi - theta;
        phiA = pi + phi;
 
        % normalization
        Nl = sqrt(2*factorial(l-2)/factorial(l+2));

        % calculate aGlm, aClm for +-polarised point source
        [Gp, Gm] = Glm(l, m, cos(thetaA));
        Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
        Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
        Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
        a_Glm(ii) =  Nl*Wlmp;
        a_Clm(ii) = -Nl*Xlmp;

      case 'pointC'
        % first convert to antipodal point 
        % (direction of wave propagation = - direction to source)
        thetaA = pi - theta;
        phiA = pi + phi;
 
        % normalization
        Nl = sqrt(2*factorial(l-2)/factorial(l+2));

        % calculate aGlm, aClm for x-polarised point source
        [Gp, Gm] = Glm(l, m, cos(thetaA));
        Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
        Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
        Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
        a_Glm(ii) = Nl*Xlmp;
        a_Clm(ii) = Nl*Wlmp;

      case 'point'
        % first convert to antipodal point 
        % (direction of wave propagation = - direction to source)
        thetaA = pi - theta;
        phiA = pi + phi;
 
        % normalization
        Nl = sqrt(2*factorial(l-2)/factorial(l+2));

        % calculate aGlm, aClm for x-polarised point source
        [Gp, Gm] = Glm(l, m, cos(thetaA));
        Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
        Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
        Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
        a_Glm(ii) = Nl*(Wlmp + Xlmp);
        a_Clm(ii) = Nl*(Wlmp - Xlmp);

      case 'statiso'
        if l==l0
          re = randn;
          im = randn;
          a_Glm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);

          re = randn;
          im = randn;
          a_Clm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      case 'statisoG'
        if l==l0
          re = randn;
          im = randn;
          a_Glm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      case 'statisoC'
        if l==l0
          re = randn;
          im = randn;
          a_Clm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      case 'iso'
        if l<l0 | l==l0
          re = randn;
          im = randn;
          a_Glm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);

          re = randn;
          im = randn;
          a_Clm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      case 'isoG'
        if l<l0 | l==l0
          re = randn;
          im = randn;
          a_Glm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      case 'isoC'
        if l<l0 | l==l0
          re = randn;
          im = randn;
          a_Clm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
        end

      otherwise
        error('unknown type')

    end

    % increment counter
    ii = ii+1;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate values of arrays for negative m
lvecm =  flipud(lvec(lmax:end));
mvecm = -flipud(mvec(lmax:end));

% note that there is no condition relating a_lm and a_l,-m at
% a single frequency f (only at +f and -f)
a_Glmm = zeros(length(lvecm),1);
a_Clmm = zeros(length(lvecm),1);
for ii=1:length(lvecm)
  
  l = lvecm(ii);
  m = mvecm(ii);

  switch type

    case 'pointP'
      % first convert to antipodal point 
      % (direction of wave propagation = - direction to source)
      thetaA = pi - theta;
      phiA = pi + phi;
 
      % normalization
      Nl = sqrt(2*factorial(l-2)/factorial(l+2));

      % calculate aGlm, aClm for +-polarised point source in direction (theta,phi)
      [Gp, Gm] = Glm(l, m, cos(thetaA));
      Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
      Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
      Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
      a_Glmm(ii) =  Nl*Wlmp;
      a_Clmm(ii) = -Nl*Xlmp;

    case 'pointC'
      % first convert to antipodal point 
      % (direction of wave propagation = - direction to source)
      thetaA = pi - theta;
      phiA = pi + phi;
 
      % normalization
      Nl = sqrt(2*factorial(l-2)/factorial(l+2));

      % calculate aGlm, aClm for x-polarised point source in direction (theta,phi)
      [Gp, Gm] = Glm(l, m, cos(thetaA));
      Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
      Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
      Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
      a_Glmm(ii) = Nl*Xlmp;
      a_Clmm(ii) = Nl*Wlmp;

    case 'point'
      % first convert to antipodal point 
      % (direction of wave propagation = - direction to source)
      thetaA = pi - theta;
      phiA = pi + phi;
 
      % normalization
      Nl = sqrt(2*factorial(l-2)/factorial(l+2));

      % calculate aGlm, aClm for x-polarised point source in direction (theta,phi)
      [Gp, Gm] = Glm(l, m, cos(thetaA));
      Nlm = sqrt( ((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)) );    
      Wlmp = conj(2*Nlm*Gp*exp(i*m*phiA));
      Xlmp = conj(2*sqrt(-1)*Nlm*Gm*exp(i*m*phiA));
      a_Glmm(ii) = Nl*(Wlmp + Xlmp);
      a_Clmm(ii) = Nl*(Wlmp - Xlmp);

    case 'statiso'
      if l==l0
        re = randn;
        im = randn;
        a_Glmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);

        re = randn;
        im = randn;
        a_Clmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    case 'statisoG'
      if l==l0
        re = randn;
        im = randn;
        a_Glmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    case 'statisoC'
      if l==l0
        re = randn;
        im = randn;
        a_Clmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    case 'iso'
      if l<l0 | l==l0
        re = randn;
        im = randn;
        a_Glmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);

        re = randn;
        im = randn;
        a_Clmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    case 'isoG'
      if l<l0 | l==l0
        re = randn;
        im = randn;
        a_Glmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    case 'isoC'
      if l<l0 | l==l0
        re = randn;
        im = randn;
        a_Clmm(ii) = (1/sqrt(2))*(re + sqrt(-1)*im);
      end

    otherwise
      error('unknown type')

  end
end

% extend arrays to have values for negative m
a_Glm = [a_Glmm; a_Glm]; 
a_Clm = [a_Clmm; a_Clm]; 

% extend vectors of l and m values to negative values of m
lvec = [lvecm; lvec];
mvec = [mvecm; mvec];

% combine grad and curl modes
a_Plm.data = [a_Glm; a_Clm];
a_Plm.lvec = [lvec; lvec];
a_Plm.mvec = [mvec; mvec];
a_Plm.lmax = lmax;

% debugging
%a_Plm.mvec
%a_Plm.lvec
%a_Plm.data

% save data for later use
filename = [type '_sph' num2str(lmax) '.mat'];
save(filename, 'a_Plm')

return

