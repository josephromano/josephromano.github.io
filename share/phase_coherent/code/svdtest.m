function svdtest(type)
%
% compare ML estimations for calculation schemes
%
% type = 'over' or 'under' for over or underdetermined system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
  case 'over'

    % N > M (overdetermined) %%%%%%%%%%%%
    N = 3; M = 2;
    fprintf('overdetermined: N=%d, M=%d\n', N, M)

    % covariance matrix, reponse matrix, gw sky
    C = [1 0 0; 0 2 0; 0 0 3];
    C = [1 0.5 0.25; 0.5 2 0.5; 0.25 0.5 3]
    R = [1 0; 0 1; 1 1]
    %R = [1 2; 3 4; 5 6]
    h = [10; 20];
    iC = inv(C);
    L = chol(iC, 'lower');
    S = L'*R

    % method 1 
    fprintf('method 1:\n')
    %F = R'*iC*R
    F = S'*S
    iF = inv(F)
    iF * F

    hest = iF * F * h

    % method 2 %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('method 2:\n')
    [U, Sigma, V] = svd(R);

    % construct pseudo inverse of Sigma
    pinvSigma = zeros(size(Sigma'));
    K = min(size(Sigma,1), size(Sigma,2));
    if K==1
      pinvSigma(1,1)=1/Sigma(1,1);
    else
    svals = diag(Sigma);
    pinvSigma(1:K,1:K) = diag(1./svals);
    end

    Sigma
    pinvSigma
    pinvSigma*Sigma

    R
    pinvR = V*pinvSigma*U'
    pinvR*R

    hest = pinvR * R * h
    
    % EXTRA BELOW!!
    % construct pseudo inverse of Sigmabar
    [Ubar, Sigmabar, Vbar] = svd(S);

    pinvSigmabar = zeros(size(Sigmabar'));
    K = min(size(Sigmabar,1), size(Sigmabar,2));
    if K==1
      pinvSigmabar(1,1)=1/Sigmabar(1,1);
    else
      svals = diag(Sigmabar);
      pinvSigmabar(1:K,1:K) = diag(1./svals);
    end

    L'
    R
    L * L'
    R * R'
    pinvR = V*pinvSigma*U'
    pinvS = Vbar*pinvSigmabar*Ubar'
    pinvR * inv(L')

  case 'under'

    % N < M (underdetermined) %%%%%%%%%%%%
    N = 2; M = 3;
    fprintf('underdetermined: N=%d, M=%d\n', N, M)

    % covariance matrix, reponse matrix, gw sky
    C = [1 0; 0 2];
    C = [1 0.5; 0.5 2]
    R = [1 0 0; 0 1 0]
    %R = [1 2 3; 4 5 6]
    h = [10; 20; 30];
    iC = inv(C);
    L = chol(iC, 'lower');
    S = L'*R

    % method 1 
    fprintf('method 1:\n')
    %F = R'*iC*R
    F = S'*S

    % inv(F) does not exist
    % need to find pseudo inverse of F using SVD of S
    [Ubar, Sigmabar, Vbar] = svd(S);

    % construct pseudo inverse of Sigmabar
    pinvSigmabar = zeros(size(Sigmabar'));
    K = min(size(Sigmabar,1), size(Sigmabar,2));
    if K==1
      pinvSigmabar(1,1)=1/Sigmabar(1,1);
    else
      svals = diag(Sigmabar);
      pinvSigmabar(1:K,1:K) = diag(1./svals);
    end

    Sigmabar
    pinvSigmabar
    pinvSigmabar*Sigmabar

    % now construct pseudo inverse of F
    pinvF = Vbar*pinvSigmabar*pinvSigmabar'*Vbar'
    pinvF * F 
    hest = pinvF * F * h

    % method 2 %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('method 2:\n')
    [U, Sigma, V] = svd(R);

    % construct pseudo inverse of Sigma
    pinvSigma = zeros(size(Sigma'));
    K = min(size(Sigma,1), size(Sigma,2));
    if K==1
      pinvSigma(1,1)=1/Sigma(1,1);
    else
    svals = diag(Sigma);
    pinvSigma(1:K,1:K) = diag(1./svals);
    end

    Sigma
    pinvSigma
    pinvSigma*Sigma

    R
    pinvR = V*pinvSigma*U'
    pinvR*R

    hest = pinvR * R * h

    % EXTRA BELOW!!
    L * L'
    R * R'
    pinvR = V*pinvSigma*U'
    pinvS = Vbar*pinvSigmabar*Ubar'
    pinvR * inv(L')

  otherwise
    error('unrecognized type\n');

end

return

