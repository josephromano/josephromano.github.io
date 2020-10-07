function [Cn, invCn] = calNoiseCovariance(Nd, Nt, deltaf);
%
% calculate noise covariance matrix for Nd detectors, Nt times
% deltaf is the size of the frequency bins (typically 0.25 Hz)
%
%###################################################################
% Assuming Hanford, Livingston, India, and Australia have same PSD.
% Sn(f=100Hz) = 1.5907433033016117e-47 Hz^{-1} for aLIGO
% Sn(f=100Hz) = 2.0630282699660411e-47 Hz^{-1} for AdVirgo
% Sn(f=100Hz) = 9.3200152368999989e-48 Hz^{-1} for KAGRA
% Assuming \delta f = 1.0 Hz
%###################################################################

if Nd==1
  % Cn = ligo
  Cn_diag = 1.591e-47*ones(Nt,1);
  Cn_diag = Cn_diag/(4*deltaf); % normalization for 1-sided spectra

elseif Nd==3
  % Cn = [ligo,ligo,virgo]
  Cn_diag = [1.591e-47*ones(Nt,1); ...
             1.591e-47*ones(Nt,1); ...
             2.063e-47*ones(Nt,1)];
  Cn_diag = Cn_diag/(4*deltaf); % normalization for 1-sided spectra

elseif Nd==6
  % Cn = [ligo,ligo,virgo,kagra,ligo,ligo]

  Cn_diag = [1.591e-47*ones(Nt,1); ...
             1.591e-47*ones(Nt,1); ...
             2.063e-47*ones(Nt,1); ...
             9.320e-48*ones(Nt,1); ...
             1.591e-47*ones(Nt,1); ...
             1.591e-47*ones(Nt,1)];
  Cn_diag = Cn_diag/(4*deltaf); % normalization for 1-sided spectra

else
  error('not implemented for N~=6 detectors!');
end

% calculate covariance and inverse covariance matrices
Cn = diag(Cn_diag);
invCn_diag = 1./Cn_diag;
invCn = diag(invCn_diag);

return

