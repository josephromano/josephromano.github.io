function [t_hPix] = translate_map(hPix, x0, f)
%
% translate a pixel-based map hPix = [hpPix, hcPix] (in row or
% column vector form) by a translation vector x0 (3x1), where f 
% is frequency in Hz
%
% the translated map is given by
%
% t_hpPix = hpPix exp(-i2pi f k.x0/c)
% t_hcPix = hcPix exp(-i2pi f k.x0/c)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const = physConstants('mks');
c = const.c;

% extract number of pixels
nPix = length(hPix)/2;
nSide = nPix2nSide(nPix);

% extract h+, hx
hpPix = hPix(1:nPix);
hcPix = hPix(nPix+1:end);

% calculate cell array of theta, phi values for each pixel
tp = pix2ang(nSide);

% loop over pixels on sphere
t_hpPix = zeros(size(hpPix));
t_hcPix = zeros(size(hcPix));

for ii = 1:1:nPix

  %if mod(ii,50)==0
  %  fprintf('working on %d of %d\n', ii, nPix)
  %end

  % extract theta, phi (direction to source)
  theta = tp{ii}(1);
  phi = tp{ii}(2);

  % first convert to antipodal point (khat = - direction to source)
  theta = pi - theta;
  phi = pi + phi;

  % calculate phase
  k = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
  phase = -sqrt(-1)*2*pi*f*k*x0/c;

  % calculate translated h+, hx
  t_hp = hpPix(ii)*exp(phase);
  t_hc = hcPix(ii)*exp(phase);
  t_hpPix(ii) = t_hp;
  t_hcPix(ii) = t_hc;

end

% package output
if size(t_hpPix,1)<size(t_hpPix,2)
  t_hPix = [t_hpPix, t_hcPix];
else
  t_hPix = [t_hpPix; t_hcPix];
end
  
return

