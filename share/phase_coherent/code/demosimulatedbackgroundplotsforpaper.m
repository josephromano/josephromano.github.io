% script demosimulatedbackgroundplotsforpaper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% point
% make projections of real, imag parts of h+. hx for tot map
makeoneprojectionfromsaveddata('pointsource_inj_pix3072.mat', 'realhp');
clims = get(gca,'Clim');
print('-djpeg', 'point_realhp.jpg');

makeoneprojectionfromsaveddata('pointsource_inj_pix3072.mat', 'imaghp');
set(gca, 'Clim', [-clims(2) clims(2)]);
print('-djpeg', 'point_imaghp.jpg');

makeoneprojectionfromsaveddata('pointsource_inj_pix3072.mat', 'realhc');
clims = get(gca,'Clim');
print('-djpeg', 'point_realhc.jpg');

makeoneprojectionfromsaveddata('pointsource_inj_pix3072.mat', 'imaghc');
set(gca, 'Clim', [-clims(2) clims(2)]);
print('-djpeg', 'point_imaghc.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grad
% make projections of real, imag parts of h+. hx for tot map
makeoneprojectionfromsaveddata('grad_inj_pix3072.mat', 'realhp');
print('-djpeg', 'grad_realhp.jpg');

makeoneprojectionfromsaveddata('grad_inj_pix3072.mat', 'imaghp');
print('-djpeg', 'grad_imaghp.jpg');

makeoneprojectionfromsaveddata('grad_inj_pix3072.mat', 'realhc');
print('-djpeg', 'grad_realhc.jpg');

makeoneprojectionfromsaveddata('grad_inj_pix3072.mat', 'imaghc');
print('-djpeg', 'grad_imaghc.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% curl
% make projections of real, imag parts of h+. hx for tot map
makeoneprojectionfromsaveddata('curl_inj_pix3072.mat', 'realhp');
print('-djpeg', 'curl_realhp.jpg');

makeoneprojectionfromsaveddata('curl_inj_pix3072.mat', 'imaghp');
print('-djpeg', 'curl_imaghp.jpg');

makeoneprojectionfromsaveddata('curl_inj_pix3072.mat', 'realhc');
print('-djpeg', 'curl_realhc.jpg');

makeoneprojectionfromsaveddata('curl_inj_pix3072.mat', 'imaghc');
print('-djpeg', 'curl_imaghc.jpg');

return

