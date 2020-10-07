setup

% generate basis skies (could for all types of recoveries)
clear all
basis_skies_pixels_ifo_noise(100, 10, 6, 400, 'both');
clear all
basis_skies_pixels_ifo_noise(100, 10, 3, 400, 'both');
clear all
basis_skies_pixels_ifo_noise(100, 10, 3, 800, 'both');

% recover
clear all
recover_skies_pixels_ifo_noise('pointsource', 'both', 'response', 768, 6, 400);
clear all
recover_skies_pixels_ifo_noise('grad', 'both', 'response', 768, 6, 400);
clear all
recover_skies_pixels_ifo_noise('curl', 'both', 'response', 768, 6, 400);

clear all
recover_skies_pixels_ifo_noise('pointsource', 'both', 'response', 768, 3, 400);
clear all
recover_skies_pixels_ifo_noise('grad', 'both', 'response', 768, 3, 400);
clear all
recover_skies_pixels_ifo_noise('curl', 'both', 'response', 768, 3, 400);

clear all
recover_skies_pixels_ifo_noise('pointsource', 'both', 'response', 768, 3, 800);
clear all
recover_skies_pixels_ifo_noise('grad', 'both', 'response', 768, 3, 800);
clear all
recover_skies_pixels_ifo_noise('curl', 'both', 'response', 768, 3, 800);

