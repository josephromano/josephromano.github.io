function [orf, f] = overlap_ifo(site1, site2, lmax, flow, fhigh)
%
% calculate overlap reduction function for an isotropic uncorrelated
% stochastic background
%
% example:
%  overlap_ifo('H1', 'L1', 30, 0, 200)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
const = physConstants('mks');

% get detector info
detector1 = getdetectorNew(site1);
detector2 = getdetectorNew(site2);
u1 = detector1.u;
v1 = detector1.v;
r1 = detector1.r;
u2 = detector2.u;
v2 = detector2.v;
r2 = detector2.r;

% discrete frequencies
Nf = 500;
f = linspace(flow, fhigh, Nf);

% calculate overlap function
orf = zeros(1, Nf);  
for l=2:lmax

  fprintf('working on %d of %d\n', l, lmax);

  for m=-l:l

    [RG1, RC1] = response_ifo(l, m, f, u1, v1, r1);
    [RG2, RC2] = response_ifo(l, m, f, u2, v2, r2);

    orf = orf + RG1.*conj(RG2) + RC1.*conj(RC2);

  end

end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plot
plot(f,real(orf))
xlabel('freq (Hz)')
ylabel('overlap(f)')

outfile = ['overlap_' site1 '_' site2 '_lmax_' num2str(lmax) '.eps'];
print('-depsc2', outfile);

return

