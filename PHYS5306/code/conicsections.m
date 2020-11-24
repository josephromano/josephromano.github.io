% plot conic sections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=1;
mu=1;
k=-1; % choose both signs of k
k=1;

% derived quantities
alpha = l^2/(mu*abs(k));

% minimum energy
if k> 0
  Emin = -0.5*mu*k^2/l^2;
  rmin = l^2/(mu*k);
  fprintf('Emin = %g\n', Emin);
  Evalues = linspace(Emin,-Emin/2, 5);
else
  Evalues = linspace(0.3, 1, 5);
end

% effective potential
%r = linspace(rmin/10, 3*rmin, 100);
%V = -k./r +0.5*l^2./(mu*r.^2);
%figure(1)
%plot(r, V);
%clear r;

% energy values

for ii = 1:length(Evalues)
  E = Evalues(ii);

  % discrete theta values
  theta = linspace(-pi,pi,100);

  epsilon = sqrt(1+(2*E*l^2)/(mu*k^2));

  % conic section
  if k>0
    r = alpha./(1+epsilon*cos(theta));
  else
    r = alpha./(-1+epsilon*cos(theta));
  end
 
  % remove negative r values (correspond to |theta| too large)
  ind = find(r>0);
  r = r(ind);
  theta = theta(ind);

  % cartesian coords
  x = r.*cos(theta);
  y = r.*sin(theta);

  figure(2)
  hold on;
  plot(x,y)

end

axis equal
xlim([-4 10])
ylim([-5 5])

