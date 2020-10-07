function show_ifos(theta, phi)
%
% add ifo locations to a plot as white circles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

N = length(theta);
for ii=1:N
  % convert theta, phi to lon, lat
  lon = phi(ii)*180/pi;  
  if lon>180
    lon = lon-360;
  end
  lat = (pi/2-theta(ii))*180/pi;  
  [x,y]=m_ll2xy(lon, lat);
  %plot(x,y,'ko'); % black circles
  %plot(x,y,'wo'); % white circles
  %plot(x,y,'ko', 'MarkerSize', 10, 'MarkerFaceColor', [1 1 1]); % white filled circles
end

% allow white circles to be printed  
%set(gcf,'color',[1 1 1]); % Set the figure frame color to white
%set(gca,'color',[1 1 1]); % Set the axis frame color to white
%set(gcf,'InvertHardCopy','off');

return

