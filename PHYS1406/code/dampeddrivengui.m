function dampeddriven()
%
% integrate F = ma for damped, driven harmonic motion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%define figure
f=figure('Position',[300 400 900 440],...
         'Name','Damped and driven oscillations');

%define sliders     
% note 'Position' [left bottom width height]

%%%%%%%%%%%%%%%%%
dampingslider=uicontrol(f,'Style','slider',...
             'Max',2.5,'Min',0,'Value',0,...
             'SliderStep',[1/50 1/50],...
             'Position',[330 45 20 100],...
             'Callback',@dampingslider_callback);
dampingvalue=uicontrol(f,'Style','text',...
             'String',num2str(get(dampingslider,'Value')),...
             'Position',[300 150 80 15],...
             'Callback',@dampingvalue_callback);
dampingtext=uicontrol(f,'Style','text','String',...
             'Damping','Position',[300 170 80 15]);

amplitudeslider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[1/10 1/10],...
             'Position',[430 45 20 100],...
             'Callback',@amplitudeslider_callback);
amplitudevalue=uicontrol(f,'Style','text',...
             'String',num2str(get(amplitudeslider,'Value')),...
             'Position',[400 150 80 15],...
             'Callback',@amplitudevalue_callback);
amplitudetext=uicontrol(f,'Style','text','String',...
             'Amplitude','Position',[400 170 80 15]);

frequencyslider=uicontrol(f,'Style','slider',...
             'Max',3,'Min',0,'Value',0,...
             'SliderStep', [1/60 1/60],...
             'Position',[530 45 20 100],...
             'Callback',@frequencyslider_callback);
frequencyvalue=uicontrol(f,'Style','text',...
             'String',num2str(get(frequencyslider,'Value')),...
             'Position',[500 150 80 15],...
             'Callback',@frequencyvalue_callback);
frequencytext=uicontrol(f,'Style','text','String',...
             'Frequency','Position',[500 170 80 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some parameters
x0 = 1; % max displacement (m)
f0 = 1; % natural frequency (Hz)
omega0 = 2*pi*f0; 

% initial values
damping = 0;
A = 0;
f_ratio = 0;

% convert slider values to relevant parameter for differential equation
beta = damping*omega0; % beta = b/2m (Hz); F = -bv (damping force)
F0overm = A*omega0^2*x0; % driving force F0 over m (m/s^2)
f = f_ratio*f0;
omega = 2*pi*f;

% display resonant frequency and maximum displacement on resonance
if omega0^2>2*beta^2
  omegaR = sqrt(omega0^2 - 2*beta^2); % resonant frequency
  %DR = F0overm/(2*beta*sqrt(omega0^2 - beta^2)); 
  DR = F0overm/sqrt((omega0^2 - omegaR^2)^2 + 4*beta^2*omegaR^2);
  fprintf('resonant frequency = %f\n', omegaR/(2*pi));
  fprintf('steady state max displacement on resonance = %f\n', DR);
else
  fprintf('resonant frequency does not exist; damping too large\n')
end

% maximum displacement
D = F0overm/sqrt((omega0^2 - omega^2)^2 + 4*beta^2*omega^2);
fprintf('steady state max displacement = %f\n', D);

% discrete times
tmin = 0;
tmax = 25; % sec
Fs = 100*max(f0, f); % sample frequency (Hz)
dt = 1/Fs;
N = floor(tmax/dt);
t = linspace(tmin, tmax, N);

% driving motion
xdriving = A*cos(omega*t);

% initialize variables
x = zeros(N,1);
v = zeros(N,1);
x(1) = x0;
v(1) = 0;

% numerically integrate F = mxddot
for ii=2:N

  % evaluate acceleration from force
  a = F0overm*cos(omega*t(ii)) - 2*beta*v(ii-1) - omega0^2*x(ii-1);
  v(ii) = v(ii-1) + a*dt;
  x(ii) = x(ii-1) + v(ii)*dt;

end

%%%%%%%%%%%%%%%
%define plot (with labels and limits)
displacementplot=axes('Units','Pixels','Position',[80 250 750 160]);
plot(t,x,'-b',t,xdriving,'-r');
grid on
xmax = max(abs(x));
ylim([-xmax xmax])
xlabel('time (seconds)')
ylabel('displacement')
legend('mass motion', 'driving motion', 'location', 'northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% damping contribution
function dampingslider_callback(hObject,eventdata)

  delete(displacementplot)

  val = get(hObject,'Value');
  set(dampingvalue,'String',num2str(val));
  damping = get(dampingslider,'Value');

  beta = damping*omega0; % beta = b/2m (Hz); F = -bv (damping force)

  % display resonant frequency and maximum displacement on resonance
  if omega0^2>2*beta^2
    omegaR = sqrt(omega0^2 - 2*beta^2); % resonant frequency
    %DR = F0overm/(2*beta*sqrt(omega0^2 - beta^2)); 
    DR = F0overm/sqrt((omega0^2 - omegaR^2)^2 + 4*beta^2*omegaR^2);
    fprintf('resonant frequency = %f\n', omegaR/(2*pi));
    fprintf('steady state max displacement on resonance = %f\n', DR);
  else
    fprintf('resonant frequency does not exist; damping too large\n')
  end

  % maximum displacement
  D = F0overm/sqrt((omega0^2 - omega^2)^2 + 4*beta^2*omega^2);
  fprintf('steady state max displacement = %f\n', D);

  % discrete times
  tmin = 0;
  tmax = 25; % sec
  Fs = 100*max(f0, f); % sample frequency (Hz)
  dt = 1/Fs;
  N = floor(tmax/dt);
  t = linspace(tmin, tmax, N);

  % driving motion
  xdriving = A*cos(omega*t);

  % initialize variables
  x = zeros(N,1);
  v = zeros(N,1);
  x(1) = x0;
  v(1) = 0;

  % numerically integrate F = mxddot
  for ii=2:N

    % evaluate acceleration from force
    a = F0overm*cos(omega*t(ii)) - 2*beta*v(ii-1) - omega0^2*x(ii-1);
    v(ii) = v(ii-1) + a*dt;
    x(ii) = x(ii-1) + v(ii)*dt;

  end

  %define plot (with labels and limits)
  displacementplot=axes('Units','Pixels','Position',[80 250 750 160]);
  %set(ff,'CurrentAxes',displacementplot)
  plot(t,x,'-b',t,xdriving,'-r');
  grid on
  xmax = max(abs(x));
  ylim([-xmax xmax])
  xlabel('time (seconds)')
  ylabel('displacement')
  legend('mass motion', 'driving motion', 'location', 'northeast')

end 
     
function dampingvalue_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude contribution 
function amplitudeslider_callback(hObject,eventdata)

  delete(displacementplot)

  val = get(hObject,'Value');
  set(amplitudevalue,'String',num2str(val));
  A = get(amplitudeslider,'Value');

  F0overm = A*omega0^2*x0; % driving force F0 over m (m/s^2)

  % display resonant frequency and maximum displacement on resonance
  if omega0^2>2*beta^2
    omegaR = sqrt(omega0^2 - 2*beta^2); % resonant frequency
    %DR = F0overm/(2*beta*sqrt(omega0^2 - beta^2));
    DR = F0overm/sqrt((omega0^2 - omegaR^2)^2 + 4*beta^2*omegaR^2);
    fprintf('resonant frequency = %f\n', omegaR/(2*pi));
    fprintf('steady state max displacement on resonance = %f\n', DR);
  else
    fprintf('resonant frequency does not exist; damping too large\n')
  end

  % maximum displacement
  D = F0overm/sqrt((omega0^2 - omega^2)^2 + 4*beta^2*omega^2);
  fprintf('steady state max displacement = %f\n', D);

  % discrete times
  tmin = 0;
  tmax = 25; % sec
  Fs = 100*max(f0, f); % sample frequency (Hz)
  dt = 1/Fs;
  N = floor(tmax/dt);
  t = linspace(tmin, tmax, N);

  % driving motion
  xdriving = A*cos(omega*t);

  % initialize variables
  x = zeros(N,1);
  v = zeros(N,1);
  x(1) = x0;
  v(1) = 0;

  % numerically integrate F = mxddot
  for ii=2:N

    % evaluate acceleration from force
    a = F0overm*cos(omega*t(ii)) - 2*beta*v(ii-1) - omega0^2*x(ii-1);
    v(ii) = v(ii-1) + a*dt;
    x(ii) = x(ii-1) + v(ii)*dt;

  end

  %define plot (with labels and limits)
  displacementplot=axes('Units','Pixels','Position',[80 250 750 160]);
  %set(ff,'CurrentAxes',displacementplot)
  plot(t,x,'-b',t,xdriving,'-r');
  grid on
  xmax = max(abs(x));
  ylim([-xmax xmax])
  xlabel('time (seconds)')
  ylabel('displacement')
  legend('mass motion', 'driving motion', 'location', 'northeast')

end 
     
function amplitudevalue_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency contribution
function frequencyslider_callback(hObject,eventdata)

  delete(displacementplot)

  val = get(hObject,'Value');
  set(frequencyvalue,'String',num2str(val));
  f_ratio = get(frequencyslider,'Value');

  f = f_ratio*f0;
  omega = 2*pi*f;

  % display resonant frequency and maximum displacement on resonance
  if omega0^2>2*beta^2
    omegaR = sqrt(omega0^2 - 2*beta^2); % resonant frequency
    %DR = F0overm/(2*beta*sqrt(omega0^2 - beta^2));
    DR = F0overm/sqrt((omega0^2 - omegaR^2)^2 + 4*beta^2*omegaR^2);
    fprintf('resonant frequency = %f\n', omegaR/(2*pi));
    fprintf('steady state max displacement on resonance = %f\n', DR);
  else
    fprintf('resonant frequency does not exist; damping too large\n')
  end

  % maximum displacement
  D = F0overm/sqrt((omega0^2 - omega^2)^2 + 4*beta^2*omega^2);
  fprintf('steady state max displacement = %f\n', D);

  % discrete times
  tmin = 0;
  tmax = 25; % sec
  Fs = 100*max(f0, f); % sample frequency (Hz)
  dt = 1/Fs;
  N = floor(tmax/dt);
  t = linspace(tmin, tmax, N);

  % driving motion
  xdriving = A*cos(omega*t);

  % initialize variables
  x = zeros(N,1);
  v = zeros(N,1);
  x(1) = x0;
  v(1) = 0;

  % numerically integrate F = mxddot
  for ii=2:N

    % evaluate acceleration from force
    a = F0overm*cos(omega*t(ii)) - 2*beta*v(ii-1) - omega0^2*x(ii-1);
    v(ii) = v(ii-1) + a*dt;
    x(ii) = x(ii-1) + v(ii)*dt;

  end

  %define plot (with labels and limits)
  displacementplot=axes('Units','Pixels','Position',[80 250 750 160]);
  %set(ff,'CurrentAxes',displacementplot)
  plot(t,x,'-b',t,xdriving,'-r');
  grid on
  xmax = max(abs(x));
  ylim([-xmax xmax])
  xlabel('time (seconds)')
  ylabel('displacement')
  legend('mass motion', 'driving motion', 'location', 'northeast')

end 
     
function frequencyvalue_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
end
