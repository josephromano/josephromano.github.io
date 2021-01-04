function amplitudemodulation()
%
% combine carrier and signal waves using amplitude
% modulation
%
% y(t) = (1 + x_m(t)) * x_c(t)
% 
% where 
%
% x_m(t) = A_m sin(2pi f_m t + phi_m) 
% x_c(t) = A_c sin(2pi f_c t + phi_c) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%define figure
f=figure('Position',[300 400 900 640],...
         'Name','Amplitude modulation');

%define sliders     
% note 'Position' [left bottom width height]

%%%%%%%%%%%%%%%%%
amp1slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[130 45 20 100],...
             'Callback',@amp1slider_callback);
amp1value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp1slider,'Value')),...
             'Position',[100 150 80 15],...
             'Callback',@amp1value_callback);
amp1text=uicontrol(f,'Style','text','String',...
             'Amplitude_c','Position',[100 170 80 15]);

freq1slider=uicontrol(f,'Style','slider',...
             'Max',10.0,'Min',0,'Value',1,...
             'SliderStep',[1/10 1/10],...
             'Position',[230 45 20 100],...
             'Callback',@freq1slider_callback);
freq1value=uicontrol(f,'Style','text',...
             'String',num2str(get(freq1slider,'Value')),...
             'Position',[200 150 80 15],...
             'Callback',@freq1value_callback);
freq1text=uicontrol(f,'Style','text','String',...
             'Freq_c (Hz)','Position',[200 170 80 15]);

phase1slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[330 45 20 100],...
             'Callback',@phase1slider_callback);
phase1value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase1slider,'Value')),...
             'Position',[300 150 80 15],...
             'Callback',@phase1value_callback);
phase1text=uicontrol(f,'Style','text','String',...
             'Phase_c','Position',[300 170 80 15]);

%%%%%%%%%%%%%%%%%
amp2slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[530 45 20 100],...
             'Callback',@amp2slider_callback);
amp2value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp2slider,'Value')),...
             'Position',[500 150 80 15],...
             'Callback',@amp2value_callback);
amp2text=uicontrol(f,'Style','text','String',...
             'Amplitude_m','Position',[500 170 80 15]);

freq2slider=uicontrol(f,'Style','slider',...
             'Max',10.0,'Min',0,'Value',1,...
             'SliderStep',[1/10 1/10],...
             'Position',[630 45 20 100],...
             'Callback',@freq2slider_callback);
freq2value=uicontrol(f,'Style','text',...
             'String',num2str(get(freq2slider,'Value')),...
             'Position',[600 150 80 15],...
             'Callback',@freq2value_callback);
freq2text=uicontrol(f,'Style','text','String',...
             'Freq_m (Hz)','Position',[600 170 80 15]);

phase2slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[730 45 20 100],...
             'Callback',@phase2slider_callback);
phase2value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase2slider,'Value')),...
             'Position',[700 150 80 15],...
             'Callback',@phase2value_callback);
phase2text=uicontrol(f,'Style','text','String',...
             'Phase_m','Position',[700 170 80 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define plots (with labels and limits)

totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
grid on
xlimmin = 0;
xlimmax = 3;
ylimmin = -3;
ylimmax = 3;
xlim([xlimmin xlimmax])
ylim([ylimmin ylimmax])
%xlabel('time (seconds)')
ylabel('summed waveform')

indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
grid on
xlimmin = 0;
xlimmax = 3;
ylimmin = -3;
ylimmax = 3;
xlim([xlimmin xlimmax])
ylim([ylimmin ylimmax])
xlabel('time (seconds)')
ylabel('individual waveforms')

%discrete times
tmin = 0;
tmax = 4; % sec
N = 1000;
t = linspace(tmin, tmax, N);

%initial values
a1 = 0;
a2 = 0;
f1 = 1;
f2 = 1;
phi1 = 0;
phi2 = 0;

% waveforms
y = zeros(1,N);
y1 = zeros(1,N);
y2 = zeros(1,N);

%%%%%%%%%%%%%%%%
% amplitude contribution from 1st component
function amp1slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(amp1value,'String',num2str(val));
  a1 = get(amp1slider,'Value');
  y1 = a1*sin(2*pi*f1*t + phi1*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function amp1value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% frequency contribution from 1st component
function freq1slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(freq1value,'String',num2str(val));
  f1 = get(freq1slider,'Value');
  y1 = a1*sin(2*pi*f1*t + phi1*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function freq1value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 1st component
function phase1slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(phase1value,'String',num2str(val));
  phi1 = get(phase1slider,'Value');
  y1 = a1*sin(2*pi*f1*t + phi1*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function phase1value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 2nd component
function amp2slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(amp2value,'String',num2str(val));
  a2 = get(amp2slider,'Value');
  y2 = a2*sin(2*pi*f2*t + phi2*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function amp2value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% frequency contribution from 2nd component
function freq2slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(freq2value,'String',num2str(val));
  f2 = get(freq2slider,'Value');
  y2 = a2*sin(2*pi*f2*t + phi2*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function freq2value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 2nd component
function phase2slider_callback(hObject,eventdata)
  delete(totwaveformplot)

  val = get(hObject,'Value');
  set(phase2value,'String',num2str(val));
  phi2 = get(phase2slider,'Value');
  y2 = a2*sin(2*pi*f2*t + phi2*pi/180);
  y = (1+y2).*y1;

  %define plot (with labels and limits)
  totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
  set(f,'CurrentAxes',totwaveformplot)
  plot(t,y,'-r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('modulated waveform')
  legend('(1+ym)*yc', 'location', 'northeast')

  indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
  set(f,'CurrentAxes',indwaveformplot)
  plot(t,y1,'-b',t,y2,'-k');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('individual waveforms')
  legend('yc', 'ym', 'location', 'northeast')

end 
     
function phase2value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
end
