function fouriersynthesize()
%
% fourier synthesizer using first eight harmonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%define figure
f=figure('Position',[300 400 950 640],...
         'Name','Fourier synthesizer');

%define sliders     
% note 'Position' [left bottom width height]

%%%%%%%%%%%%%%%%%
amp1slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[130 245 20 100],...
             'Callback',@amp1slider_callback);
amp1value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp1slider,'Value')),...
             'Position',[100 350 80 15],...
             'Callback',@amp1value_callback);
amp1text=uicontrol(f,'Style','text','String',...
             'Amplitude 1','Position',[100 370 80 15]);

phase1slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[130 75 20 100],...
             'Callback',@phase1slider_callback);
phase1value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase1slider,'Value')),...
             'Position',[100 180 80 15],...
             'Callback',@phase1value_callback);
phase1text=uicontrol(f,'Style','text','String',...
             'Phase 1','Position',[100 200 80 15]);

%%%%%%%%%%%%%%%%%
amp2slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[230 245 20 100],...
             'Callback',@amp2slider_callback);
amp2value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp1slider,'Value')),...
             'Position',[200 350 80 15],...
             'Callback',@amp2value_callback);
amp2text=uicontrol(f,'Style','text','String',...
             'Amplitude 2','Position',[200 370 80 15]);

phase2slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[230 75 20 100],...
             'Callback',@phase2slider_callback);
phase2value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase2slider,'Value')),...
             'Position',[200 180 80 15],...
             'Callback',@phase2value_callback);
phase2text=uicontrol(f,'Style','text','String',...
             'Phase 2','Position',[200 200 80 15]);

%%%%%%%%%%%%%%%%%
amp3slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[330 245 20 100],...
             'Callback',@amp3slider_callback);
amp3value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp3slider,'Value')),...
             'Position',[300 350 80 15],...
             'Callback',@amp3value_callback);
amp3text=uicontrol(f,'Style','text','String',...
             'Amplitude 3','Position',[300 370 80 15]);

phase3slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[330 75 20 100],...
             'Callback',@phase3slider_callback);
phase3value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase3slider,'Value')),...
             'Position',[300 180 80 15],...
             'Callback',@phase3value_callback);
phase3text=uicontrol(f,'Style','text','String',...
             'Phase 3','Position',[300 200 80 15]);

%%%%%%%%%%%%%%%%%
amp4slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[430 245 20 100],...
             'Callback',@amp4slider_callback);
amp4value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp4slider,'Value')),...
             'Position',[400 350 80 15],...
             'Callback',@amp4value_callback);
amp4text=uicontrol(f,'Style','text','String',...
             'Amplitude 4','Position',[400 370 80 15]);

phase4slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[430 75 20 100],...
             'Callback',@phase4slider_callback);
phase4value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase4slider,'Value')),...
             'Position',[400 180 80 15],...
             'Callback',@phase4value_callback);
phase4text=uicontrol(f,'Style','text','String',...
             'Phase 4','Position',[400 200 80 15]);

%%%%%%%%%%%%%%%%%
amp5slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[530 245 20 100],...
             'Callback',@amp5slider_callback);
amp5value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp5slider,'Value')),...
             'Position',[500 350 80 15],...
             'Callback',@amp5value_callback);
amp5text=uicontrol(f,'Style','text','String',...
             'Amplitude 5','Position',[500 370 80 15]);

phase5slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[530 75 20 100],...
             'Callback',@phase5slider_callback);
phase5value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase5slider,'Value')),...
             'Position',[500 180 80 15],...
             'Callback',@phase5value_callback);
phase5text=uicontrol(f,'Style','text','String',...
             'Phase 5','Position',[500 200 80 15]);

%%%%%%%%%%%%%%%%%
amp6slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[630 245 20 100],...
             'Callback',@amp6slider_callback);
amp6value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp6slider,'Value')),...
             'Position',[600 350 80 15],...
             'Callback',@amp6value_callback);
amp6text=uicontrol(f,'Style','text','String',...
             'Amplitude 6','Position',[600 370 80 15]);

phase6slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[630 75 20 100],...
             'Callback',@phase6slider_callback);
phase6value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase6slider,'Value')),...
             'Position',[600 180 80 15],...
             'Callback',@phase6value_callback);
phase6text=uicontrol(f,'Style','text','String',...
             'Phase 6','Position',[600 200 80 15]);

%%%%%%%%%%%%%%%%%
amp7slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[730 245 20 100],...
             'Callback',@amp7slider_callback);
amp7value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp7slider,'Value')),...
             'Position',[700 350 80 15],...
             'Callback',@amp7value_callback);
amp7text=uicontrol(f,'Style','text','String',...
             'Amplitude 7','Position',[700 370 80 15]);

phase7slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[730 75 20 100],...
             'Callback',@phase7slider_callback);
phase7value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase7slider,'Value')),...
             'Position',[700 180 80 15],...
             'Callback',@phase7value_callback);
phase7text=uicontrol(f,'Style','text','String',...
             'Phase 7','Position',[700 200 80 15]);

%%%%%%%%%%%%%%%%%
amp8slider=uicontrol(f,'Style','slider',...
             'Max',1.0,'Min',0,'Value',0,...
             'SliderStep',[0.1 0.1],...
             'Position',[830 245 20 100],...
             'Callback',@amp8slider_callback);
amp8value=uicontrol(f,'Style','text',...
             'String',num2str(get(amp8slider,'Value')),...
             'Position',[800 350 80 15],...
             'Callback',@amp8value_callback);
amp8text=uicontrol(f,'Style','text','String',...
             'Amplitude 8','Position',[800 370 80 15]);

phase8slider=uicontrol(f,'Style','slider',...
             'Max',180.0,'Min',-180.0,'Value',0,...
             'SliderStep', [1/36 1/36],...
             'Position',[830 75 20 100],...
             'Callback',@phase8slider_callback);
phase8value=uicontrol(f,'Style','text',...
             'String',num2str(get(phase8slider,'Value')),...
             'Position',[800 180 80 15],...
             'Callback',@phase8value_callback);
phase8text=uicontrol(f,'Style','text','String',...
             'Phase 8','Position',[800 200 80 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define plot (with labels and limits)
waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
grid on
xlimmin = 0;
xlimmax = 4;
ylimmin = -3;
ylimmax = 3;
xlim([xlimmin xlimmax])
ylim([ylimmin ylimmax])
xlabel('time (seconds)')
ylabel('summed waveform')

%discrete times
tmin = 0;
tmax = 4; % sec
N = 1000;
t = linspace(tmin, tmax, N);

%harmonics
f1 = 1; % Hz
f2 = 2*f1;
f3 = 3*f1;
f4 = 4*f1;
f5 = 5*f1;
f6 = 6*f1;
f7 = 7*f1;
f8 = 8*f1;

%initial amplitudes
a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
a8 = 0;

%initial phase offsets
phi1 = 0;
phi2 = 0;
phi3 = 0;
phi4 = 0;
phi5 = 0;
phi6 = 0;
phi7 = 0;
phi8 = 0;

% waveforms
y = zeros(1,N);
y1 = zeros(1,N);
y2 = zeros(1,N);
y3 = zeros(1,N);
y4 = zeros(1,N);
y5 = zeros(1,N);
y6 = zeros(1,N);
y7 = zeros(1,N);
y8 = zeros(1,N);

%%%%%%%%%%%%%%%%
% amplitude contribution from 1st harmonic
function amp1slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp1value,'String',num2str(val));
  a1 = get(amp1slider,'Value');
  y1 = a1*sin(2*pi*f1*t + phi1*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp1value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 1st harmonic
function phase1slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase1value,'String',num2str(val));
  phi1 = get(phase1slider,'Value');
  y1 = a1*sin(2*pi*f1*t + phi1*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase1value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 2nd harmonic
function amp2slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp2value,'String',num2str(val));
  a2 = get(amp2slider,'Value');
  y2 = a2*sin(2*pi*f2*t + phi2*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp2value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 2nd harmonic
function phase2slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase2value,'String',num2str(val));
  phi2 = get(phase2slider,'Value');
  y2 = a2*sin(2*pi*f2*t + phi2*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase2value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 3rd harmonic
function amp3slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp3value,'String',num2str(val));
  a3 = get(amp3slider,'Value');
  y3 = a3*sin(2*pi*f3*t + phi3*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp3value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 3rd harmonic
function phase3slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase3value,'String',num2str(val));
  phi3 = get(phase3slider,'Value');
  y3 = a3*sin(2*pi*f3*t + phi3*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase3value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 4th harmonic
function amp4slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp4value,'String',num2str(val));
  a4 = get(amp4slider,'Value');
  y4 = a4*sin(2*pi*f4*t + phi4*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp4value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 4th harmonic
function phase4slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase4value,'String',num2str(val));
  phi4 = get(phase4slider,'Value');
  y4 = a4*sin(2*pi*f4*t + phi4*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase4value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 5th harmonic
function amp5slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp5value,'String',num2str(val));
  a5 = get(amp5slider,'Value');
  y5 = a5*sin(2*pi*f5*t + phi5*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp5value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 5th harmonic
function phase5slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase5value,'String',num2str(val));
  phi5 = get(phase5slider,'Value');
  y5 = a5*sin(2*pi*f5*t + phi5*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase5value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 6th harmonic
function amp6slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp6value,'String',num2str(val));
  a6 = get(amp6slider,'Value');
  y6 = a6*sin(2*pi*f6*t + phi6*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp6value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 6th harmonic
function phase6slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase6value,'String',num2str(val));
  phi6 = get(phase6slider,'Value');
  y6 = a6*sin(2*pi*f6*t + phi6*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase6value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 7th harmonic
function amp7slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp7value,'String',num2str(val));
  a7 = get(amp7slider,'Value');
  y7 = a7*sin(2*pi*f7*t + phi7*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp7value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 7th harmonic
function phase7slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase7value,'String',num2str(val));
  phi7 = get(phase7slider,'Value');
  y7 = a7*sin(2*pi*f7*t + phi7*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase7value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

%%%%%%%%%%%%%%%%
% amplitude contribution from 8th harmonic
function amp8slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(amp8value,'String',num2str(val));
  a8 = get(amp8slider,'Value');
  y8 = a8*sin(2*pi*f8*t + phi8*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function amp8value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

% phase contribution from 8th harmonic
function phase8slider_callback(hObject,eventdata)
  delete(waveformplot)

  val = get(hObject,'Value');
  set(phase8value,'String',num2str(val));
  phi8 = get(phase8slider,'Value');
  y8 = a8*sin(2*pi*f8*t + phi8*pi/180);
  y = y1+y2+y3+y4+y5+y6+y7+y8;

  %define plot (with labels and limits)
  waveformplot=axes('Units','Pixels','Position',[80 450 800 160]);
  set(f,'CurrentAxes',waveformplot)
  %plot(t,y);
  plot(t,y1,'b',t,y2,'b',t,y3,'b',t,y4,'b',t,y5,'b',t,y6,'b',t,y7,'b',t,y8,'b',t,y,'r');
  grid on
  xlim([xlimmin xlimmax])
  ylim([ylimmin ylimmax])
  xlabel('time (seconds)')
  ylabel('summed waveform')

end 
     
function phase8value_callback(hObject,eventdata)
  val = str2double(get(hObject,'String'));
end 

end
