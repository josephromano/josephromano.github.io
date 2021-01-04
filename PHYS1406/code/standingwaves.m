function standingwaves(n)
%
% demonstrate standing waves on a string fixed at
% both ends (variable freq)
%
% n = 1, 2, 3 ... is number of harmonic
% 
% n different from an integer doesn't produce a standing wave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%define figure
f=figure('Position',[300 400 900 640],...
         'Name','Standing waves');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot parameters
L = 1; % length of string
totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);

xlimmin = 0;
xlimmax = L;
ylimmin = -12;
ylimmax = 12;

%%%%%%%%%%%%%%%%%%%%%%%
% parameters for string
v = 1; % velocity of sound in string

% frequency
f0 = n*v/(2*L);

% discrete times
t_tot = 15*L/v;
delta_t = (L/v)/100; % 100 delta_T to travel down the arm
num_t = 1+floor(t_tot/delta_t);
t = transpose(linspace(0, t_tot, num_t));

% discrete x's
delta_x = v*delta_t; 
num_x = 1+floor(L/delta_x);
x = transpose(linspace(0, L, num_x));

% damping factor 
d = 2/3;
d = 1;

% waveforms
a = 1;
N = 6;
yl = zeros(num_x,N);
yr = zeros(num_x,N);
y = zeros(num_x,1);

for k=1:num_t

    delete(totwaveformplot)
    delete(indwaveformplot)
   
    % plot right- and left-moving waves
    indwaveformplot=axes('Units','Pixels','Position',[80 250 750 160]);
    set(f,'CurrentAxes',indwaveformplot)
    plot(x,yr(:,1),x,yl(:,1),x,yr(:,2),x,yl(:,2),x,yr(:,3),x,yl(:,3),x,yr(:,4),x,yl(:,4),x,yr(:,5),x,yl(:,5),x,yr(:,6),x,yl(:,6));
    grid on
    xlim([xlimmin xlimmax])
    ylim([ylimmin ylimmax])
    %ylim([-1 1])
    xlabel('x (m)')

    % plot summed waveform
    totwaveformplot=axes('Units','Pixels','Position',[80 450 750 160]);
    set(f,'CurrentAxes',totwaveformplot)
    temp = zeros(num_x,1);
    for jj=1:N
      temp = temp + yl(:,jj) + yr(:,jj);
    end
    y = temp;
    plot(x,y,'-r');
    grid on
    xlim([xlimmin xlimmax])
    ylim([ylimmin ylimmax])
    titstr = ['elapsed time = ' num2str(t(k)) ' sec'];
    title(titstr)
 
    % create new wave profiles (x-axis)
    yr(:,1) = [a*sin(2*pi*f0*t(k)); yr(1:end-1,1)];
    yl(:,1) = [yl(2:end,1); -yr(end,1)];
    for jj=2:N
      yr(:,jj) = [-d*yl(1,jj-1); yr(1:end-1,jj)];
      yl(:,jj) = [yl(2:end,jj); -yr(end,jj)];
    end

    pause(.00001)

end 
    
return
 
