% find discrete frequencies for free bar (no tension)
close all

N = 10001;
x = linspace(0,12,N); % x = kL
y1 = tanh(x/2);
y2 = -tanh(x/2);
y3 = tan(x/2);

% find first intersection
ndx = find(y3>y2 & y3<0);
x1 = x(ndx(1));
fprintf('first intersection: kL*2/pi = %f\n', x1*2/pi);

% find second intersection
ndx = find(y3<y1 & y3>0);
x2 = x(ndx(end));
fprintf('second intersection: kL*2/pi = %f\n', x2*2/pi);

% make plot
plot(x,y1,'r',x,y2,'r',x,y3,'b');
ylim([-2 2])
xlabel('kL')
grid on

return

