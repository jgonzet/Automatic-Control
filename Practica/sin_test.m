close all;
t=-20*pi:pi/200:20*pi;
a=3;
a=2;
b=-2;
w=.5;



o=a*sin(w*t)+b*cos(w*t);
u=2.83*sin(w*t-0.785);

plot(t,o)
hold
plot(t,u)