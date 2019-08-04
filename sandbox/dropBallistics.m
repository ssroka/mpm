









clear;clc;

theta = 90;
U10 = 20;
t = linspace(0,10);
y0 = 0;
x0 = 0;
a  = -9.81;
v0  = [U10*cosd(theta);U10*sind(theta)];



y = 0.5*a*t.^2 + v0(2)*t + y0;
x = v0(1)*t + x0;

t(y<0)=[];
x(y<0)=[];
y(y<0)=[];


plot(t,y,'r')
hold on
plot(t,x,'b')



















