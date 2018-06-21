function []= plotRPP(theta, r,alpha)

a1=pi; %set baseline
a2=a1-theta; %set conic sector 
a3=a1+theta;

%plot bottom half of conic sector
t=linspace(a1,a2);
x=r*cos(t);
y=r*sin(t);

%plot top half of conic sector
t1=linspace(a1,a3);
x1=r*cos(t1);
y1=r*sin(t1);


hold on
plot([0,x,0],[0,y,0])
plot([0,x1,0],[0,y1,0])
line([-alpha -alpha], [-r r])
hold off

end
