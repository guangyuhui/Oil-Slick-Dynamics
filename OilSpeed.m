clear all
close all

%parameters
g=1;
niu=1;
delta=1;

%domain
x=linspace(0,20,500);
dx=x(2)-x(1);
y=linspace(0,0.5,100);
dy=y(2)-y(1);

%initial state
for i=1:length(x)
    if x(i)>=0 && x(i)<=2/sqrt(pi)
        %h(i)=3+3*sin(x(i));
        %h(i)=sqrt(10-0.1*x(i)^2);
        h(i)=sqrt(4/pi-x(i)^2);
    else 
        h(i)=0;
    end
end
%h(1)=1;

u=zeros(length(y),length(x));
q=u(1,:);




