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
Sigma1=zeros(length(y)-1,length(x));
Sigma2=zeros(length(y)-1,length(x));

%Lax-Friedrichs
e=ones(length(x),1);
d3=spdiags([e,0*e,e],-1:1,length(x),length(x));
d3(2,1)=2;

d1=spdiags([e,0*e,-1*e],-1:1,length(x),length(x));
d1(2,1)=2;

d=spdiags([-1*e,1*e],0:1,length(x),length(x));

  

% Second derivative matrix
e=ones(length(y)-1,1);
D=spdiags([e,-2*e,e],-1:1,length(y)-1,length(y)-1);
D(1,1)=-2; D(1,2)=2;
I=eye(length(y)-1);
 

 

 
 t=0;


for m=1:500
    
dt=0.2*dx/(max(abs(q))+1);


% Shear Stress at the interface
sigma1=g*delta/(2*niu)*d*((h').^2)/dx; 
Sigma1(1,:)=sigma1;



 DDy=D*u(1:length(y)-1,:)/dy^2;


% thickness at next time step

h=h*d3/2-dt/(2*dx)*(q.*h)*d1;




% Volume of the slick
V=sum(h*dx);





% update the shear stress
sigma2=g*delta/(2*niu)*d*((h').^2)/dx;
Sigma2(1,:)=sigma2;



t=t+dt;


  


 
 A=I-0.5*niu*dt*D/dy^2;
 B=u(1:length(y)-1,:)+niu*dt*(0.5*DDy-0.5*2*Sigma1/dy-0.5*2*Sigma2/dy);
 
 u(1:length(y)-1,:)=A\B;
 

 u(length(y),:)=0;
 q=u(1,:);
 

 [mq,lo]=max(q);
 ratio(m)=(lo-1)*dx/mq;
 

 

 
 %Boundary layer
 for i=1:length(x)
     I1=find(u(:,i)<=0.8*max(u(:,i)));
     I2=find(u(:,i)<=0.5*max(u(:,i)));
     I3=find(u(:,i)<=0.1*max(u(:,i)));
     boundary1(i)=-dy*(I1(1)-1);
     boundary2(i)=-dy*(I2(1)-1);
     boundary3(i)=-dy*(I3(1)-1);
 end
 
 figure(1)
 plot(x,h)
 axis([0,20,0,sqrt(4/pi)])
 M1(m)=getframe;
 
 T(m)=t;
  
 figure(2)
 plot(x/(t+1)^(3/8), h*(t+1)^(3/8))
 axis([0,20,0,sqrt(4/pi)])
 M2(m)=getframe;
 
 
 figure(3)
 plot(x,boundary1,'r',x,boundary2,'b',x,boundary3,'g')
 axis([0,20,-0.5,0])
 M3(m)=getframe;
 
 
end



 
 
     
 
 
 
 
 


