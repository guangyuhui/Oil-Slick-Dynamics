clear all
close all

%parameters
g=1;
niu=1;
delta=1;

%domain
x=linspace(-20,20,1000);
dx=x(2)-x(1);
y=linspace(0,0.5,100);
dy=y(2)-y(1);

%initial state
for i=1:length(x)
    if x(i)>=-2/sqrt(pi) && x(i)<=2/sqrt(pi)
        %h(i)=3+3*sin(x(i));
        %h(i)=sqrt(10-0.1*x(i)^2);
        h(i)=sqrt(4/pi-x(i)^2);
    else 
        h(i)=0;
    end
end
%h(1)=1;

% Volume of the slick
V=sum(h*dx);

u=10*repmat(x+60,length(y),1);

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
 

%Export results onto videos
%vidobj1=VideoWriter('Boundarylayer');
%open(vidobj1);

 
 t=0;


for m=1:3000
    
dt=0.4*dx/(max(abs(q))+1);


% Shear Stress at the interface
sigma1=g*delta/(2*niu)*d*((h').^2)/dx; 
Sigma1(1,:)=sigma1;



DDy=D*u(1:length(y)-1,:)/dy^2;


% thickness at next time step

h=h*d3/2-dt/(2*dx)*(q.*h)*d1;
%h1=max(h-0.005*max(h),0);
%h=h1*V/(sum(h1)*dx);











% update the shear stress
sigma2=g*delta/(2*niu)*d*((h').^2)/dx;
Sigma2(1,:)=sigma2;



t=t+dt;


  


 
 A=I-0.5*niu*dt*D/dy^2;
 B=u(1:length(y)-1,:)+niu*dt*(0.5*DDy-0.5*2*Sigma1/dy-0.5*2*Sigma2/dy);
 
 u(1:length(y)-1,:)=A\B;
 

 u(length(y),:)=10*(x+60);
 q=u(1,:);
 
 

 [mq,lo]=max(q);
 Lo(m)=(lo-1)*dx;
 
 %edge(m)=(find(h==0,1)-1)*dx;

 

 
 %Boundary layer
%  for i=1:length(x)
%      I1=find(u(:,i)<=0.8*max(u(:,i)));
%      I2=find(u(:,i)<=0.5*max(u(:,i)));
%      I3=find(u(:,i)<=0.1*max(u(:,i)));
%      boundary1(i)=-dy*(I1(1)-1);
%      boundary2(i)=-dy*(I2(1)-1);
%      boundary3(i)=-dy*(I3(1)-1);
%  end
 
 figure(1)
 plot(x,h)
 axis([-20,20,0,2])
 xlabel('x')
 ylabel('h')
 M1(m)=getframe;
%  writeVideo(vidobj1,M1(m));
 
 T(m)=t;
  
 figure(2)
 plot(x/(t+1)^(3/8), h*(t+1)^(3/8))
 axis([-20,20,0,sqrt(4/pi)])
 xlabel('scaled x')
 ylabel('scaled h')
 M2(m)=getframe;
%  writeVideo(vidobj1,M2(m))
%  
%  
%  figure(3)
%  plot(x,boundary1,'r',x,boundary2,'b',x,boundary3,'g')
%  legend('80% boundary','50% boundary','10% boundary')
%  axis([-20,20,-0.5,0])
%  xlabel('x')
%  ylabel('y')
%  M3(m)=getframe;
 %writeVideo(vidobj1,M3(m))
 
 
end

%close(vidobj1)



 
 
     
 
 
 
 
 


