clear all
close all

%parameters
g=1;
niu=1;
delta=1;
niu_oil=1;

%domain
x=linspace(-40,40,1000);
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

% Free stream 
u_free=0*repmat(x+60,length(y),1);

%Initial conditions
u=u_free;

q=u(1,:);
Sigma1=zeros(length(y),length(x));
Sigma2=zeros(length(y),length(x));

%Lax-Friedrichs
e=ones(length(x),1);
d3=spdiags([e,0*e,e],-1:1,length(x),length(x));
d3(2,1)=2;

d1=spdiags([e,0*e,-1*e],-1:1,length(x),length(x));
d1(2,1)=2;

d=spdiags([-1*e,1*e],0:1,length(x),length(x));

%Upwind 
d2=spdiags([1*e,-1*e],0:1,length(x),length(x));

  

% Second derivative matrix
e=ones(length(y),1);
D=spdiags([e,-2*e,e],-1:1,length(y),length(y));
D(1,1)=-2; D(1,2)=2;
D(length(y),length(y))=-1;
I=eye(length(y));
 

%Export results onto videos
%vidobj1=VideoWriter('Flowingwater_velocity');
%open(vidobj1);

 
 t=0;


for m=1:10
    
dt=min(0.4*dx/(max(abs(q))+1),dx^2/4);


% Shear Stress at the interface
sigma1=g*delta/(2*niu)*d*((h').^2)/dx-niu_oil/niu*d*(h'.*(d*q')/dx)/dx; 
Sigma1(1,:)=sigma1;



DDy=D*u/dy^2;


% thickness at next time step

h=h*d3/2-dt/(2*dx)*(q.*h)*d1;
h1=max(h-0.00003*max(h),0);
h=h1*V/(sum(h1)*dx);











% update the shear stress
sigma2=g*delta/(2*niu)*d*((h').^2)/dx-niu_oil/niu*d*(h'.*(d*q')/dx)/dx;
Sigma2(1,:)=sigma2;



t=t+dt;


  


 
 A=I-0.5*niu*dt*D/dy^2;
 B=u+niu*dt*(0.5*DDy-0.5*2*Sigma1/dy-0.5*2*Sigma2/dy);
 
 u=A\B;
 u=u-(dt/dx)*u.*(u*d2)+(dt/dx)*u_free.*(u_free*d2);
 

 u(length(y),:)=u_free(length(y),:);
 %u(:,1)=u_free(:,1);
 q=u(1,:);
 
 

 [mq,lo]=max(q);
 Lo(m)=(lo-1)*dx;
 
 [hmax,xmax]=max(h);
 center(m)=(xmax-0.5*length(x))*dx;
 edge(m)=(find(h,1,'last')-0.5*length(x))*dx;

 

 
 %Boundary layer
%  for i=1:length(x)
%      I1=find(u(:,i)<=0.8*max(u(:,i)));
%      I2=find(u(:,i)<=0.5*max(u(:,i)));
%      I3=find(u(:,i)<=0.1*max(u(:,i)));
%      boundary1(i)=-dy*(I1(1)-1);
%      boundary2(i)=-dy*(I2(1)-1);
%      boundary3(i)=-dy*(I3(1)-1);
%  end
 
%  figure(1)
%  plot(x,h,'Linewidth',2)
%  axis([x(1),x(length(x)),0,2])
%  xlabel('x','FontSize',20)
%  ylabel('h','FontSize',20)
%  M1(m)=getframe(gcf);
%  writeVideo(vidobj1,M1(m));
 
 T(m)=t;
 
 figure(2)
 p=plot(x,q,x,u(length(y),:))
 p(1).LineWidth=2;
 p(2).LineWidth=2;
 xlabel('x','FontSize',20)
 ylabel('y','FontSize',20)
 leg=legend('Oil velocity','Free stream velocity')
 set(leg,'FontSize',8)
 M2(m)=getframe(gcf);
 %writeVideo(vidobj1,M2(m));
  
%  figure(2)
%  plot(x/(t+1)^(3/8), h*(t+1)^(3/8),'Linewidth',2)
%  axis([0,20,0,2])
%  xlabel('scaled x','FontSize',20)
%  ylabel('scaled h','FontSize',20)
%  M2(m)=getframe(gcf);
%  writeVideo(vidobj1,M2(m))
%  
%  
%  figure(3)
%  plot(x,boundary1,'r',x,boundary2,'b',x,boundary3,'g')
%  legend('80% boundary','50% boundary','10% boundary')
%  axis([0,20,-0.5,0])
%  xlabel('x')
%  ylabel('y')
%  M3(m)=getframe;
%  writeVideo(vidobj1,M3(m))
 
 
end

%close(vidobj1)



 
 
     
 
 
 
 
 


