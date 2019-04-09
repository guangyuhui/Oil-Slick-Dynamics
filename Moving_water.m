clear all
close all

%parameters
g=1;
niu=1;
delta=1;
niu_oil=1;
eta=1;
alpha = 0; % evaporation rate

%domain
x=linspace(0,80,2000);
dx=x(2)-x(1);
y=linspace(0,0.5,100);
dy=y(2)-y(1);

%initial state
for i=1:length(x)
    if x(i)<=4/sqrt(pi)
        %h(i)=3+3*sin(x(i));
        %h(i)=sqrt(10-0.1*x(i)^2);
        %h(i)=sqrt(4/pi-x(i)^2);
        h(i) = sqrt(1/pi-x(i)^2/16);
    else 
        h(i)=0;
    end
end
%h(1)=1;

% Volume of the slick
V=sum(h*dx);

% Free stream 
u_free=eta*repmat(sin(x),length(y),1);

%Initial conditions
u=u_free;
 

q=0;
Sigma1=zeros(length(y),length(x));
Sigma2=zeros(length(y),length(x));

%Lax-Friedrichs
e=ones(length(x),1);
d3=spdiags([e,0*e,e],-1:1,length(x),length(x));
d3(2,1)=2;

d1=spdiags([e,0*e,-1*e],-1:1,length(x),length(x));
d1(2,1)=2;

d=spdiags([-1*e,1*e],0:1,length(x),length(x));
d(1,1)=0;

%Upwind 
d2=spdiags([1*e,-1*e],0:1,length(x),length(x));
d(2,1)=1;

  

% Second derivative matrix
e=ones(length(y),1);
D=spdiags([e,-2*e,e],-1:1,length(y),length(y));
D(1,1)=-1; D(1,2)=1;
%D(length(y),length(y))=-1;
I=eye(length(y));
 

%Export results onto videos
%vidobj1=VideoWriter('Flowingwater_velocity');
%open(vidobj1);

 
 t=0;


for m=1:235
    
    dt=0.05;%0.5*dx/(max(abs(q))+1);%,dx^2/4);
    
    u_old = u;


    if m==1
    
    % Shear Stress at the interface
    %sigma1=g*delta/(2*niu)*d*((h').^2)/dx-niu_oil/niu*d*(h'.*(d*q')/dx)/dx;
         sigma1=h/niu.*(1/2*(u(1,:).^2*d')/dx);
         Sigma1(1,:)=sigma1;
    end
% 
% 
% 
%     DDy=D*u/dy^2;
% 
% 
    % thickness at next time step

    h=h*d3/2-dt/(2*dx)*(q.*h)*d1-dt*alpha*h;
    
    h1=max(h-0.0000003*max(h),0);
    h=h1*V/(sum(h1)*dx);











%     % update the shear stress
%     %sigma2=g*delta/(2*niu)*d*((h').^2)/dx-niu_oil/niu*d*(h'.*(d*q')/dx)/dx;
%     sigma2=g*delta/niu*h'.*(d*h')/dx;
%     Sigma2(1,:)=sigma2;



    t=t+dt;


  
     d_q_sq=1/2*(u(1,:).^2*d')/dx;
     
     for j=1:length(x)
         I_uj=I;
         I_uj(1,1)=1+h(j)/(2*dy);
         A=I_uj-(niu*dt/(2*dy^2))*D;
         %A=I_uj-0.25*niu*dt*D/dy^2;
         B=I_uj*u(:,j)+0.5*niu*dt*D/dy^2*u(:,j)-0.5*niu*dt/dy*Sigma1(:,j)-dt/(2*dy)*h(j)*[d_q_sq(j);zeros(length(y)-1,1)];
         u(:,j)=A\B;
     end

     u=u-(0.5*dt/dx)*((u.^2)*d2);%+(dt/dx)*u_free.*(u_free*d2);
     
     
     % shear stress at the interface
     sigma1=h/niu.*(1/2*(u(1,:).^2*d')/dx+(u(1,:)-u_old(1,:))/dt);
     Sigma1(1,:)=sigma1;


     
     u(length(y),:)=u_free(length(y),:);
     u(:,1)=u_free(:,1);
     q=u(1,:);
     
     acceleration(m,:)=(u_old(1,:)-q)/dt;
     advection(m,:) = 1/2*(q.^2)*d1;
     
 



     [mq,lo]=max(q);
     Lo(m)=(lo-1)*dx;

     [hmax,xmax]=max(h);

     % half slick
     %center(m)=(xmax-1)*dx;
     %edge(m)=(find(h,1,'last')-1)*dx;

     % whole slick
    %  center(m)=(xmax-0.5*length(x))*dx;
    %  edge(m)=(find(h,1,'last')-0.5*length(x))*dx;
    
     expectation(m) = sum((x.^1).*h*dx);
     edge(m) = floor(expectation(m)/dx);

     vis_w = niu*sigma1;
     gravity = g*delta*0.5*(h.^2)*d1;
     ratio_g_vw(m)= gravity(1)./vis_w(1);
     %max_ratio_g_vw(m) = max(ratio_g_vw(1:edge(m)));


%      %Boundary layer
%    if m>5
%      for i=2:length(x)
%          I1=find(u(:,i)>0.95*max(u(:,i)));
%          %I2=find(u(:,i)<=0.5*max(u(:,i)));
%          %I3=find(u(:,i)<=0.1*max(u(:,i)));
%          boundary1(i)=-dy*(I1(1)-1);
%          %boundary2(i)=-dy*(I2(1)-1);
%          %boundary3(i)=-dy*(I3(1)-1);
%      end
%    end

     boundary(m,:)=(u_free(1,:)-q)./sigma1;
     
     
     
%      figure(1)
%      plot(x,h,'Linewidth',2)
%      axis([x(1),x(length(x)),0,2])
%      xlabel('x','FontSize',20)
%      ylabel('h','FontSize',20)
%      M1(m)=getframe(gcf);
    %  writeVideo(vidobj1,M1(m));

     T(m)=t;

    %  figure(2)
    %  p=plot(x,q,x,u(length(y),:))
    %  p(1).LineWidth=2;
    %  p(2).LineWidth=2;
    %  xlabel('x','FontSize',20)
    %  ylabel('y','FontSize',20)
    %  leg=legend('Oil velocity','Free stream velocity')
    %  set(leg,'FontSize',8)
    %  M2(m)=getframe(gcf);
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
    
%     if m>5
%      figure(3)
%      plot(x,boundary1,'r','LineWidth',2)
%      legend('95% boundary')
%      axis([0,20,-0.5,0])
%      xlabel('x')
%      ylabel('y')
%      M3(m)=getframe;
%     %  writeVideo(vidobj1,M3(m))
%     end

     if edge(m)==80
         break;
     end
         
 
 
end



% [order_t,order_y]=ode45(@(t,y) moving_order_eq(t,y,niu,V,eta),[0,t],[2/sqrt(pi);0]); 
% 
% 
% 
% figure(2)
% 
% plot(T,log(expectation),'Linewidth',2,'color','r')
% xlabel('T','FontSize',20)
% ylabel('log(R(t))','FontSize',20)
% hold on 
% plot(order_t,log(order_y(:,1)),'Linewidth',2,'color','b')


%close(vidobj1)
