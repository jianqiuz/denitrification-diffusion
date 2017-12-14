%%%diffusion on a static 1D domain

P(1)=0.0244;%diffusion coefficient in cm^2/s
P(2)=1.23;%initial concentration at 0 cm, in ng N2O/cm^3
P(3)=1.33*10^(-4); %gas production rate p in ng/cm^3/s (from simple flux estimation model based on field measurement)
P(4)=760*10^(-4);% Max. rate of gas uptake Umax in ng/cm^3/s (from simple flux estimation model based on field measurement)
P(5)=35;%Michaelis constant of uptake rate Km in ng/cm^3/ from Gu's, Conrad's, Laudone's Hos
P(6)=0.61; % Henry's law constant at 25 degree C ( Caq/Cgas)
L=25; %length of the domain in cm
maxt=43200; % Max. simulation time in s 48 hours
m=0; %shape of the PDE
t=linspace(0,maxt,100)% tspan5
x=linspace(0,L,100)%xmesh


%call of PDEPE
sol=pdepe(m,@DiffusionPDEfun, @DiffusionICfun, @DiffusionBCfun, x,t,[],P);
u=sol;

%%%surface plot
figure(1)
surf(x,t,FN,'edgecolor','none');
xlabel('Depth (cm)','fontsize',12,'fontweight','b','fontname','times new roman')
ylabel('Time (s)','fontsize',12,'fontweight','b','fontname','times new roman')
zlabel('N2O(ngN2O/cm^3)','fontsize',12,'fontweight','b','fontname','times new roman')
axis([0 L 0 maxt 0 40])
set(gcf(), 'Renderer','painters')
set(gca, 'fontsize',10,'fontweight','b','fontname','times new roman')

surf(x,t,SG12,'edgecolor','none');

subplot(1,2,1)
p(1)=surf(x,t,FWTp,'edgecolor','none');
xlabel('Depth (cm)','fontsize',16,'fontweight','b','fontname','times new roman')
ylabel('Time (h)','fontsize',16,'fontweight','b','fontname','times new roman')

axis([0 L 0 maxt 0 40])
subplot(1,2,2)
p(2)=surf(x,t,FWp,'edgecolor','none');
xlabel('Depth (cm)','fontsize',16,'fontweight','b','fontname','times new roman')
ylabel('Time (h)','fontsize',16,'fontweight','b','fontname','times new roman')

axis([0 L 0 maxt 0 40])


linkaxes(p,t);




%%2-D line plot
figure(2)
hold all
 n=linspace(1,length(t),10);
     n(2)=4;
     n(3)=8;
     n(4)=13;
     n(5)=25;
     
     
     plot(sol(n,:),x, 'LineWidth',1)
 
 

 
 xlabel('N2O concentraction(ngN2O/cm^3)','fontsize',20,'fontweight','b','fontname','arial')
 ylabel('Depth (cm)','fontsize',20,'fontweight','b','fontname','arial')
 axis([0 50 0 L])
 set(gca, 'fontsize',15,'fontweight','b','fontname','arial')
 
 
%%2-D line plot
figure(3)
hold all
 for n=linspace(1,length(t),10)
     plot(t,sol(:,n),'LineWidth',2)
 end
 xlabel('Time t','fontsize',20,'fontweight','b','fontname','arial')
 ylabel('Species u','fontsize',20,'fontweight','b','fontname','arial')
 axis([0 maxt 0 50])
 set(gca, 'fontsize',15,'fontweight','b','fontname','arial')

 %%%net production
 fid=fopen('PDE solution1.csv');
 out=textscan(fid,'%s%f$f','delimiter',',');
 fclose(fid);
 M=csvread('PDE solution1.csv');
 
 surf(x,t,M,'edgecolor','none');
xlabel('Depth (cm)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Time (s)','fontsize',20,'fontweight','b','fontname','arial')
zlabel('N2O concentraction(ngN2O/cm^3)','fontsize',20,'fontweight','b','fontname','arial')
axis([0 L 0 maxt 0 100])
set(gcf(), 'Renderer','painters')
set(gca, 'fontsize',15,'fontweight','b','fontname','arial')

for i=1:99
    if M(i+1,i)<M(i,i)
        M(i+1,i)=M(i,i)+0.001;
  
    end
end


for i=1:99
    if M(i,i+1)<M(i,i)
        M(i,i+1)=M(i,i)+0.001;
  
    end
end

 
 
 Y=zeros(100);
 
 
 for i=2:99
     for j=1:99
         Y(j,i+1)=Y(j,i)+x(2)*((M(j+1,i)-M(j,i))/t(2)-P(1)*(M(j+1,i+1)-2*M(j+1,i)+M(j+1,i-1)+M(j,i+1)-2*M(j,i)+M(j,i-1))/(2*x(2)*x(2)));
     end
 end
 
 
 
 for i=1:100
     for j=1:100
    if Y(i,j)<0
        Y(i,j)=0;
    end
        
  
    end
end


for i=1:99
    if M(i,i+1)<M(i,i)
        M(i,i+1)=M(i,i)+0.001;
  
    end
end
 
 figure(1)
surf(x,t,Y,'edgecolor','none');
xlabel('Depth (cm)','fontsize',20,'fontweight','b','fontname','arial')
ylabel('Time (s)','fontsize',20,'fontweight','b','fontname','arial')
zlabel('N2O production(ngN2O/g soil/s)','fontsize',20,'fontweight','b','fontname','arial')
axis([0 L 0 maxt 0 0.06])
set(gcf(), 'Renderer','painters')
set(gca, 'fontsize',15,'fontweight','b','fontname','arial')
         
     
