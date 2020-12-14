clc; clear;

%Constants set (for all Methods)

Cm=0.01; % Membrane Capcitance uF/cm^2
dt=0.04; % Time Step ms
t=0:dt:25; %Time Array ms
I=0.1; %External Current Applied
ENa=55.17; % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003 % mS/cm^2 Leakage conductance
V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value

for i=1:length(t)-1
%Euler method to find the next m/n/h value
m(i+1)=m(i)+dt*((am(V(i))*(1-m(i)))-(bm(V(i))*m(i)));
n(i+1)=n(i)+dt*((an(V(i))*(1-n(i)))-(bn(V(i))*n(i)));
h(i+1)=h(i)+dt*((ah(V(i))*(1-h(i)))-(bh(V(i))*h(i)));
gNa=gbarNa*m(i)^3*h(i);
gK=gbarK*n(i)^4;
gl=gbarl;
INa=gNa*(V(i)-ENa);
IK=gK*(V(i)-EK);
Il=gl*(V(i)-El);
%Euler method to find the next voltage value
V(i+1)=V(i)+(dt)*((1/Cm)*(I-(INa+IK+Il)));
end

%Store variables for graphing later

FE=V;
FEm=m;
FEn=n;
FEh=h;
clear V m n h;

%-----------------------------

% Runge-Kutta Method

V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value

for i=1:length(t)-1 % Loop through each step until time is finished
%4 step method of Runge-Kutta
K1=dt*HH(i,[V(i); n(i); m(i); h(i)]);
k1=K1(1,1);n1=K1(2,1);m1=K1(3,1);h1=K1(4,1);% obtain 4 k variables (V,m,n,h) from HH function
K2=dt*HH(i+(0.5*dt),[V(i)+(0.5*k1);n(i)+(0.5*n1);m(i)+(0.5*m1);h(i)+(0.5*h1)]);
k2=K2(1,1);n2=K2(2,1);m2=K2(3,1);h2=K2(4,1);
K3=dt*HH(i+(0.5*dt),[V(i)+(0.5*k2);n(i)+(0.5*n2);m(i)+(0.5*m2);h(i)+(0.5*h2)]);
k3=K3(1,1);n3=K3(2,1);m3=K3(3,1);h3=K3(4,1);
K4=dt*HH(i+dt,[V(i)+k3;n(i)+n3;m(i)+m3;h(i)+h3]);
k4=K4(1,1);n4=K4(2,1);m4=K4(3,1);h4=K4(4,1);
%create next step for each variable
V(i+1)=V(i)+1/6*(k1+2*k2+2*k3+k4);
n(i+1)=n(i)+1/6*(n1+2*n2+2*n3+n4);
m(i+1)=m(i)+1/6*(m1+2*m2+2*m3+m4);
h(i+1)=h(i)+1/6*(h1+2*h2+2*h3+h4);
end

%set variables for graphing later

RK=V;
RKm=m;
RKn=n;
RKh=h;
clear V m n h;

%-------------------

% Predictor-Corrector Method

V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value

%First four steps are found using the Runge-Kutta Method

for i=1:3
K1=dt*HH(i,[V(i); n(i); m(i); h(i)]);
k1=K1(1,1);n1=K1(2,1);m1=K1(3,1);h1=K1(4,1);
K2=dt*HH(i+(0.5*dt),[V(i)+(0.5*k1);n(i)+(0.5*n1);m(i)+(0.5*m1);h(i)+(0.5*h1)]);
k2=K2(1,1);n2=K2(2,1);m2=K2(3,1);h2=K2(4,1);
K3=dt*HH(i+(0.5*dt),[V(i)+(0.5*k2);n(i)+(0.5*n2);m(i)+(0.5*m2);h(i)+(0.5*h2)]);
k3=K3(1,1);n3=K3(2,1);m3=K3(3,1);h3=K3(4,1);
K4=dt*HH(i+dt,[V(i)+k3;n(i)+n3;m(i)+m3;h(i)+h3]);
k4=K4(1,1);n4=K4(2,1);m4=K4(3,1);h4=K4(4,1);
V(i+1)=V(i)+1/6*(k1+2*k2+2*k3+k4);
n(i+1)=n(i)+1/6*(n1+2*n2+2*n3+n4);
m(i+1)=m(i)+1/6*(m1+2*m2+2*m3+m4);
h(i+1)=h(i)+1/6*(h1+2*h2+2*h3+h4);
end

for i = 4:length(t)-1 % P-C Method
%predictor
yp=[V(i);n(i);m(i);h(i)]+(dt/24)*(55*HH(t(i),[V(i);n(i);m(i);h(i)])-59*(HH(t(i-1),[V(i-1);n(i-1);m(i-1);h(i-1)]))+37*(HH(t(i-2),[V(i-2);n(i-2);m(i-2);h(i-2)]))-9*(HH(t(i-3),[V(i-3);n(i-3);m(i-3);h(i-3)])));
%corrector
yc=[V(i);n(i);m(i);h(i)]+(dt/24)*(9*(HH(t(i+1),yp))+19*(HH(t(i),[V(i);n(i);m(i);h(i)]))-5*HH(t(i-1),[V(i-1);n(i-1);m(i-1);h(i-1)])+HH(t(i-2),[V(i-2);n(i-2);m(i-2);h(i-2)]));
C=yc+(19/270)*(yp-yc);
V(i+1)=C(1,1);
n(i+1)=C(2,1);
m(i+1)=C(3,1);
h(i+1)=C(4,1);
end

%Store variables for graphing
PC=V;
PCm=m;
PCn=n;
PCh=h;
clear V m n h;

%------------------------

% ODE45 Method

V=-60; % Initial Membrane voltage
m=am(V)/(am(V)+bm(V)); % Initial m-value
n=an(V)/(an(V)+bn(V)); % Initial n-value
h=ah(V)/(ah(V)+bh(V)); % Initial h-value
y0=[V;n;m;h];
tspan = [0,max(t)];
%Matlab's ode45 function
[time,V] = ode45(@HH,tspan,y0);
OD=V(:,1);
ODn=V(:,2);
ODm=V(:,3);
ODh=V(:,4);

clear V;

%-------

%Plotting the functions

plot(t,FE,t,RK,t,PC,time,OD);

legend('Forward Euler','Runge-Kutta','Predictor Corrector','ODE45');

xlabel('Time (ms)');

ylabel('Voltage (mV)');

title('Voltage Change for Hodgkin-Huxley Model');

figure

plot(t,FEn,'b',t,RKn,'b:',t,PCn,'b-.',time,ODn,'b--',t,FEm,'g',t,RKm,'g:',t,PCm,'g-.',time,ODm,'g--',t,FEh,'r',t,RKh,'r:',t,PCh,'r-.',time,ODh,'r--');

ylabel('Gaining Variables')

xlabel('Time (ms)')

axis([0 5 0 1])

legend('n Euler','n Runge-Kutta','n Predictor Corrector','n 0ODE45','m Euler','m Runge-Kutta','m Predictor Corrector','m ODE45','h Euler','h Runge-Kutta','h Predictor Corrector','h ODE45');