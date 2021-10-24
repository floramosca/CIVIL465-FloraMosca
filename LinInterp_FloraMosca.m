%Title: Response of a four floors frame to the Tolmezzo earthquake using
%linear interpolation 
%Content: MDOF analysis and linear interpolation of results
%This code is developed based on the inputs of the Tolmezzo earthquake 
%ground acceleration from Pr. Pierino Lestuzzi course CIVIL420 
%This code uses the simplification that the response of a MDOF can be
%appoximated by taking into account only the first mode of vibration of the
%system
%Author: Flora Mosca
%Created: 9/01/2021
%Last update: 24/10/2021
clc
clear
close all

%% frame

%data
h=4;

M1=20e3;
M2=16e3;
M3=24e3;
M4=20e3;

I1=136.7e-6;
I2=77.6e-6;

E=210e9;

dumping=5/100;

load Tolm.txt
Tolm=Tolm';

%mass matrix
M=[M1 0 0 0;0 M2 0 0;0 0 M3 0;0 0 0 M4];

%height matrix
H=[h;2*h;3*h;4*h];

%stiffness matrix
k1=3*E*I1/h^3;
k2=12*E*I2/h^3;

k=k1+k2; %floor stiffness
K=[2 -1 0 0;-1 2 -1 0; 0 -1 2 -1;0 0 -1 1]*k;

%natural frequences
[V,W2]=eig(K,M);
w=sqrt(diag(W2));
phimax=(max(abs(V))'*diag(eye(4))')';
phi=V./phimax;
A=phi;
f=w/2/pi;

%****************************************************************
%plotting the modes' shapes
%****************************************************************

figure()
height=[0;H];
shape=[0;0;0;0;0];
firstModeShape=[0;-A(:,1)];
secondModeShape=[0;-A(:,2)];
thirdModeShape=[0;-A(:,3)];
fourthModeShape=[0;-A(:,4)];
p=plot(shape,height,'k','DisplayName','Undeformed shape');
p.LineWidth = 1.5;
hold on
plot(firstModeShape,height,'DisplayName','First mode')
plot(secondModeShape,height,'DisplayName','Second mode')
plot(thirdModeShape,height,'DisplayName','Third mode')
plot(fourthModeShape,height,'DisplayName','Fourth mode')

ylabel('h [m]')
xlabel('d [-]')
legend('Location','southeast')
title('Normalised modes shape')

%****************************************************************
%performing the modal analysis
%****************************************************************

%modal matrices
Mstar=A'*M*A;
Kstar=A'*K*A;

%participation factors vector
r=A*diag(eye(4));

%participation factors
rsm=r./diag(Mstar);

%effective modal mass
Mmod=rsm.^2.*diag(Mstar);


%Accelerogram in Seism
%
Seism=Tolm.*rsm(1);
%time interval: dt
dt=1/100;
%
%nb of points: nb
nb=length(Seism);
t=0:dt:(nb-1)*dt;

%FIRST MODE
%Linear interpolation
%****************************************************************
%Computing of coefficients A,B,C,D,Ap,Bp,Cp,Dp
%****************************************************************
%
w1=w(1);
%preliminairy calculations
e=exp(-dumping*w1*dt);
wd=w1*sqrt(1-dumping^2);
wdt=wd*dt;
Rz=sqrt(1-dumping^2);
%
a=e*(dumping/Rz*sin(wdt)+cos(wdt));
B=e/wd*sin(wdt);
C=-1/w1^2*(2*dumping/w1/dt+e*(((1-2*dumping^2)/wdt-dumping/Rz)*sin(wdt)-(1+2*dumping/w1/dt)*cos(wdt)));
D=-1/w1^2*(1-2*dumping/(w1*dt)+e*((2*dumping^2-1)/wdt*sin(wdt)+2*dumping/w1/dt*cos(wdt)));
Ap=-e*(w1/Rz*sin(wdt));
Bp=e*(cos(wdt)-dumping/Rz*sin(wdt));
Cp=-1/w1^2*(-1/dt+e*((w1/Rz+dumping/dt/Rz)*sin(wdt)+1/dt*cos(wdt)));
Dp=-1/w1^2/dt*(1-e*(dumping/Rz*sin(wdt)+cos(wdt)));
%
%****************************************************************
%Computing of displacements and velocities
%****************************************************************
%
u(1)=0;
v(1)=0;
%
for I=1:nb-1
	u(I+1)=a*u(I)+B*v(I)+C*Seism(I)+D*Seism(I+1);
	v(I+1)=Ap*u(I)+Bp*v(I)+Cp*Seism(I)+Dp*Seism(I+1);
end
%
z1=u;
z2=z1*A(2,1)/A(1,1);
z3=z1*A(3,1)/A(1,1);
z4=z1*A(4,1)/A(1,1);
z=[z1;z2;z3;z4];
x=A*z;

%****************************************************************
%plotting the response
%****************************************************************

figure()
plot(t,z1,'DisplayName','First floor displacement')
grid on
legend('Location','southeast')
ylabel('z [mm]')
xlabel('t [s]')
title('First mode modal displacement of first floor')


figure()
plot(t,x(1,:),'DisplayName','First floor displacement')
hold on 
plot(t,x(2,:),'DisplayName','Second floor displacement')
plot(t,x(3,:),'DisplayName','Third floor displacement')
plot(t,x(4,:),'DisplayName','Fourth floor displacement')
ylabel('x [mm]')
xlabel('t [s]')
legend('Location','southeast')
title('Displacement of the frame due to Tolmezzo earthquake')
grid on

