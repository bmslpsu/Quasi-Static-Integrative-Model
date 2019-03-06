%Notes:
%to have an idea of how the coordinate system is set up please refer to the
%lagrangian model
clear all
clc
close all
%% loads and assigns data
digits(4); % sets decimal point accuracy
load('AnglesInter.mat') %loads previously generated data
%% variables
test=0; %plots the test graphs in the code. useful for debugging
c=0.6/1000; %turns chrod length to m
n=10; %number of wing elements
global time
global wing_length
time=xx/220;
wing_length=2/1000;
%% Extracts wings angles from data
[phi_f, psi_f, beta_f,phi,psi,beta]=ExtractAngles(xx,yy1,yy2,yy3);
%%
[phi_dotf,psi_dotf,beta_dotf] =GetEulerAngleVelocity(phi_f, psi_f, beta_f, phi,psi,beta);

%% Stationary wing reference frame
ex=[1;0;0];
ey=[0;1;0];
ez=[0;0;1];
%% Euler Rotation
R_inv=EulerRotation();

%% moving vectors in terms of stationary frame
ex1=R_inv*ex;
ey1=R_inv*ey;
ez1=R_inv*ez;


%% find angular velocity of wing with respect to stationary frame
[omega, omega_mag,omega_rad]  =GetWingAngVel(ex1,ey1,ez1,phi,psi,beta,phi_dotf,psi_dotf,beta_dotf);
figure
plot(omega_mag)
title('Angular velocity of a wing with respect to a nonmoving frame at the base of the wing')

%% find the location of the center of pressure for each wing element
element=FindDistanceOfCOP(phi_f,psi_f,beta_f,n,c,R_inv);

%% find linear velocity of each element for each time step
element =FindLinearVelocity(element, omega)

%% find the forces acting on each blade element

%% Functions---------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function element=LiftAndDragForces(element)

j=1;
rho=1.25; %density of air
for i=1:length(element(j).linear_vel)
    C_L=0.225+1.58*sind(2.13*beta_f(i)-7.28);
    C_D=1.92-1.55*cosd(2.04*befa_f(i)-9.82);
end

end

function element=FindLinearVelocity(element, omega)
%finds the linear velocity of each element throughout a full wing stroke
j=1;
disp('calculating the linear velocity for each wing')
for i=1:length(omega)
    V_linear(1:3,i)=cross(omega(1:3,i)*180/pi,element(j).location_cop(1:3,i));
    i
end
element(j).linear_vel=V_linear;
disp('done calculating linear velocity')
end

function element=FindDistanceOfCOP(phi_f,psi_f,beta_f,n,c,R_inv)
%% Find the location of center of pressure
global wing_length
alpha=0:0.01:pi; %angle used in this equation is in radians. I check by comparing to paper
x_test=(0.82*alpha/pi+0.05); %normalized with respect to chord length
x_cp=c*(0.82*abs(psi_f*pi/180)/pi+0.05); %location of center of pressure in x-axis
delz=wing_length/n;

%% finds the distance vector to each center of pressure for a each element
for j=1:n
    disp('Finding COP for 1st element')
    for i=1:length(x_cp)-1
        r_cpp(1:3,i)=[-x_cp(i); 0; delz*j];
        beta1=beta_f(i);
        phi1=phi_f(i);
        psi1=psi_f(i);
        r_cp(1:3,i)=vpa(subs(R_inv*r_cpp(1:3,i)));
        i
    end
    disp('Complete for this element')
    element(j).location_cop=r_cp;
    element(j).locationInMovingFrame=r_cpp;
end
disp('Complete for entire wing')

%% test plot for the center of pressure
test=0;
if test==1
    for j=1:length(element)
        figure
        for i=1:length(element(j).location_cop)
            plot3(element(j).location_cop(1,i),element(j).location_cop(2,i),element(j).location_cop(3,i),'*')
            hold on
            title('location of point of pressure for a single element throughout a wingstroke in the statonary wing frame')
        end
    end
    
    figure
    %this plot below is not complete. will only plot the last element
    for i=1:length(r_cpp)
        plot3(r_cpp(1,i),r_cpp(2,i),r_cpp(3,i),'*')
        hold on
    end
    title('location of point of pressure for a single element throughout a wingstroke in the moving wing frame')
    figure
    plot(x_cp)
end
end

function [phi_f, psi_f, beta_f,phi,psi,beta]=ExtractAngles(xx,yy1,yy2,yy3)
phi=yy2; %stroke angle
psi=yy1; %deviation angle
beta=yy3; %rotation angle

%% filtering the position data due to noise by me
[b, a] = butter(4, 3.5/(250/2),'low');
phi_f=filtfilt(b, a, phi);
psi_f=filtfilt(b, a, psi);
beta_f=filtfilt(b, a, beta);

%% plots
% these figures show a comparison between filtered and unfiltered data
figure
plot(phi_f)
hold on
plot(phi)
figure
plot(psi_f)
hold on
plot(psi)
figure
plot(beta_f)
hold on
plot(beta)
end

function [phi_dotf,psi_dotf,beta_dotf] =GetEulerAngleVelocity(phi_f, psi_f, beta_f, phi,psi,beta)
%% angular velocities
global time
phi_dot=diff(phi)/(time(2)-time(1));
psi_dot=diff(psi)/(time(2)-time(1));
beta_dot=diff(beta)/(time(2)-time(1));

phi_dotf=diff(phi_f)/(time(2)-time(1));
psi_dotf=diff(psi_f)/(time(2)-time(1));
beta_dotf=diff(beta_f)/(time(2)-time(1));

figure
plot(phi_dotf)
hold on
plot(phi_dot)
figure
plot(psi_dotf)
hold on
plot(psi_dot)
figure
plot(beta_dotf)
hold on
plot(beta_dot)
end

function R_inv=EulerRotation()
%% rotations
syms beta1 phi1 psi1
Rx = [1 0 0; 0 cosd(beta1) -sind(beta1); 0 sind(beta1) cosd(beta1)];
Ry = [cosd(phi1) 0 sind(phi1); 0 1 0; -sind(phi1) 0 cosd(phi1)];
Rz = [cosd(psi1) -sind(psi1) 0; sind(psi1) cosd(psi1) 0; 0 0 1];
R=Rz*Rx*Ry; % complete rotation from the stationary wing base frame to the moving wing frame
R_inv=inv(R); %from moving frame to stationary frame
end

function [omega, omega_mag,omega_rad]  =GetWingAngVel(ex1,ey1,ez1,phi,psi,beta,phi_dotf,psi_dotf,beta_dotf)
%% angular velocity in vector form
disp('Calculating the angular velocity in vector form')
for i=1:length(phi_dotf)
    beta1=beta(i);
    phi1=phi(i);
    psi1=psi(i);
    omega(1:3,i)=phi_dotf(i)*vpa(subs(ey1))+psi_dotf(i)*vpa(subs(ez1))+beta_dotf(i)*vpa(subs(ex1));
    i
end
disp('done with vector ang vel')
%% Magnitude of ang vel in deg/s
disp('calculating magnitude of angular velocity')
for i=1:length(omega)
    omega_mag(i)=norm(omega(1:3,i));
end
disp('Finished calulating ang vel magnitude')
omega_rad=omega_mag*pi/180; % magnitude in radians
end