%Notes:
%to have an idea of how the coordinate system is set up please refer to the
%lagrangian model
clear all
clc
close all
%% variables
test=1; %plots the test graphs in the code. useful for debugging
c=0.6/1000; %turns chrod length to m
n=20; %number of wing elements
%% loads and assigns data
digits(4); % sets decimal point accuracy
load('AnglesInter.mat') %loads previously generated data

phi=yy2; %stroke angle
psi=yy1; %deviation angle
beta=yy3; %rotation angle
time=xx/220;
%% filtering the position data due to noise by me
[b, a] = butter(4, 3.5/(250/2),'low');
phi_f=filtfilt(b, a, phi);
psi_f=filtfilt(b, a, psi);
beta_f=filtfilt(b, a, beta);
if test==1
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
%% angular velocities
phi_dot=diff(phi)/(time(2)-time(1));
psi_dot=diff(psi)/(time(2)-time(1));
beta_dot=diff(beta)/(time(2)-time(1));

phi_dotf=diff(phi_f)/(time(2)-time(1));
psi_dotf=diff(psi_f)/(time(2)-time(1));
beta_dotf=diff(beta_f)/(time(2)-time(1));
if test==1
    % these figures show a comparison between filtered and unfiltered data
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

%% Stationary wing reference frame
ex=[1;0;0];
ey=[0;1;0];
ez=[0;0;1];

%% rotations
syms beta1 phi1 psi1
Rx = [1 0 0; 0 cosd(beta1) -sind(beta1); 0 sind(beta1) cosd(beta1)];
Ry = [cosd(phi1) 0 sind(phi1); 0 1 0; -sind(phi1) 0 cosd(phi1)];
Rz = [cosd(psi1) -sind(psi1) 0; sind(psi1) cosd(psi1) 0; 0 0 1];

R=Rz*Rx*Ry; % complete rotation from the stationary wing base frame to the moving wing frame

R_inv=inv(R); %from moving frame to stationary frame
%% moving vectors in terms of stationary frame
ex1=R_inv*ex;
ey1=R_inv*ey;
ez1=R_inv*ez;
%% angular velocity in vector form
disp('Calculating the angular velocity in vector form')
for i=1:length(phi_dotf)
    beta1=beta_f(i);
    phi1=phi_f(i);
    psi1=psi_f(i);
    omega(1:3,i)=phi_dotf(i)*vpa(subs(ey1))+psi_dotf(i)*vpa(subs(ez1))+beta_dotf(i)*vpa(subs(ex1));
    i
end
disp('done with vector ang vel')
%% Magnitude of ang vel in deg/s
disp('calculating magnitude of angular velocity')
for i=1:length(omega)
    omega_mag(i)=norm(omega(1:3,i));
    i
end
disp('Finished calulating ang vel magnitude')
omega_rad=omega_mag*pi/180; % magnitude in radians

%% Find the location of center of pressure
alpha=0:0.01:pi; %angle used in this equation is in radians. I check by comparing to paper
x_test=(0.82*alpha/pi+0.05); %normalized with respect to chord length

x_cp=c*(0.82*abs(psi_f*pi/180)/pi+0.05); %location of center of pressure in x-axis
delz=0.002/10;
%% finds the distance vector to each center of pressure for a single element

for i=1:length(x_cp)
    r_cpp(1:3,i)=[-x_cp(i); 0; delz];
    beta1=beta_f(i);
    phi1=phi_f(i);
    psi1=psi_f(i);
    r_cp(1:3,i)=vpa(subs(R_inv*r_cpp(1:3,i)));
    i
end

%% test plot for the center of pressure
    for i=1:991
        plot3(r_cp(1,i),r_cp(2,i),r_cp(3,i),'*')
        hold on
    end
    title('location of point of pressure for a single element throughout a wingstroke in the moving wing frame')
   %%
    figure
    for i=1:length(r_cpp)
        plot3(r_cpp(1,i),r_cpp(2,i),r_cpp(3,i),'*')
        hold on
    end
    title('location of point of pressure for a single element throughout a wingstroke in the stationary wing frame')
    figure
    plot(x_cp)
    
%% Functions---------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



