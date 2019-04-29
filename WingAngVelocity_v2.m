%Notes:
%to have an idea of how the coordinate system is set up please refer to the
%lagrangian model
%% note on the axese
%z-axis is along the length of the wing
%x-axis is along the chord of the wing starting at the wing base and
%is parallel to the abdomen of the fly
%y-axis is perpendicular to the surface of the wing
clear all
clc
close all
%% loads and assigns data
digits(4); % sets decimal point accuracy
load('AnglesInter.mat') %loads previously generated data
%% variables
test=0; %plots the test graphs in the code. useful for debugging


n=10; %number of wing elements
global time
global wing_length rho
rho=1.255;
global c 
global C_r 
C_r=1.55;
c=0.6/1000; %turns chrod length to m
time=xx/220;
wing_length=2/1000; % winglength in meters
%% Extracts wings angles from data
[phi_f, psi_f, beta_f,phi,psi,beta]=ExtractAngles(xx,yy1,yy2,yy3);
%beta 
%psi 
%phi
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
%omega_mag is in deg
%omega is in deg
[omega, omega_mag,omega_rad]  =GetWingAngVel(ex1,ey1,ez1,phi,psi,beta,phi_dotf,psi_dotf,beta_dotf);
figure
plot(omega_mag)
title('Angular velocity of a wing with respect to a nonmoving frame at the base of the wing')

%% find the location of the center of pressure for each wing element
element=FindDistanceOfCOP(phi_f,psi_f,beta_f,n,c,R_inv);

%% find linear velocity of each element for each time step
element =FindLinearVelocity(element, omega);

%% find the lift and drag forces acting on each blade element
del_r=wing_length/n;
element =LiftAndDragForces(element,phi_f,del_r);

%%  Find the added mass force acting on each wing
% note: the input anglur velocity is in deg/s. it is converted to rad/s in
% the function. 
element1=FindRotationalForce(element,beta_dotf,del_r,c);

%% Added mass force
element2=FindAddedMass(element1,beta_dotf,beta_f,del_r,c,time);

%% find force directions
element3=FindForceVectors(element2,R_inv,ey1,ex1,beta_f,phi_f,psi_f);

%% find the total forces in x y and z directions
[force_x, force_y, force_z]=Find_forces_XYZ(element3);
%% Functions---------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [force_x, force_y, force_z]=Find_forces_XYZ(element)

force_Total=0; % the total force in x y and z directions
for j=1:length(element)
   force_Total=element(j).force_Rot_vec + element(j).force_AM_vec + element(j).force_lift_vec + element(j).force_drag_vec+force_Total;
end
   force_x=force_Total(1,:);
   force_y=force_Total(2,:);
   force_z=force_Total(3,:);

end

function element=FindForceVectors(element,R_inv,ey1,ex1,beta_f,phi_f,psi_f)
%lift and drag are assumed vertical and horizontal
%rotation and added mass forces are perpendicular to wing surface
syms psi1 beta1 phi1

%the normal to the surface of the wing defined in the wing reference frame
%will be in the xy plane of the system
e_WingNormal=[cosd(90-psi1); sind(90-psi1); 0]; %define the vector of the added mass and rot force
tic
for j=1:length(element)
    disp(['calculating vector force for element ' num2str(j)])
    for i=1:length(beta_f)-2 %i had to use -2 here because the addedmass force has an acceleration component
        %and using the diff function reduces the length of the vector by 1
        beta1=beta_f(i);
        phi1=phi_f(i);
        psi1=psi_f(i);
        Rot_matrix=vpa(subs(R_inv));
        f_lift_vec(:,i)=Rot_matrix*element(j).force_Lift(i)*ey1;
        f_drag_vec(:,i)=Rot_matrix*element(j).force_Drag(i)*ex1;
        f_addedMass_vec(:,i)=Rot_matrix*element(j).force_AddedMass(i)*e_WingNormal;
        f_Rot_vec(:,i)=Rot_matrix*element(j).force_Rotation(i)*e_WingNormal;

        if i==floor(length(beta_f)/4) 
            disp(['Calculations for element ' num2str(j) ' are 1/4 done'])
        elseif i==floor(length(beta_f)/2)
            disp(['Calculations for element ' num2str(j) ' are 1/2 done'])
        end
        
    end
    element(j).force_lift_vec=vpa(subs(f_lift_vec));
    element(j).force_drag_vec=vpa(subs(f_drag_vec));
    element(j).force_AM_vec=vpa(subs(f_addedMass_vec));
    element(j).force_Rot_vec=vpa(subs(f_Rot_vec));
end
toc
end

function element=FindAddedMass(element,beta_dotf,beta_f,del_r,c,time)
global rho
for j=1:length(element)
    disp(['Calcualting added mass force for element ' num2str(j)])
    linear_acc=diff(eval(element(j).linear_vel)')'/(time(2)-time(1));
    % due to present noise, I have to filter the acceleration as well
    [b, a] = butter(4, 20.5/(250/2),'low');
    linear_acc_filt1=filtfilt(b, a, linear_acc(1,:));
    linear_acc_filt2=filtfilt(b, a, linear_acc(2,:));
    linear_acc_filt3=filtfilt(b, a, linear_acc(3,:));
    figure
    plot(linear_acc(1,:))
    hold on
    plot(linear_acc_filt1)
    title('Linear acceleration with filtered data')
    linear_acc_filt=[linear_acc_filt1; linear_acc_filt2;linear_acc_filt3];
    %issue with lengths: when using the derivative, the length of the new array
    %is one element smaller. In order to continue with the analysis, the larger
    %elements must be cropped to match the same size
    for i=1:length(linear_acc_filt1)
        part1=rho*pi*c^2/4*del_r;
        part2=(dot(element(j).linear_vel(:,i),linear_acc_filt(:,i))*sind(beta_f(i)))/element(j).linear_vel_norm(i);
        part3=eval(element(j).linear_vel_norm(i)*beta_dotf(i)*cosd(beta_f(i)));
        f_addedMass(i)=part1*(part2+part3);
    end
    disp(['Done calculating added mass force for element' num2str(j)])
    beep
    element(j).force_AddedMass=f_addedMass;
end
end

function element =FindRotationalForce(element,beta_dotf,del_r,c)
%del_r (m)
%c (m)
%C_r from lit (1.55)
%beta_dotf (deg/s) converted to rad/s in the function for force
%element: the struct that has the information for each element in the wing
%throughout a wing stroke
global C_r
global rho
for j=1:length(element)
    disp(['Calculating rotation force for element' num2str(j)]);
    for i=1:length(element(j).linear_vel_norm)
        F_rot(i)=eval(C_r*rho*c^2*del_r*element(j).linear_vel_norm(i)*(beta_dotf(i)*pi/180)^2);
        
        if i==floor(length(element(j).linear_vel_norm)/2)
            disp('calculation of rotational force is half done')
        end
    end
    element(j).force_Rotation=F_rot;
end
end

function element=LiftAndDragForces(element,phi_f,del_r)
%this function will find only the magnitude not the direction

for j=1:length(element)
global rho

disp('calculating force for one element')
for i=1:length(element(j).linear_vel)
    C_L=0.225+1.58*sind(2.13*phi_f(i)-7.28);
    C_D=1.92-1.55*cosd(2.04*phi_f(i)-9.82);
    element(j).force_Drag(i)=0.5*rho*norm(element(j).linear_vel(:,i))^2*C_D*del_r;
    element(j).force_Lift(i)=0.5*rho*norm(element(j).linear_vel(:,i))^2*C_L*del_r;
end
disp(['Finished force for element' num2str(j)])
end
disp('Done for entire wing')
end

function element=FindLinearVelocity(element, omega)
%finds the linear velocity of each element throughout a full wing stroke
for j=1:length(element)
disp(['calculating the linear velocity for element' num2str(j)])
for i=1:length(omega)
    V_linear(1:3,i)=cross(omega(1:3,i)*pi/180,element(j).location_cop(1:3,i));
    V_linear_Norm(i)=norm(V_linear(1:3,i));
end
element(j).linear_vel=V_linear;
element(j).linear_vel_norm=V_linear_Norm;
disp(['done calculating linear velocity for element ' num2str(j)])
end
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
    disp(['Finding COP for element' num2str(j)])
    for i=1:length(x_cp)-1
        r_cpp(1:3,i)=[-x_cp(i); 0; delz*j];
        beta1=beta_f(i);
        phi1=phi_f(i);
        psi1=psi_f(i);
        r_cp(1:3,i)=vpa(subs(R_inv*r_cpp(1:3,i)));
        if i==floor(length(x_cp)/2)
            disp(['COP distance calculation for element ' num2str(j) ' is half done'])
        end
    end
    disp(['Complete COP for element ' num2str(j) 'out of ' num2str(n)])
    element(j).location_cop=r_cp;
    element(j).locationInMovingFrame=r_cpp;
    element(j).Distance_COP=norm(r_cp);
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
%this function finds the rotation matrix from the fly stationary wing axis
%to the moving wing axis. 
%R is the rotation from the stationary to the moving. Therefore finding its
%inverse is required

% important note: as the wing is always moving, we cannot have one rotation
% matrix. to avoid having many matricies, i used a symbolic representation
% for this rotation. R_ivn; the output is a symbolic matrix which is
% evaluated for every datapoint of the euler angles when needed. however,
% it is not saved
syms beta1 phi1 psi1
Rx = [1 0 0; 0 cosd(beta1) -sind(beta1); 0 sind(beta1) cosd(beta1)];
Ry = [cosd(phi1) 0 sind(phi1); 0 1 0; -sind(phi1) 0 cosd(phi1)];
Rz = [cosd(psi1) -sind(psi1) 0; sind(psi1) cosd(psi1) 0; 0 0 1];
R=Rz*Rx*Ry; % complete rotation from the stationary wing base frame to the moving wing frame
R_inv=inv(R); %from moving frame to stationary frame
end

function [omega, omega_mag,omega_rad]  =GetWingAngVel(ex1,ey1,ez1,phi,psi,beta,phi_dotf,psi_dotf,beta_dotf)
%% angular velocity in vector form
% omega (deg)
% omega_mag(deg)
%omega_rad (rad)
disp('Calculating the angular velocity in vector form')
for i=1:length(phi_dotf)
    beta1=beta(i);
    phi1=phi(i);
    psi1=psi(i);
    omega(1:3,i)=phi_dotf(i)*vpa(subs(ey1))+psi_dotf(i)*vpa(subs(ez1))+beta_dotf(i)*vpa(subs(ex1));
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