%Notes:
%to have an idea of how the coordinate system is set up please refer to the
%lagrangian model

% v4 has vector calculations in symbolic at the start of the code to reduce
% calculation time for the rest of the code. this worked and the code is
% much faster
%% note on the axese
%z-axis is along the length of the wing
%x-axis is along the chord of the wing starting at the wing base and
%is parallel to the abdomen of the fly
%y-axis is perpendicular to the surface of the wing
clear all
clc
close all
tic
%% loads and assigns data
digits(4); % sets decimal point accuracy
load('AnglesInter.mat') %loads previously generated data
%% variables
test=0; %plots the test graphs in the code. useful for debugging
n=25; %number of wing elements
global time
global wing_length rho
rho=1.255;
global c
global C_r
C_r=1.55;
c=0.59/1000; %turns chrod length to m
time=(xx-xx(1))/220;
wing_length=2/1000; % winglength in meters
del_r=wing_length/n; % the length of each element along the span
%% Extracts wings angles from data
[phi_f, psi_f, beta_f,phi,psi,beta]=ExtractAngles(xx,yy1,yy2,yy3);
%% find the angular velocity of the wing for each euler angle
[phi_dotf,psi_dotf,beta_dotf] =GetEulerAngleVelocity(phi_f, psi_f, beta_f, phi,psi,beta);
%% Stationary wing reference frame
ex=[1;0;0];
ey=[0;1;0];
ez=[0;0;1];

%% Euler Rotation
R_inv=EulerRotation();

%% finding the location of the wing center of mass wrt to the body the fly
% this code  was added here as part of the gant jean is writing. It plays
% no part in the QS code and the force analysis.
syms a_1 b_1 % the coordinates of the wing center of mass in the wing reference frame
% note as the wing ref frame moves with the wing, the location of the wing
% COM remains constant in that frame

COM_cord=[a_1;0;b_1]; % since the wing is in the xz plane of the ref frame

COM_Ground=R_inv*COM_cord; % takes us from the wing ref to the wing base ref which is stationary wrt to the fly
R_z=[cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1];
COM_Ground_2=R_z*COM_Ground;
%gravity acts along the y-axis of the global frame
COM_height=COM_Ground_2(2);

%% moving vectors in terms of stationary frame (no longer needed)
ex1=R_inv*ex;
ey1=R_inv*ey;
ez1=R_inv*ez;
%% moving vectors in stationary frame 2: hope to reduce runtime in later code
[ex11, ey11, ez11,R_inv2]=Find_vectors(R_inv,phi_f,psi_f,beta_f);

%% find angular velocity of wing with respect to stationary frame
%omega_mag is in deg
%omega is in deg
[omega, omega_mag,omega_rad]=GetWingAngVel(ex11,ey11,ez11,phi_dotf,psi_dotf,beta_dotf);
figure
plot(omega_mag)
%'note: this gives the absolute value of the magnitude not the sign'
title('Angular velocity of a wing with respect to a nonmoving frame at the base of the wing')

%% find the location of the center of pressure for each wing element
element=FindDistanceOfCOP(psi_f,n,c,R_inv2);

%% find linear velocity of each element for each time step
element =FindLinearVelocity(element, omega ); %omega instead of phi_f

%% find the lift and drag forces acting on each blade element
element =LiftAndDragForces(element,phi_f,del_r,1);
%%  Find the added mass force acting on each wing
% note: the input anglur velocity is in deg/s. it is converted to rad/s in
% the function.
element1=FindRotationalForce(element,beta_dotf,del_r,c);

%% Added mass force
element2=FindAddedMass(element1,beta_dotf,beta_f,del_r,c,time);

%% find force directions
element3=FindForceVectors(element2,R_inv2,ey11,ex11,ez11,beta_f,phi_f,psi_f);

%% find the total forces in x y and z directions
[force_x, force_y, force_z]=Find_forces_XYZ(element3);

%% plots used to analyze the data
Create_Plots(phi_f,force_y,force_x,time)

%% find the third moment of inertia
S_3=Find_Third_Moment(wing_length,c);
%% save script for the data
save_name=['WingRemaining_' num2str(wing_length/0.002*100) '_' num2str(c/0.0006*100) '_Percent.mat'];
save_root='C:\Users\was29\Documents\MATLAB\QS_Data\';
save([save_root save_name],'element3','force_x','force_y','force_z','c','wing_length','time','S_3')
%% end of code timer
toc
%% Functions---------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function S_3=Find_Third_Moment(wing_length,c)
n_elements=500;
x_cor=linspace(wing_length/(2*n_elements),wing_length-wing_length/(2*n_elements),n_elements);
y_cor=linspace(c/(2*n_elements),c-c/(2*n_elements),n_elements);
area_el=(c/n_elements)*wing_length/n_elements;
S_3=0;
for j=1:length(y_cor)
    for i=1:length(x_cor)
        S_3=S_3+x_cor(i)^3*(area_el);
    end
end

end
function []=Create_Plots(phi_f,force_y,force_x,time)
figure
subplot(3,1,1)
plot(time,phi_f)
title('Stroke angle')
subplot(3,1,2)
plot(time(1:end-2),force_y*2)
title('Upwards force')
subplot(3,1,3)
plot(time(1:end-2),force_x*2)
title('Force in x-direction (drag+addded mass +rot)')
end

function [force_x, force_y, force_z]=Find_forces_XYZ(element)

force_Total=0; % the total force in x y and z directions for each wing
for j=1:length(element)
    force_Total=element(j).force_Rot_vec + element(j).force_AM_vec + element(j).force_lift_vec + element(j).force_drag_vec+force_Total;
end
force_x=force_Total(1,:);
force_y=force_Total(2,:);
force_z=force_Total(3,:);

end

function element=FindForceVectors(element,R_inv2,ey11,ex11,ez1,beta_f,phi_f,psi_f)
%lift and drag are assumed vertical and horizontal
%rotation and added mass forces are perpendicular to wing surface
syms psi1 beta1 phi1

%the normal to the surface of the wing defined in the wing reference frame
%will be in the xy plane of the system
%e_WingNormal=[cosd(90-psi1); sind(90-psi1); 0]; %define the vector of the added mass and rot force
for j=1:length(element)
    disp(['calculating vector force for element ' num2str(j)])
    for i=1:length(beta_f)-2 %i had to use -2 here because the addedmass force has an acceleration component
        %and using the diff function reduces the length of the vector by 1
        e_WingNormal=[cosd(90-psi_f(i)); sind(90-psi_f(i)); 0]; %define the vector of the added mass and rot force
        Rot_matrix=R_inv2(:,:,i);
        % f_lift_vec(:,i)=Rot_matrix*element(j).force_Lift(i)*ey11(1:3,i);
        % f_drag_vec(:,i)=Rot_matrix*element(j).force_Drag(i)*ex11(1:3,i);
        f_lift_vec(:,i)=element(j).force_Lift(i)*[0;1;0]; %removed the rotation matrix because the direction of these two forces is known 
        f_drag_vec(:,i)=element(j).force_Drag(i)*[1;0;0];
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
    disp(['Calculations for element ' num2str(j) ' are done'])
    
end
end

function element=FindAddedMass(element,beta_dotf,beta_f,del_r,c,time)
global rho
for j=1:length(element)
    disp(['Calcualting added mass force for element ' num2str(j)])
    linear_acc=diff((element(j).linear_vel)')'/(time(2)-time(1));
    % due to present noise, I have to filter the acceleration as well
    [b, a] = butter(4, 0.1,'low');
    linear_acc_filt1=filtfilt(b, a, linear_acc(1,:));
    linear_acc_filt2=filtfilt(b, a, linear_acc(2,:));
    linear_acc_filt3=filtfilt(b, a, linear_acc(3,:));
    linear_acc_filt=[linear_acc_filt1; linear_acc_filt2;linear_acc_filt3];
    %issue with lengths: when using the derivative, the length of the new array
    %is one element smaller. In order to continue with the analysis, the larger
    %elements must be cropped to match the same size
    for i=1:length(linear_acc_filt1)
        part1=rho*pi*c^2/4*del_r;
        part2=(dot(element(j).linear_vel(:,i),linear_acc_filt(:,i))*sind(beta_f(i)))/element(j).linear_vel_norm(i);
        part3=(element(j).linear_vel_norm(i)*beta_dotf(i)*cosd(beta_f(i)));
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
F_rot=zeros(1,length(element(1).linear_vel_norm));
for j=1:length(element)
    disp(['Calculating rotation force for element' num2str(j)]);
    for i=1:length(element(j).linear_vel_norm)
        F_rot(i)=C_r*rho*c^2*del_r*element(j).linear_vel_norm(i)*(beta_dotf(i)*pi/180); %modified this. should reduce force
        
        if i==floor(length(element(j).linear_vel_norm)/2) % if statement to keep track of where the code is
            disp(['calculation of rotational force for element' num2str(j) 'are half done'])
        end
    end
    element(j).force_Rotation=F_rot;
    disp(['calculation of rotational force for element' num2str(j) 'are done'])
end
end

function element=LiftAndDragForces(element,phi_f,del_r,test)
%this function will find only the magnitude not the direction
%test is used to plot some debugging graphs
for j=1:length(element)
    global rho c
    
    disp('calculating force for one element')
    for i=1:length(element(j).linear_vel)
        C_L=0.225+1.58*sind(2.13*phi_f(i)-7.28);
        C_D=1.92-1.55*cosd(2.04*phi_f(i)-9.82);
        element(j).force_Drag(i)=0.5*c*rho*norm(element(j).linear_vel(:,i))^2*C_D*del_r;
        element(j).force_Lift(i)=0.5*c*rho*norm(element(j).linear_vel(:,i))^2*C_L*del_r;
    end
    disp(['Finished force for element' num2str(j)])
end
disp('Done for entire wing')

if test==1
    Lift_wing_force=zeros(1,length(element(1).force_Lift));
    for j=1:length(element(1).force_Lift)
        Lift_wing_dummy=0;
        for i=1:length(element)
            Lift_wing_dummy=Lift_wing_dummy+element(i).force_Lift(j);
        end
        Lift_wing_force(j)=Lift_wing_dummy;
    end
end
figure 
plot(Lift_wing_force)
title('lift force for the whole wing throughout a stroke')
        
end


function element=FindLinearVelocity(element, omega)
%finds the linear velocity of each element throughout a full wing stroke
V_linear=zeros(3,length(omega));
V_linear_Norm=zeros(1,length(omega));
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

function [phi_f, psi_f, beta_f,phi,psi,beta]=ExtractAngles(xx,yy1,yy2,yy3)
phi=yy2; %stroke angle
psi=yy1; %deviation angle
beta=yy3; %rotation angle

%% filtering the position data due to noise by me
[b, a] = butter(2, 0.02,'low');
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
%% filtering also to be done on the derivative of the angle
[b, a] = butter(2, 0.02,'low');

phi_dotf=diff(phi_f)/(time(2)-time(1));
psi_dotf=diff(psi_f)/(time(2)-time(1));
beta_dotf=diff(beta_f)/(time(2)-time(1));
phi_dotf=filtfilt(b, a, phi_dotf);
psi_dotf=filtfilt(b, a, psi_dotf);
beta_dotf=filtfilt(b, a, beta_dotf);
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
R=Rz*Ry*Rx; % complete rotation from the stationary wing base frame to the moving wing frame
%R_inv=inv(R); %from moving frame to stationary frame
R_inv=inv(R);
end