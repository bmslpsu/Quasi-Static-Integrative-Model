clear all 
clc
digits(4)
%this function analyzes the QS data but uses a slightly different method to
%do so
%% get the intact wing torque
root='C:\Users\was29\Documents\MATLAB\QS_Data\';
[file_intactwing,path_intact] = uigetfile(root,'Select the data for the intact wing');
load([path_intact file_intactwing])
%% start torque calculations for intact wing
[torque_intact_total]=Find_Wing_Torque_2(element3);
S_intact=S_3;
%% load the chordwise damaged data
third_moment=[];
[file_chordwise,path_chordwise] = uigetfile(root,'Select the data for the chordwise damaged wing','MultiSelect', 'on');
for i=1:length(file_chordwise)
    load([path_intact file_chordwise{i}])
    [Total_Torque(i)]=Find_Wing_Torque_2(element3);
    third_moment=[third_moment S_3];
    i
end
 plot(third_moment/S_intact,(torque_intact_total-Total_Torque)/(10^-6*0.002*9.81))
 title('Torque vs (S3_D/S3_I) for chrodwise wing damage (reduced wing span length)')
 xlabel('S3 Damaged/ S3 Intact')
 ylabel('Normalized torque T/(mgl)')
 ylim([-0.05 0.15])
 
 %% load the spanwise damaged data
[file_spanwise,path_span] = uigetfile(root,'Select the data for the spanwise damaged wing','MultiSelect', 'on');
third_moment_2=[];
for i=1:length(file_spanwise)
    load([path_span file_spanwise{i}])
    [Total_Torque_2(i)]=Find_Wing_Torque_2(element3);
    third_moment_2=[third_moment_2 S_3];
    i
end
hold on
plot(third_moment_2/S_intact,(torque_intact_total-Total_Torque_2)/(10^-6*0.002*9.81))
title('Torque vs (S3_D/S3_I) for spanwise wing damage (reduced wing chord length)')
xlabel('S3 Damaged/ S3 Intact')
ylabel('Normalized torque T/(mgl)')
ylim([-0.05 0.15])
%% plot 
 plot(third_moment/S_intact,(torque_intact_total-Total_Torque)/(10^-6*0.002*9.81))
 hold on
plot(third_moment_2/S_intact,(torque_intact_total-Total_Torque_2)/(10^-6*0.002*9.81))
 title('Torque vs (S3_D/S3_I) for chrodwise and spanwise wing damage')
 xlabel('S3 Damaged/ S3 Intact')
 ylabel('Normalized torque T/(mgl)')
 ylim([-0.05 0.15])
 legend('Chordwise damage','Spanwise damage')
%-----------------------------------------------------------
%% Functions-----------------------------------------------
function [Torque_mag]=Find_Wing_Torque(element3)
Torque_element=[];
for j=1:length(element3)
    Torque_element(j)=mean(mean(element3(j).force_AM_vec(1,:)+element3(j).force_drag_vec(1,:)+...,
        element3(j).force_Rot_vec(1,:))*element3(j).location_cop(3,:));
   %Torque_element_2(j)=mean(mean(element3(j).force_AM_vec(3,:)+element3(j).force_drag_vec(3,:)+...,
    %    element3(j).force_Rot_vec(3,:))*element3(j).location_cop(1,:));
end
Torque_mag_1=sum(Torque_element); %mean torque through each wingstroke for each element
%Torque_Mag_2=sum(Torque_element_2); 
Torque_mag=Torque_mag_1; % +Torque_Mag_2;
end

function [Torque_mag]=Find_Wing_Torque_2(element3)
Torque_element=[];
for j=1:length(element3)
    Torque_element(j)=mean(mean(element3(j).force_AM_vec(1,:)+element3(j).force_drag_vec(1,:)+...,
        element3(j).force_Rot_vec(1,:))*element3(j).location_cop(3,:));
    Torque_element_2(j)=mean(mean(element3(j).force_AM_vec(3,:)+element3(j).force_drag_vec(3,:)+...,
        element3(j).force_Rot_vec(3,:))*element3(j).location_cop(1,:));
end
Torque_mag_1=sum(Torque_element); %mean torque through each wingstroke for each element
Torque_Mag_2=sum(Torque_element_2); 
Torque_mag=Torque_mag_1+Torque_Mag_2;
end