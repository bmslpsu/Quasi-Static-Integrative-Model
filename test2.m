mean_f=0;
for i=1:20
    mean_f=mean_f+ mean(mean(element3(i).force_drag_vec(1,:))*element3(i).location_cop(3,:));
end
%%
Torque_element=[];
for j=1:length(element3)
   
    Torque_element(j)=mean(mean(element3(j).force_AM_vec(1,:)+element3(j).force_drag_vec(1,:)+...,
        element3(j).force_Rot_vec(1,:))*element3(j).location_cop(3,:));
    %Torque_element_2(j)=mean(mean(element3(j).force_AM_vec(3,:)+element3(j).force_drag_vec(3,:)+...,
      %  element3(j).force_Rot_vec(3,:))*element3(j).location_cop(1,:));
end
Torque_Mag_element=sum(Torque_element) %mean torque through each wingstroke for each element
%Torque_Mag_element=sum(Torque_element_2) 
