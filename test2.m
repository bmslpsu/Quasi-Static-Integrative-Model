mean_f=0;
for i=1:20
   mean_f=mean_f+ mean(mean(element3(i).force_drag_vec(1,:))*element3(i).location_cop(3,:));
end