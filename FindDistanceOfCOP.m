function element=FindDistanceOfCOP(phi_f,psi_f,beta_f,n,c,R_inv)
%% Find the location of center of pressure
global wing_length
alpha=0:0.01:pi; %angle used in this equation is in radians. I check by comparing to paper
x_test=(0.82*alpha/pi+0.05); %normalized with respect to chord length
x_cp=c*(0.82*abs(psi_f*pi/180)/pi+0.05); %location of center of pressure in x-axis
delz=wing_length/n;

% finds the distance vector to each center of pressure for a each element
for j=1:n
    disp(['Finding COP for element' num2str(j)])
    for i=1:length(x_cp)-1
        r_cpp(1:3,i)=[-x_cp(i); 0; delz/2+delz*(j-1)];
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