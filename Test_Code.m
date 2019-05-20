
%% This code was made to verify if there is a difference between using subscribed angles vs actual data
%based on the ang vel, i found out that both methods give us the same
%angular velocity. Based on this I can assume angular velocity values are
%not wrong
time_test=linspace(0,1,991)/220;
phi=-65*cos(2*3.14*220*time_test)+(65-47.5);
plot(time,phi)
hold on
plot(time,phi_f)
figure
phi_dot=diff(phi)/(time_test(2)-time_test(1));
plot(time_test(1:end-1),phi_dot)
hold on
plot(time(1:end-1),phi_dotf)
%% find the lift force based on this data
F_lift_test =1.22*0.5*(phi_dot*1.5/2000*3.14/180).^2*1.5/1000;
figure
plot(F_lift_test)
hold on
F_lift_test_2 =1.22*0.5*(phi_dotf*1.5/2000*3.14/180).^2*1.5/1000;
plot(F_lift_test_2)