
%% This code was made to verify if there is a difference between using subscribed angles vs actual data
%based on the ang vel, i found out that both methods give us the same
%angular velocity. Based on this I can assume angular velocity values are
%not wrong
time=(0:0.01:2)/220;
phi=-65*cos(3.14*220*time);
plot(time,phi)
figure
phi_dot=diff(phi)/(time(2)-time(1));
plot(time(1:end-1),phi_dot)
