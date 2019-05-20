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

    omega(1:3,i)=phi_dotf(i)*vpa(subs(ey1))+psi_dotf*vpa(subs(ez1))+beta_dotf*vpa(subs(ex1));
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