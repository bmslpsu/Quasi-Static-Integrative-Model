function [ex1, ey1, ez1]=Find_vectors(R_inv,phi_f,psi_f,beta_f)
%this function finds the vectors in hopes of increasing the speed of
%the code
ex=[1;0;0];
ey=[0;1;0];
ez=[0;0;1];
for i=1:length(phi_f)
    beta1=beta_f(i);
    psi1=psi_f(i);
    phi1=phi_f(i);
    ex1(1:3,i)=vpa(subs(R_inv*ex));
    ey1(1:3,i)=vpa(subs(R_inv*ey));
    ez1(1:3,i)=vpa(subs(R_inv*ez));
end